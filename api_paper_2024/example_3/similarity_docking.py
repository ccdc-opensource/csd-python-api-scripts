#!/usr/bin/env python
#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2024-05-02: created by the Cambridge Crystallographic Data Centre
#

"""
    similarity_docking.py - search the CSD for similar structures to a known
    ligand, dock them, saving structures which score better than the known ligand.

    This has been adapted from the example program provided in the CSD Python API to allow for
    use of other data sources.
    
    Original program written by Richard Sykes
    Modified for use with Chembl by Jason C. Cole and Alex Moldovan
"""
######################################################################

import argparse
import math
import os
import sys

from ccdc import conformer
from ccdc.docking import Docker
from ccdc.entry import Entry
from ccdc.io import EntryReader, EntryWriter
from ccdc.molecule import Molecule
from ccdc.search import SimilaritySearch
from requests.exceptions import RequestException


def import_chembl_dependencies():
    try:
        from chembl_webresource_client.new_client import new_client
        from chembl_webresource_client.http_errors import HttpApplicationError
        from requests.exceptions import ConnectionError
    except ImportError:
        return False, \
            "To use the CHEMBl option you must have chembl_webresource_client installed. \n" \
            + "you can install it by running \n\npip install chembl_webresource_client\n\n" \
            + "See https://github.com/chembl/chembl_webresource_client for more information"
    return True, ""


class ChemblParseError(RuntimeError):
    def __init__(self, smiles, identifier, what):
        super().__init__(f"Molecule generation failed for {smiles} id {identifier}- {what}")


class ChemblSearchError(RuntimeError):
    def __init__(self, identifier, search_type, what):
        super().__init__(f"Chembl {search_type} search failed for id {identifier}- {what}")


class ChemblHit:
    def __init__(self, smiles, identifier):

        self.identifier = identifier

        try:
            m = Molecule.from_string(smiles)

            conformer_generator = conformer.ConformerGenerator()
            conformer_generator.settings.max_conformers = 1

            conformers = conformer_generator.generate(m)

        except (TypeError, RuntimeError) as e:
            raise ChemblParseError(smiles, identifier, f"{e}")

        if conformers is None or len(conformers) == 0:
            raise ChemblParseError(smiles, identifier, f"No conformers generated ")

        self.molecule = conformers[0].molecule.copy()
        self.molecule.identifier = identifier
        if len(self.molecule.atoms) == 0:
            raise ChemblParseError(smiles, identifier, "Molecule generated contains zero atoms")
        for atom in self.molecule.atoms:
            if atom.coordinates is None:
                raise ChemblParseError(smiles, identifier, f"Atom {atom.label} has no coordinate")


def id_to_chembl_id(identifier):
    return identifier.split('|')[0]


########################################################################################


class SimilarityDocking(argparse.ArgumentParser):
    """Defines arguments and runs the job."""

    def __init__(self):
        """Defines the arguments."""
        super(self.__class__, self).__init__(description=__doc__)

        # Input and output options
        self.add_argument(
            'conf_file_name',
            help='The configuration file to use.'
        )
        self.add_argument(
            '-l', '--ligand', default=None,
            help='Ligand to use in the docking. This could be a filename, CSD code, ChemblID or a SMILES code.'
        )

        self.add_argument(
            '-o', '--output-directory', default=os.path.abspath('output'),
            help='Directory in which results will be stored.'
        )

        # Similarity search parameters
        self.add_argument(
            '-m', '--max_hit_structures', default=None, type=int,
            help='The maximum number of CSD structures to try and dock (default: unlimited).'
        )

        # Docking parameters
        self.add_argument(
            '-f', '--fitness-function', default=None,
            help='The fitness function to use.'
        )
        self.add_argument(
            '-n', '--ndocks', default=10, type=int,
            help='The number of docking attempts to make.'
        )
        self.add_argument(
            '-s', '--speed', default=None,
            choices=('very fast', 'fast', 'medium', 'slow', 'very slow'),
            help=('How thorough the docking attempts should be. '
                  'Values correspond to 10%%, 25%%, 50%%, 75%% and 100%% '
                  'auto-scaling of docking parameters.')
        )
        self.add_argument(
            '-cmb', '--search-chembl', action='store_true',
            help='Whether to search chembl as the similarity data source rather than the CSD'
        )
        self.add_argument(
            '-mc', '--maximum-cycles', default=100, type=int,
            help='Maximum number of steps away from the original ligand to try'
        )
        self.add_argument(
            '-ms', '--minimum-similarity', default=50, type=int, choices=range(40, 91, 10),
            metavar="[40-90]",
            help='Mininum similarity to allow (%% value)'
        )
        self.add_argument(
            '-mw', '--maximum-weight', default=None, type=float,
            help='Maximum molecular weight of any new ligand to try'
        )
        self.add_argument(
            '-sp', '--start-point-threshold', default=0.9, type=float,
            help='Score threshold to include a new start point: If a new docking scores within '
                 '[start-point-threshold] * [score] of the best so far, it will be a new start point - default 0.9'
        )
        self.add_argument(
            '-msp', '--maximum-start-points', default=None, type=int,
            help='Maximum number of start points to retain. By default, all will be retained')

        self.best_poses = []

    def get_arguments(self):
        """Extract and process the arguments."""
        self.args = self.parse_args()
        if self.args.search_chembl:
            status, msg = import_chembl_dependencies()
            if status is False:
                self.error(msg)

            print("Searching CHEMBl via RESTful API")
        else:
            print("Searching in the CSD")

        self.settings = Docker.Settings.from_file(self.args.conf_file_name)
        if self.args.output_directory:
            if os.path.exists(self.args.output_directory):
                self.error(f"The output directory {self.args.output_directory} already exists - wont overwrite")

            self.settings.output_directory = self.args.output_directory

        if self.args.fitness_function is not None:
            self.settings.fitness_function = self.args.fitness_function
        if self.args.speed is not None:
            self.settings.autoscale = {
                'very fast': 10,
                'fast': 25,
                'medium': 50,
                'slow': 75,
                'very slow': 100
            }[self.args.speed]
        self.ligand_info = self.settings.ligand_files
        if self.args.ligand is None:
            self.args.ligand = self.ligand_info[0].file_name

        self.settings.clear_ligand_files()
        self.settings.set_hostname(ndocks=self.args.ndocks)
        self.docker = Docker(settings=self.settings)
        if not os.path.isdir(self.settings.output_directory):
            os.makedirs(self.settings.output_directory)
        self.winner_writer = EntryWriter(
            os.path.join(self.settings.output_directory, 'winners.mol2')
        )
        self.better_writer = EntryWriter(
            os.path.join(self.settings.output_directory, 'better.mol2')
        )
        self.fail_writer = EntryWriter(
            os.path.join(self.settings.output_directory, 'fails.mol2')
        )

        self.start_point_writer = EntryWriter(
            os.path.join(self.settings.output_directory, 'start_points.mol2')
        )

    def run(self):
        """Runs the job."""
        self.get_arguments()
        self.session = self.docker.dock(
            os.path.join(self.settings.output_directory, 'settings.conf'),
            'interactive'
        )
        self.run_substrate()
        self.done = set()
        dock_index = 1
        cycles = 0
        while cycles < self.args.maximum_cycles:
            print(f"Start cycle {cycles}")
            hits = self.search()
            if not hits:
                break

            cycles += 1
            for h in hits:
                self.done.add(h.identifier)
                mol = h.molecule

                if mol.all_atoms_have_sites is False:
                    print('%-8s: No 3d structure' % h.identifier)
                    self.fail_writer.write(mol)
                    continue
                mol.assign_bond_types()
                mol.add_hydrogens()

                # Skip this one as its beyond the weight threshold
                if self.args.maximum_weight and mol.molecular_weight > self.args.maximum_weight:
                    continue

                entry = Entry.from_molecule(mol)
                poses = self.session.dock(entry)
                dock_index += 1
                if not poses:
                    print('%-8s: Failed to dock' % h.identifier)
                    continue
                score, pose = self.best(poses)
                self.write_poses(poses=poses, dock_index=dock_index)
                self.write_scores(score=score, hit=h, pose=pose)

    def search(self):
        # If we have too many start points, chop off the lowest scoring ones
        self.best_poses = sorted(self.best_poses, key=lambda x: x[0])  # Only sort on
        if self.args.maximum_start_points is not None and len(self.best_poses) > self.args.maximum_start_points:
            self.best_poses = self.best_poses[len(self.best_poses) - self.args.maximum_start_points:]

        if self.args.search_chembl:
            return self.search_chembl()
        else:
            return self.search_csd()

    def search_chembl(self):
        try:
            return self._search_chembl()
        except ChemblSearchError as e:
            print(f"Chembl Searching failed with an exception {e} - script will terminate")
            return None

    def search_csd(self):
        try:
            return self._search_csd()
        except RuntimeError as e:
            print(f"CSD Searching failed with an exception {e} - script will terminate")
            return None

    def get_chembl_molecule(self, chemblid):
        """ Get a SMILEs entry from its ID

        :param chemblid:
        :return:
        """

        molecule = new_client.molecule
        try:
            m1 = molecule.filter(chembl_id=f'{chemblid}').only(['molecule_chembl_id', 'molecule_structures'])
        except (HttpApplicationError, ConnectionError, RequestException, RuntimeError, IOError) as e:
            raise ChemblSearchError(chemblid, "identifier", e)

        smiles = m1[0]['molecule_structures']['canonical_smiles']
        identifier = m1[0]['molecule_chembl_id']

        return ChemblHit(smiles, identifier)

    def _search_chembl(self):
        """
        Similarity search ChemBL. Will search with sequentially more relaxes settings to fins
        new as-yet-untried ligands up to a similarity threshold specified on the command line.
        :return: found hits
        """

        min_threshold = self.args.minimum_similarity

        hits = []
        similarity = new_client.similarity

        while not hits and len(self.best_poses) > 0:
            mol_score = self.best_poses[-1][0]
            mol = self.best_poses[-1][1].molecule
            mol.remove_unknown_atoms()
            mol.add_hydrogens()
            # mol.standardise_aromatic_bonds()
            threshold = 100.0

            while not hits and round(threshold, 6) > min_threshold:
                threshold -= 10.0
                print(f"Running {mol.identifier}, score {mol_score}, similarity threshold {threshold}")
                try:
                    # This minimises interaction with RDKit and SMILES - basically RDKit is strict about aromaticity
                    # so unfortunately fails, for example if a code contains carboxylates specified with aromatic bonds,
                    #  which is a standard used in mol2 files ...
                    if mol.identifier.startswith("CHEMBL"):
                        chemblid = id_to_chembl_id(mol.identifier)
                        print(f"Searching Chembl with identifier {chemblid}")
                        res = similarity.filter(chembl_id=f"{chemblid}", similarity=threshold).only(
                            ['molecule_structures', 'molecule_chembl_id', 'similarity'])
                    else:
                        smiles = mol.smiles
                        if smiles is None:
                            print('Cant get a good SMILES for entry')
                            break
                        print(f"Searching Chembl with SMILES {smiles}")
                        res = similarity.filter(smiles=f"{smiles}", similarity=threshold).only(
                            ['molecule_structures', 'molecule_chembl_id', 'similarity'])
                    print('Found %s structures similar to %s with threshold %.2f.' % (
                        len(res), mol.identifier, threshold))
                except (HttpApplicationError, ConnectionError, RequestException, RuntimeError, IOError) as e:
                    raise ChemblSearchError(chemblid, "similarity", e)

                if len(res) > 0:
                    done_chembl = set([id_to_chembl_id(d) for d in self.done])
                    for r in res:
                        if round(float(r['similarity']), 6) > 99.9:
                            continue  # Dont want exact matches

                        if r['molecule_chembl_id'] not in done_chembl:
                            try:
                                hits.append(
                                    ChemblHit(r['molecule_structures']['canonical_smiles'], r['molecule_chembl_id']))

                            except ChemblParseError as e:
                                print(f"Skipping {r['molecule_chembl_id']} - {e}")
                                self.done.add(r['molecule_chembl_id'])

                        if len(hits) == self.args.max_hit_structures:
                            break
            if round(threshold, 6) <= min_threshold:
                self.best_poses.pop()
        return hits

    def _search_csd(self):
        '''
        Similarity search the csd. Will search with sequentially more relaxes settings to fins
        new as-yet-untried ligands up to a similarity threshold as specified on the command line.

        :return: found hits
        '''

        threshold = 1.0
        min_threshold = self.args.minimum_similarity / 100.0
        hits = []

        while not hits and len(self.best_poses) > 0:
            mol = self.best_poses[-1][1].molecule
            mol.standardise_aromatic_bonds()
            while not hits and round(threshold, 6) > min_threshold:
                # initially just search for CSD structures similar to the ligand
                threshold -= 0.1
                simsearcher = SimilaritySearch(mol, threshold)
                simsearcher.settings.only_organic = True
                simsearcher.settings.has_3d_coordinates = True
                hits = [h for h in simsearcher.search(max_hit_structures=self.args.max_hit_structures)
                        if h.identifier not in self.done]
                print('Found %s structures similar to %s with threshold %.2f.' %
                      (len(hits), mol.identifier, threshold))
                if round(threshold, 6) == min_threshold:
                    self.best_poses.pop()

                if not hits:
                    print('No hits with threshold %.2f against %s' % (threshold, mol.identifier))
                    continue

                # then further specify the search by finding structures
                # similar to individual components of hits
                def match_components(h):
                    if len(h.molecule.components) > 1:
                        for m in h.molecule.components:
                            m.identifier = h.identifier
                            hh = simsearcher.search(m)
                            if hh:
                                return hh[0]
                        return None
                    return h

                hits = [x for x in [match_components(h) for h in hits] if x is not None]
                print('%d new hits with threshold %.2f against %s' %
                      (len(hits), threshold, mol.identifier))

                return hits

    def run_substrate(self):

        """Runs the substrate, for the comparison."""

        if os.path.exists(self.args.ligand):
            self.substrate = EntryReader(self.args.ligand)[0]
        elif (self.args.ligand.lower().startswith("chembl")):
            mol = self.get_chembl_molecule(self.args.ligand).molecule
            self.substrate = Entry.from_molecule(mol)
        else:
            self.substrate = Entry.from_molecule(ChemblHit(self.args.ligand, self.args.ligand).molecule)

        self.start_point_writer.write(self.substrate)

        self.substrate.identifier = self.substrate.identifier.replace(':', '_')
        poses = self.session.dock(self.substrate)

        if not poses:
            print('FAILED TO DOCK SUBSTRATE')
            sys.exit()
        self.best_score, best_pose = self.best(poses)
        self.best_poses.append((self.best_score, best_pose))
        self.substrate_score = self.best_score
        for inx, p in enumerate(poses):
            outfile_path = os.path.join(self.settings.output_directory,
                                        'gold_soln_CSD_API_m1_%d.mol2' % (inx + 1))
            with EntryWriter(outfile_path) as output_file:
                output_file.write(p)
        print('%-9s: %.3f' % ('Substrate', self.best_score))

    def best(self, poses):
        """Extracts best score and best structure from the docked poses."""
        fitness_tag = dict(
            plp='Gold.PLP.Fitness',
            chemscore='Gold.Chemscore.Fitness',
            goldscore='Gold.Goldscore.Fitness',
            asp='Gold.ASP.Fitness'
        )[self.settings.fitness_function]
        scores = [
            float(e.attributes[fitness_tag]) for e in poses
        ]
        best_score = max(scores)
        best_pose = poses[scores.index(best_score)]

        # Post normalise the score by some method. Here we just divide by N-heavy-atoms - a sort of crude
        # 'ligand efficiency' measure. > to avoid counting both dummy atoms and hydrogens

        return best_score / math.sqrt(len([a for a in best_pose.molecule.atoms if a.atomic_number > 1])), best_pose

    def write_poses(self, poses, dock_index):
        for inx, p in enumerate(poses):
            with EntryWriter(
                    os.path.join(
                        self.settings.output_directory,
                        'gold_soln_CSD_API_m%d_%d.mol2' % (dock_index, inx + 1)
                    )
            ) as w:
                w.write(p)

    def write_scores(self, score, hit, pose):
        if score > self.best_score:
            print('%-8s: %.3f Winner!' % (hit.identifier, score))
            self.best_score, best_pose = score, pose
            self.best_poses.append((score, best_pose))
            self.winner_writer.write(best_pose)
        else:
            print('%-8s: %.3f' % (hit.identifier, score))
            if score > self.substrate_score:
                self.better_writer.write(pose)

            if score > self.args.start_point_threshold * self.best_score:
                # Insert it before the best one
                print(
                    f'Adding as potential candidate for next cycle - within {self.args.start_point_threshold * 100.0}% of best scoring pose.')
                self.best_poses.append((score, pose))
                self.start_point_writer.write(pose)


########################################################################################


if __name__ == '__main__':
    similarity_docking = SimilarityDocking()
    similarity_docking.run()
