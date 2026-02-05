#!/usr/bin/env python3
#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2025-02-03: created by the Cambridge Crystallographic Data Centre

import argparse
import copy
import csv
import math

from ccdc import conformer, io


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('conformer_file',
                        metavar='<input file>',
                        help='Input file (multi-molecule file)')

    parser.add_argument('-csv', '--write-csv',
                        dest='write_csv',
                        action='store_true',
                        help='Write a csv file for all the analysed conformers.')

    return parser.parse_args()


class ProbabilityScorer:
    """
    Use the ConformerGenerator and GeometryAnalyser to score conformers based on their conformational
    probabilities and unusual torsions.
    """
    def __init__(self, user_conformer_generator_settings=None, skip_minimisation=True,
                 user_mogul_analysis_settings=None):

        self._generator = self._create_conformer_generator(user_conformer_generator_settings, skip_minimisation)
        self._mogul_analysis_engine = self._create_analyser(user_mogul_analysis_settings)

    def _create_analyser(self, user_mogul_analyser_settings):
        """Create a GeometryAnalyser engine to analyse the conformers."""
        engine = conformer.GeometryAnalyser()

        # Tweak the 'default defaults'
        # By default, in this use case, we do not want to use generalisation.
        engine.settings.generalisation = False
        # By default, only the organic subset, i.e. exclude organometallic
        engine.settings.organometallic_filter = 'Organic'

        # Ensure user settings are kept:
        if user_mogul_analyser_settings is not None:
            engine.settings = copy.copy(user_mogul_analyser_settings)

        # Only analyse torsions:
        engine.settings.bond.analyse = False
        engine.settings.angle.analyse = False
        engine.settings.ring.analyse = False

        return engine

    def _create_conformer_generator(self, user_generator_settings, skip_minimisation=True):
        """Create a ConformerGenerator engine to generate conformers for the molecules."""
        settings = conformer.ConformerSettings()
        if user_generator_settings is not None:
            settings = copy.copy(user_generator_settings)

        # Mandatory setting here: this must work in use_input_torsion_distributions mode.
        settings.use_input_torsion_distributions = True

        engine = conformer.ConformerGenerator(settings, skip_minimisation)
        return engine

    @staticmethod
    def _combined_local_density(all_local_densities):
        if all_local_densities:
            return sum(all_local_densities) / len(all_local_densities)
        return 100.0

    def probability_analysis(self, molecule, bad_normalised_score=1.0, bad_probability_score=0.0,
                             bad_rmsd_score=9999.9):

        # Approximation excluding rings:
        is_rigid = sum(bond.is_rotatable for bond in molecule.bonds) == 0

        conformers = self._generator.generate(molecule)

        normalised_score = None
        rmsd = None
        ln_probability = None

        if conformers:
            # First conformer is the most likely, so return its score.
            normalised_score = conformers[0].normalised_score
            rmsd = conformers[0].rmsd()
            ln_probability = conformers[0].probability

            if is_rigid:
                if ln_probability is None:
                    ln_probability = 0.0
                if normalised_score is None:
                    normalised_score = 0.0

        if normalised_score is None:
            normalised_score = bad_normalised_score

        if rmsd is None:
            rmsd = bad_rmsd_score

        probability = bad_probability_score
        if ln_probability is not None:
            probability = math.exp(ln_probability)

        return normalised_score, probability, rmsd

    def unusual_torsions_analysis(self, molecule, max_hist_size=15):
        checked_mol = self._mogul_analysis_engine.analyse_molecule(molecule)
        unusual_count = 0
        local_densities = []
        unusual_local_densities = []
        no_data_torsions = 0
        for tor in checked_mol.analysed_torsions:
            hist_size = sum(tor.histogram())
            if hist_size > max_hist_size:
                local_densities.append(tor.local_density)

                if tor.unusual:
                    unusual_local_densities.append(tor.local_density)
                    unusual_count += 1
            else:
                no_data_torsions += 1

                # Expected number of torsions within 10 degrees assuming even distribution of observations.
                # Number of Mogul bins would be strictly better than assuming 18.
                uniform_torsion_dist = 100 / 18
                local_densities.append(uniform_torsion_dist)
        return unusual_count, self._combined_local_density(local_densities), self._combined_local_density(
            unusual_local_densities), no_data_torsions

    def process_molecule(self, molecule):
        molecule.remove_unknown_atoms()
        molecule.assign_bond_types()

        unusual_count, avg_local_density, avg_unusual_local_density, no_data_torsions = self.unusual_torsions_analysis(
            molecule)
        normalised_score, probability, rmsd = self.probability_analysis(molecule)

        return {
            "identifier": molecule.identifier,
            "number of unusual torsions": unusual_count,
            "average local density": avg_local_density,
            "average unusual local density": avg_unusual_local_density,
            "number of torsions with no data": no_data_torsions,
            "normalised probability score": normalised_score,
            "probability score": probability,
            "RMSD to input conformation": rmsd,
        }


def run():
    args = parse_args()

    p = ProbabilityScorer()

    all_data = []

    with io.MoleculeReader(args.conformer_file) as reader:
        for i, molecule in enumerate(reader):
            data = p.process_molecule(molecule)
            if i == 0:
                print(", ".join(data.keys()))
            print(", ".join(str(value) for value in data.values()))
            all_data.append(data)

    if args.write_csv:
        keys = all_data[0].keys()
        with open('filtered_poses_analysis.csv', 'w', newline='') as output_file:
            dict_writer = csv.DictWriter(output_file, keys)
            dict_writer.writeheader()
            dict_writer.writerows(all_data)


if __name__ == "__main__":
    run()
