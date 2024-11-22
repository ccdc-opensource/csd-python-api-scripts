#! /usr/bin/env python
########################################################################################################################
#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2024-11-22: created by the Cambridge Crystallographic Data Centre
#
########################################################################################################################

from ccdc import conformer, descriptors, io, molecule
from ccdc.search import SubstructureSearch, SMARTSSubstructure


def read(molecule_file: str) -> molecule:
    print(f'Reading file: {molecule_file} ... ', end='')
    mol_reader = io.MoleculeReader(molecule_file)
    mol = mol_reader[0]
    print('done.')

    return mol


def generate_conformers(molecule: molecule, max_conformers: int = 50) -> conformer.ConformerHitList:
    """
    Generate conformers for a molecule.

    :param molecule: The Molecule (ccdc Molecule object) to generate conformers for.
    :param max_conformers: The maximum number of conformers to generate.

    :returns: ccdc.conformer.ConformerHitList
    """

    # Set up the ConformerGenerator
    confgen = conformer.ConformerGenerator()
    confgen.settings.max_conformers = max_conformers
    # confgen.settings.superimpose_conformers_onto_reference = True

    # Generate conformers and assign identifiers to them before returning
    conformers = confgen.generate(molecule)

    print(f'Generating conformers, maximum of {max_conformers} ... ', end='')
    for i, conf in enumerate(conformers):
        conf.molecule.identifier = '{}_{:04}'.format(conf.molecule.identifier, i + 1)
    print(f'done, generated {len(conformers)} conformers.')

    return conformers


def analyse(conformers: conformer.ConformerHitList) -> molecule:
    """
    Perform some basic analysis of the conformers generated.
    :param conformers: Conformers generated from ConfGen
    :return: The best molecule of all the conformers generated.
    """
    print(f'Sampling limit reached? {"Yes." if conformers.sampling_limit_reached else "No."}')

    print(f'How many rotamers had no observations? {conformers.n_rotamers_with_no_observations}.')

    most_probable_conformer = conformers[0]

    print(f'Normalised score of most probable conformer: {round(most_probable_conformer.normalised_score, 5)}.')
    print(f'Most probable conformer RMSD wrt input: {round(most_probable_conformer.rmsd(), 3)}; '
          f'wrt minimised: {round(most_probable_conformer.rmsd(wrt="minimised"), 3)}.')

    print(f'Scores of top 10 conformers: ', end='')

    top_ten = conformers[:10]
    for i in range(len(top_ten)):
        if i < len(top_ten) - 1:
            print(f'{round(top_ten[i].normalised_score, 3):.3f}, ', end='')
        else:
            print(f'{round(top_ten[i].normalised_score, 3):.3f}.')

    return most_probable_conformer.molecule


def overlay(conformers, query: str, output_filename: str) -> None:
    """
    Overlay conformers based on a SMARTS substructure pattern
    :param conformers: Conformers generated from ConfGen
    :param query: SMARTS pattern which the conformers will overlay on top of.
        Should be consistent across all conformers, e.g. benzene ring.
    """
    print(f'Overlaying conformers ... ', end='')
    conformers_mols = [c.molecule for c in conformers]
    ss_search = SubstructureSearch()
    substructure = SMARTSSubstructure(query)
    ss_search.add_substructure(substructure)
    hits = ss_search.search(conformers_mols, max_hits_per_structure=1)
    ref_ats = hits[0].match_atoms()
    print('done.')

    print('Writing file superimposed ... ', end='')
    with io.MoleculeWriter(output_filename) as writer:
        for hit in hits:
            hit_ats = hit.match_atoms()
            atoms = zip(ref_ats, hit_ats)
            ov = descriptors.MolecularDescriptors.Overlay(hits[0].molecule, hit.molecule, atoms)
            superimposed_hit = ov.molecule
            writer.write(superimposed_hit)
    print('done.')


def write_conformers_to_file(conformers: conformer.ConformerHitList, filename: str) -> None:
    """
    Write conformers to a file without any addition overlaying.
    :param conformers: Conformer generated from ConfGen.
    :param filename: The name of the output file.
    """

    with io.MoleculeWriter(filename) as writer:
        for conformer in conformers:
            writer.write(conformer.molecule)


if __name__ == '__main__':

    input_filename = 'AZD9291.mol2'
    # Read example molecule
    mol = read(input_filename)

    # Generate conformers
    conformers = generate_conformers(mol, 20)

    # Provide summary of analysis
    analyse(conformers)

    # Overlay structures based on common substructure
    query = 'c1cncnc1'
    output_filename = f'superimposed_{input_filename}'
    overlay(conformers, query, output_filename)
