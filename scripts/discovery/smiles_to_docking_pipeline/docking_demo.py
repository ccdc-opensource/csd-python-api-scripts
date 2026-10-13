#! /usr/bin/env python
########################################################################################################################
#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2026-10-08: created by the Cambridge Crystallographic Data Centre
#
########################################################################################################################

from ccdc import conformer, docking, io
from ccdc.descriptors import MolecularDescriptors


example_ligands = {
    "Gefitinib": r"C1COCCN1CCCOc2c(OC)cc3ncnc(c3c2)Nc4cc(Cl)c(F)cc4",
    "Erlotinib": r"COCCOc1cc2c(cc1OCCOC)ncnc2Nc3cccc(c3)C#C",
    "Afatinib": r"CN(C)C\C=C\C(=O)Nc3cc1c(Nc(cc2Cl)ccc2F)ncnc1cc3OC4COCC4",
}


def dock_single_ligand(protein_file, ligand_file, results_dir='results'):
    """Docks a single ligand from a file to a protein from a file using GOLD."""
    docker = docking.Docker()

    settings = docker.settings

    settings.add_protein_file(protein_file)
    settings.add_ligand_file(ligand_file)

    ALA149 = [
        res for res in settings.proteins[0].residues
        if res.identifier == 'A:ALA149'
    ][0]

    settings.binding_site = settings.BindingSiteFromResidue(
        settings.proteins[0], ALA149, 8.0
    )

    settings.output_directory = results_dir

    settings.fitness_function = 'ChemScore'

    results = docker.dock()
    return results


def smiles_to_3d(smiles, name, n_conformers=1):
    """Writes a 3D molecule file from a SMILES string."""
    mol = io.Molecule.from_string(smiles)
    mol.identifier = name

    conf_gen = conformer.ConformerGenerator()
    conf_gen.settings.max_conformers = n_conformers

    conformers = conf_gen.generate(mol)

    with io.MoleculeWriter(f'{name}.mol') as writer:
        for conf in conformers:
            writer.write(conf.molecule)


def overlay_separate_files(*mol_files):
    """Overlays similar ligands from separate molecule files."""
    overlay = MolecularDescriptors.LigandOverlay('.')

    overlay.add_ligands(
        [io.MoleculeReader(mol_file)[0] for mol_file in mol_files]
    )

    overlay.overlay_ligands()


if __name__ == '__main__':
    for name, smiles in example_ligands.items():
        smiles_to_3d(smiles, name)
        dock_single_ligand('data/4hzr_prepared.cif', f'{name}.mol', f'{name}_results')
    overlay_separate_files('Gefitinib.mol', 'Erlotinib.mol', 'Afatinib.mol')
