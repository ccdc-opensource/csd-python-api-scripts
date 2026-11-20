"""
This script can be used for any purpose without limitation subject to the
conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
This permission notice and the following statement of attribution must be
included in all copies or substantial portions of this script.

"17/11/2026": created by Chris Ringrose for the Cambridge Crystallographic Data Centre
"""

import os

from ccdc import docking, io
from ccdc.descriptors import MolecularDescriptors

HERE = os.path.dirname(os.path.abspath(__file__))

PROTEIN_FILE = os.path.join(HERE, '2VU3_protein.mol2')
TEMPLATE_FILE = os.path.join(HERE, '2VTQ_reference.mol2')   # analogue, crystal pose
LIGAND_FILE = os.path.join(HERE, '2VU3_reference.mol2')     # AT7519
OUTPUT_DIRECTORY = os.path.join(HERE, 'output')


def shared_core(template, ligand):
    """The maximum common substructure, in the template's crystal coordinates."""
    mcs = MolecularDescriptors.MaximumCommonSubstructure()
    mcs.settings.check_bond_type = False   # the two pyrazoles are perceived differently
    mcs.settings.ignore_hydrogens = True
    matched_pairs, _ = mcs.search(template, ligand)

    core = template.copy()
    keep = {template_atom.index for template_atom, _ in matched_pairs}
    core.remove_atoms([a for a in core.atoms if a.index not in keep])
    core.identifier = 'shared_core'
    return core


template = io.MoleculeReader(TEMPLATE_FILE)[0]
ligand = io.MoleculeReader(LIGAND_FILE)[0]

core = shared_core(template, ligand)
heavy_atoms = sum(1 for a in core.atoms if a.atomic_symbol != 'H')
print(f'Scaffold: {heavy_atoms} heavy atoms')

settings = docking.Docker.Settings()
settings.add_protein_file(PROTEIN_FILE)

protein = settings.proteins[0]
settings.binding_site = settings.BindingSiteFromLigand(protein, ligand, 6.0)

settings.add_ligand_file(LIGAND_FILE, ndocks=10)
settings.reference_ligand_file = LIGAND_FILE

settings.add_constraint(settings.ScaffoldMatchConstraint(core, weight=5.0))

settings.fitness_function = 'plp'
settings.autoscale = 100.
settings.output_directory = OUTPUT_DIRECTORY
settings.output_format = 'mol2'
settings.output_file = 'docked_poses.mol2'

docker = docking.Docker(settings=settings)
docker.dock(file_name=os.path.join(OUTPUT_DIRECTORY, 'scaffold_constraint.conf'))

print(f'{"Pose":<28} {"Fitness":>10} {"RMSD":>10}')
for pose in docker.results.ligands:
    rmsd = float(pose.attributes['Gold.Reference.RMSD'])
    print(f'{pose.identifier:<28} {pose.fitness():>10.2f} {rmsd:>10.2f}')
