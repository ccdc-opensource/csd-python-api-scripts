"""
This script can be used for any purpose without limitation subject to the
conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
This permission notice and the following statement of attribution must be
included in all copies or substantial portions of this script.

"17/11/2026": created by Chris Ringrose for the Cambridge Crystallographic Data Centre
"""

import os

from ccdc import docking, io

HERE = os.path.dirname(os.path.abspath(__file__))

PROTEIN_FILE = os.path.join(HERE, '1KE7_protein.mol2')
REFERENCE_FILE = os.path.join(HERE, '1KE7_reference.mol2')
OUTPUT_DIRECTORY = os.path.join(HERE, 'output')

settings = docking.Docker.Settings()
settings.add_protein_file(PROTEIN_FILE)

protein = settings.proteins[0]
reference = io.MoleculeReader(REFERENCE_FILE)[0]
settings.binding_site = settings.BindingSiteFromLigand(protein, reference, 6.0)

settings.add_ligand_file(REFERENCE_FILE, ndocks=10)
settings.reference_ligand_file = REFERENCE_FILE

# The point is placed on the hinge hydrogen bond: in 1KE7 the Leu83 backbone
# NH donates to the LS3 oxygen O5 at the coordinates below (N...O 2.68 A).
#
# Place the point where you want the ligand ATOM to sit, not on the protein
# partner, and make the radius large enough that a good pose can reach it. A
# point set a few angstroms off a real acceptor position with a tight radius
# will quietly push docking away from the correct answer rather than towards it.
#
# Valid types: 'acceptor', 'donor_acceptor', 'donor', 'aromatic_center'
#              (or 'aromatic_centre'), 'ring_center' (or 'ring_centre').
hinge_acceptor = settings.PharmacophoreConstraint(
    origin=(-10.910, 52.354, 11.763),
    type='acceptor',
    radius=1.5,
    function_shape='gaussian',   # reward decays with distance, rather than cutting off
    use_as_fitting_point=True,   # also bias where GOLD tries to place the ligand
    weight=10.0,
)
settings.add_constraint(hinge_acceptor)

settings.fitness_function = 'plp'
settings.autoscale = 100.
settings.output_directory = OUTPUT_DIRECTORY
settings.output_format = 'mol2'
settings.output_file = 'docked_poses.mol2'

docker = docking.Docker(settings=settings)
docker.dock(file_name=os.path.join(OUTPUT_DIRECTORY, 'pharmacophore_constraint.conf'))

print(f'{"Pose":<28} {"Fitness":>10} {"RMSD":>10}')
for pose in docker.results.ligands:
    rmsd = float(pose.attributes['Gold.Reference.RMSD'])
    print(f'{pose.identifier:<28} {pose.fitness():>10.2f} {rmsd:>10.2f}')
