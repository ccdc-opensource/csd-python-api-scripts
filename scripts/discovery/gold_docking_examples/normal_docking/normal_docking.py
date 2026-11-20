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

PROTEIN_FILE = os.path.join(HERE, '1Y6B_protein.mol2')
REFERENCE_FILE = os.path.join(HERE, '1Y6B_reference.mol2')
OUTPUT_DIRECTORY = os.path.join(HERE, 'output')

settings = docking.Docker.Settings()
settings.add_protein_file(PROTEIN_FILE)

# The binding site is defined from the crystallographic ligand: every protein
# atom within 6 A of it becomes part of the site GOLD searches.
protein = settings.proteins[0]
reference = io.MoleculeReader(REFERENCE_FILE)[0]
settings.binding_site = settings.BindingSiteFromLigand(protein, reference, 6.0)

# Here the ligand being docked is the crystal ligand itself; GOLD discards the
# input conformation and rebuilds it.
settings.add_ligand_file(REFERENCE_FILE, ndocks=10)

# Setting a reference ligand makes GOLD record the RMSD of every pose to it,
# as the Gold.Reference.RMSD attribute read back below.
settings.reference_ligand_file = REFERENCE_FILE

settings.fitness_function = 'plp'
settings.autoscale = 100.  # search effort, as a percentage of the default
settings.output_directory = OUTPUT_DIRECTORY
settings.output_format = 'mol2'
settings.output_file = 'docked_poses.mol2'

docker = docking.Docker(settings=settings)
docker.dock(file_name=os.path.join(OUTPUT_DIRECTORY, 'normal_docking.conf'))

print(f'{"Pose":<28} {"Fitness":>10} {"RMSD":>10}')
for pose in docker.results.ligands:
    rmsd = float(pose.attributes['Gold.Reference.RMSD'])
    print(f'{pose.identifier:<28} {pose.fitness():>10.2f} {rmsd:>10.2f}')
