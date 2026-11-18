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

PROTEIN_FILE = os.path.join(HERE, '3HS4_protein.mol2')
REFERENCE_FILE = os.path.join(HERE, '3HS4_reference.mol2')
LIGAND_FILE = os.path.join(HERE, 'acetazolamide_anion.mol2')
OUTPUT_DIRECTORY = os.path.join(HERE, 'output')


def deprotonate_sulfonamide(molecule):
    """Turn a neutral sulfonamide NH2 into the anion that coordinates zinc."""
    nitrogen = next(
        atom for atom in molecule.atoms
        if atom.atomic_symbol == 'N'
        and any(n.atomic_symbol == 'S' for n in atom.neighbours)
        and sum(1 for n in atom.neighbours if n.atomic_symbol == 'H') == 2
    )
    hydrogen = next(n for n in nitrogen.neighbours if n.atomic_symbol == 'H')
    molecule.remove_atom(hydrogen)
    nitrogen.formal_charge = -1
    return molecule


def report_metal_geometry(log_file):
    """
    Print the metal coordination geometry reported by GOLD in the protein log.
    This confirms the metal ion was used in docking.
    """
    wanted = ('Coordinating residues', 'Matched coordination geometry',
              'Metal fitting points added')
    with open(log_file) as handle:
        lines = handle.readlines()
    for index, line in enumerate(lines):
        if any(key in line for key in wanted):
            print(line.rstrip())
            # These entries put their value on the following line.
            if line.rstrip().endswith(':') and index + 1 < len(lines):
                print(f'    {lines[index + 1].strip()}')


ligand = deprotonate_sulfonamide(io.MoleculeReader(REFERENCE_FILE)[0])
ligand.identifier = 'acetazolamide_anion'
with io.MoleculeWriter(LIGAND_FILE) as writer:
    writer.write(ligand)

settings = docking.Docker.Settings()
settings.add_protein_file(PROTEIN_FILE)

protein = settings.proteins[0]
print(f'Metals kept: {[a.atomic_symbol for a in protein.metals]}')

settings.binding_site = settings.BindingSiteFromLigand(protein, ligand, 6.0)
settings.add_ligand_file(LIGAND_FILE, ndocks=10)

# RMSD is measured against the crystallographic pose, which still has its
# original protonation. GOLD compares heavy atoms, so the missing hydrogen
# does not affect the comparison.
settings.reference_ligand_file = REFERENCE_FILE

settings.fitness_function = 'plp'
settings.autoscale = 100.
settings.output_directory = OUTPUT_DIRECTORY
settings.output_format = 'mol2'
settings.output_file = 'docked_poses.mol2'

docker = docking.Docker(settings=settings)
docker.dock(file_name=os.path.join(OUTPUT_DIRECTORY, 'metalloprotein_docking.conf'))

report_metal_geometry(os.path.join(OUTPUT_DIRECTORY, 'gold_protein.log'))

print(f'{"Pose":<28} {"Fitness":>10} {"RMSD":>10}')
for pose in docker.results.ligands:
    rmsd = float(pose.attributes['Gold.Reference.RMSD'])
    print(f'{pose.identifier:<28} {pose.fitness():>10.2f} {rmsd:>10.2f}')
