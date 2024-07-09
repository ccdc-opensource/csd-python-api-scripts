#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2016-12-15: created by S. B. Wiggin, the Cambridge Crystallographic Data Centre
# 2024-07-02: minor update to include using ccdc utilities to find the solvent file

"""
Script to identify and remove bound solvent molecules from a MOF structure.

Solvents are identified using a defined list.
Output in CIF format includes only framework component with all monodentate solvent removed.
"""
#######################################################################

import os
import glob
import argparse

from ccdc import io
from ccdc import utilities

#######################################################################

arg_handler = argparse.ArgumentParser(description=__doc__)
arg_handler.add_argument(
    'input_file',
    help='CSD .gcd file from which to read MOF structures'
)
arg_handler.add_argument(
    '-o', '--output-directory',
    help='Directory into which to write stripped structures'
)
arg_handler.add_argument(
    '-m', '--monodentate', default=False, action='store_true',
    help='Whether or not to strip all unidenate (or monodentate) ligands from the structure'
)
arg_handler.add_argument(
    '-s', '--solvent-file',
    help='Location of solvent file'
)

args = arg_handler.parse_args()
if not args.output_directory:
    args.output_directory = os.path.abspath(os.path.dirname(args.input_file))

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

# Define the solvent smiles patterns
if not args.solvent_file:
    args.solvent_file = utilities.Resources().get_ccdc_solvents_dir()

if os.path.isdir(args.solvent_file):
    solvent_smiles = [
        io.MoleculeReader(f)[0].smiles
        for f in glob.glob(os.path.join(args.solvent_file, '*.mol2'))
    ]
else:
    solvent_smiles = [m.smiles for m in io.MoleculeReader(args.solvent_file)]


#######################################################################


def is_multidentate(c, mol):
    """
    Check for components bonded to metals more than once.
    If monodentate is not specified in the arguments, skip this test.
    """
    if not args.monodentate:
        return True
    got_one = False
    for a in c.atoms:
        orig_a = mol.atom(a.label)
        if any(x.is_metal for b in orig_a.bonds for x in b.atoms):
            if got_one:
                return True
            got_one = True
    return False


def is_solvent(c):
    """Check if this component is a solvent."""
    return c.smiles == 'O' or c.smiles in solvent_smiles


def has_metal(c):
    """Check if this component has any metals."""
    return any(a.is_metal for a in c.atoms)


# Iterate over entries
try:
    for entry in io.EntryReader(args.input_file):
        if entry.has_3d_structure:
            # Ensure labels are unique
            mol = entry.molecule
            mol.normalise_labels()
            # Use a copy
            clone = mol.copy()
            # Remove all bonds containing a metal atom
            clone.remove_bonds(b for b in clone.bonds if any(a.is_metal for a in b.atoms))
            # Work out which components to remove
            to_remove = [
                c
                for c in clone.components
                if not has_metal(c) and (not is_multidentate(c, mol) or is_solvent(c))
            ]
            # Remove the atoms of selected components
            mol.remove_atoms(
                mol.atom(a.label) for c in to_remove for a in c.atoms
            )
            # Write the CIF
            entry.crystal.molecule = mol
            with io.CrystalWriter('%s/%s_stripped.cif' % (args.output_directory, entry.identifier)) as writer:
                writer.write(entry.crystal)
except RuntimeError:
    print('File format not recognised')
