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

from ccdc import io
from ccdc import utilities
from mercury_interface import MercuryInterface

#######################################################################

helper = MercuryInterface()
solvent_smiles = []

# Define the solvent smiles patterns
solvent_file = utilities.Resources().get_ccdc_solvents_dir()

if os.path.isdir(solvent_file):
    solvent_smiles = [
        io.MoleculeReader(f)[0].smiles
        for f in glob.glob(os.path.join(solvent_file, '*.mol2'))
    ]

else:
    html_file = helper.output_html_file
    f = open(html_file, "w")
    f.write('<br>')
    f.write('Sorry, unable to locate solvent files in the CCDC directory')
    f.write('<br>')
    f.close()
# a user-defined solvent directory could be added here instead

#######################################################################


def is_solvent(c):
    """Check if this component is a solvent."""
    return c.smiles == 'O' or c.smiles in solvent_smiles


def has_metal(c):
    """Check if this component has any metals."""
    return any(a.is_metal for a in c.atoms)


entry = helper.current_entry
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
        if not has_metal(c) and is_solvent(c)
        ]
    # Remove the atoms of selected components
    mol.remove_atoms(
        mol.atom(a.label) for c in to_remove for a in c.atoms
    )
    # Write the CIF
    entry.crystal.molecule = mol
    with (io.CrystalWriter('%s/%s_stripped.cif' % (helper.options['working_directory_path'], entry.identifier)) as
          writer):
        writer.write(entry.crystal)
        html_file = helper.output_html_file
    f = open(html_file, "w")
    f.write('<br>')
    f.write('Cif file containing MOF framework without monodentate solvent written to your output directory')
    f.write('<br>')
    f.close()
else:
    html_file = helper.output_html_file
    f = open(html_file, "w")
    f.write('<br>')
    f.write('Sorry, this script will only work for CSD entries containing atomic coordinates')
    f.write('<br>')
    f.close()
