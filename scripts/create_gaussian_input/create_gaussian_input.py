#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2015-06-17: created by the Cambridge Crystallographic Data Centre
#

from __future__ import division, absolute_import, print_function

"""
This script will generate a generic Gaussian input file
Input: CSD Identifier as a string or .mol2 
Output: GJF input file
Author: Andy Maloney 2015
"""

import sys
import os
import ccdc.io
import glob


def fatal(*args):
    """Generates an error message if necessary to smoothly exit the program."""
    print('ERROR:', ' '.join(map(str, args)))
    sys.exit(1)


def file_writer(molecule, name):
    """Writes a standard Gaussian input file for all molecules contained in the structure files."""
    if not mol.all_atoms_have_sites:
        fatal(entry_id, 'has some atoms without coordinates')
    mol.normalise_hydrogens()

    for i, component in enumerate(molecule.components):

        file_name = '%s_molecule%d.gjf' % (name, i)
        f = open(file_name, 'w')

        f.write('#B3LYP/6-31G**  opt\n')
        f.write('\n')
        f.write('Standard Gaussian Input File for %s, molecule %d\n' % (name, i))
        f.write('\n')
        f.write('0 1\n')

        for atom in component.atoms:
            f.write('%2s  %9.6f  %9.6f  %9.6f\n' % (atom.atomic_symbol,
                                                    atom.coordinates.x,
                                                    atom.coordinates.y,
                                                    atom.coordinates.z))

        f.write('\n')
        f.write('--Link1--')
        f.write('\n')
        f.write('\n')
        f.write('\n')

        f.close()


if __name__ == '__main__':
    # Get the relevant structure typed by user
    if len(sys.argv) != 2:
        fatal('you must supply a structure identifier.')
    entry_id = sys.argv[1]

    # Checking the current directory for user cif file
    filepath = '%s' % entry_id
    if os.path.isfile(filepath):
        reader = ccdc.io.MoleculeReader(filepath)
        for mol in reader:
            identifier = mol.identifier
            file_writer(mol, identifier)

    else:
        # Read molecule from database
        reader = ccdc.io.MoleculeReader('CSD')
        identifier = entry_id
        try:
            mol = reader.molecule(entry_id)
            file_writer(mol, identifier)
        except RuntimeError:
            fatal(entry_id, '- structure not found')