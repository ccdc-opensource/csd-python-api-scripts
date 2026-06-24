#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2015-06-17: created by Andrew Maloney the Cambridge Crystallographic Data Centre
#

"""
This script will generate a generic Gaussian input file
Input: CSD Identifier as a string or .mol2 
Output: GJF input file
"""

import os
import sys

import ccdc.io


def file_writer(molecule, name):
    """Writes a standard Gaussian input file for all molecules contained in the structure files."""
    if not mol.all_atoms_have_sites:
        raise RuntimeError(f'{entry_id} has some atoms without coordinates')
    mol.normalise_hydrogens()

    for i, component in enumerate(molecule.components):

        with open(f'{name}_molecule{i}.gjf', 'w') as f:
            f.write('#B3LYP/6-31G**  opt\n\n')
            f.write(f'Standard Gaussian Input File for {name}, molecule {i}\n\n')
            f.write('0 1\n')

            for atom in component.atoms:
                f.write('%2s  %9.6f  %9.6f  %9.6f\n' % (atom.atomic_symbol,
                                                        atom.coordinates.x,
                                                        atom.coordinates.y,
                                                        atom.coordinates.z))

            f.write('\n--Link1--\n\n\n')


if __name__ == '__main__':
    # Get the relevant structure typed by user
    if len(sys.argv) != 2:
        raise RuntimeError('you must supply a structure identifier.')
    entry_id = sys.argv[1]

    # Checking the current directory for user cif file
    if os.path.isfile(entry_id):
        reader = ccdc.io.MoleculeReader(entry_id)
        for mol in reader:
            identifier = mol.identifier
            file_writer(mol, identifier)

    else:
        # Read molecule from database
        reader = ccdc.io.MoleculeReader('CSD')
        try:
            mol = reader.molecule(entry_id)
            file_writer(mol, entry_id)
        except RuntimeError:
            raise RuntimeError(f'{entry_id} - structure not found')
