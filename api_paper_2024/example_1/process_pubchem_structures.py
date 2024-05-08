#!/usr/bin/env python
#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2024-05-02: created by the Cambridge Crystallographic Data Centre
#

import sys

from ccdc.entry import Entry
from ccdc.io import EntryReader, EntryWriter


def process_structures(input_file, output_file):
    with EntryReader(input_file) as er, EntryWriter(output_file) as ew:
        for e in er:
            attribs = e.attributes
            molecule = e.molecule
            molecule.assign_bond_types('all')
            ne = Entry.from_molecule(molecule, attributes=attribs)
            ew.write(ne)


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_structures(input_file, output_file)
