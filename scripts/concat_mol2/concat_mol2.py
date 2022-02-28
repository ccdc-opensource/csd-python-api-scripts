#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2014-08-11: created by Peter Galek, The Cambridge Crystallographic Data Centre
# 2022-02-24: updated by Alex Moldovan, The Cambridge Crystallographic Data Centre
#

import glob
import argparse
import os


def main(delete_separate_files):
    mol2_files = glob.glob('*.mol2')
    count = 0
    with open('concat.mol2', 'w') as outfile:
        for f in mol2_files:
            if f == 'concat.mol2':
                continue
            with open(f, 'r') as infile:
                outfile.write(infile.read())
                count += 1
    print(f"{count} files concatenated.")

    if delete_separate_files == True:
        count = 0
        for f in mol2_files:
            if f == 'concat.mol2':
                continue
            os.remove(f)
            count += 1
        print(f"{count} files removed.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-d', '--delete_contributors', action='store_true',
                        help='Remove contributing individual mol2 files after concatenation')

    args = parser.parse_args()
    main(args.delete_contributors)
