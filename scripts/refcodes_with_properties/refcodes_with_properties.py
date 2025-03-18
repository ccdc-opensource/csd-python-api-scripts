#!/usr/bin/env python
#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2025-03-14: created by Jason C. Cole, The Cambridge Crystallographic Data Centre

'''
Filter a refcode list to the subset that have the desired properties
'''

#########################################################################

import argparse
import csv
import sys

from ccdc import io

import entry_property_calculator

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-r', '--refcode_file', help='input file containing the list of refcodes', default=None)
    parser.add_argument('-d', '--database_file', help='input file containing the list of refcodes', default=None)
    parser.add_argument('-c', '--control_file', help='configuration file containing the desired properties\n\n %s' % (
        entry_property_calculator.helptext()))
    parser.add_argument('-v', '--get_values', action="store_true",
                        help='calculate and print descriptor values where possible rather than filter\n\n %s' % (
                            entry_property_calculator.helptext()))
    parser.add_argument('-o', '--output_file', default=None,
                        help='output CSV file for results\n\n %s' % (entry_property_calculator.helptext()))

    args = parser.parse_args()

    refcode_file = args.refcode_file
    database_file = args.database_file
    control_file = args.control_file
    print_values = args.get_values

    outfile = sys.stdout
    if args.output_file != None:
        outfile = open(args.output_file, 'wb')

    filterer = entry_property_calculator.parse_control_file(open(control_file, "r").readlines())

    reader = None
    if refcode_file:
        reader = io.EntryReader(refcode_file, format='identifiers')
    elif database_file:
        reader = io.EntryReader(database_file)
    else:
        reader = io.EntryReader('CSD')

    if args.get_values:

        csvwriter = None
        for entry in reader:
            values = filterer.values(entry)
            if csvwriter == None:
                fieldnames = ["identifier"] + values.keys()
                csvwriter = csv.DictWriter(outfile, fieldnames=fieldnames)
                csvwriter.writeheader()
            values["identifier"] = entry.identifier
            csvwriter.writerow(values)

    else:
        for entry in reader:
            if filterer.evaluate(entry):
                outfile.write(entry.identifier + "\n")
