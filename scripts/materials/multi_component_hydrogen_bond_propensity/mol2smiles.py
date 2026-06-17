#!/usr/bin/env python
#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#

from ccdc.io import MoleculeReader
import os
import csv
import argparse
from utilities import file_list, string_scrubber, read_experimental_csv


def read_mol_file(directory, file):
    '''Returns: identifier, smiles'''
    mol_reader = MoleculeReader(os.path.join(directory, file))
    mol = mol_reader[0]
    
    return mol.identifier, mol.heaviest_component.smiles


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument(
        "--input_dir",
        type=str,
        required=True,
        help="Directory containing API folders."
    )
    parser.add_argument(
        "--output_filename",
        type=str,
        required=True,
        help="Filename of formatted .csv file containing SMILES."
    )
    parser.add_argument(
        "--experimental_csv",
        type=str,
        required=False,
        help="Filename of formatted .csv file containing identifier names and experimental bool."
    )
    parser.add_argument(
        '--clean_id',
        action='store_true',
        help='Removes special characters from ids that may be problematic.'
    )
    
    args = parser.parse_args()

    output_path = os.path.join(args.input_dir, args.output_filename)
    with open(output_path, 'w', newline='', encoding="utf-8") as output_file:
        csvwriter = csv.writer(output_file, delimiter=',', quotechar='|')
        csvwriter.writerow(['identifier', 'n_components',
                           'component_a', 'component_b', 'neutral_a', 'neutral_b'])
    
    if args.experimental_csv:
        experimental_dict = read_experimental_csv(args.experimental_csv)

    # API group directories contain one or more API files and a directory of coformers
    API_groups = [name for name in os.listdir(
        args.input_dir) if os.path.isdir(os.path.join(args.input_dir, name))]
    
    exp_replaced = combo_count = 0
    
    with open(output_path, 'a+', newline='', encoding="utf-8") as output_file:
        
        for API_group in API_groups:
            csvwriter = csv.writer(output_file, delimiter=',', quotechar='|')
            API_group_path = os.path.join(args.input_dir, API_group)
            
            for API_file in file_list(API_group_path):
                api_id, api_smiles = read_mol_file(
                    API_group_path, API_file)
                print(api_id)
                coformer_dir_path = os.path.join(API_group_path, 'coformers')
                
                for coformer_file in file_list(coformer_dir_path):
                    coformer_id, coformer_smiles = read_mol_file(
                        coformer_dir_path, coformer_file)
                    
                    combo_count += 1
                    exp_bool = "?"
                    # Try to look up the experimental boolean in dictionary, if provided
                    if args.experimental_csv:
                        for (x, y) in [(api_id, coformer_id), (coformer_id, api_id)]:
                            if (x, y) in list(experimental_dict.keys()):
                                exp_bool = experimental_dict[(x, y)]
                                exp_replaced += 1
                    
                    # Clean the ids if the option is turned on
                    if args.clean_id:
                        api_id = string_scrubber(api_id)
                        coformer_id = string_scrubber(coformer_id)
                    
                    n_components = 2
                    if api_smiles == coformer_smiles:
                        n_components = 1
                    
                    combo_id = ".".join([api_id, coformer_id, str(exp_bool)])
                    csvwriter.writerow([f'"{combo_id}"', n_components, api_smiles, coformer_smiles, "", ""])
    
    if args.experimental_csv:
        print(f"Found experimental labels for {exp_replaced} out of {combo_count} combinations")


if __name__ == '__main__':
    main()
