#!/usr/bin/env python3
#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2026-02-03: created by the Cambridge Crystallographic Data Centre

import argparse
import csv

from ccdc import conformer, io


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('inmolfn',
                        metavar='<input molecule file>',
                        help='Input file (single- or multi-molecule file)')

    parser.add_argument('-m', '--mode',
                        choices=['absolute', 'relative'],
                        default='absolute',
                        help='Limit mode: absolute (fixed threshold) or relative '
                             '(threshold based on molecule with fewest unusual torsions). '
                             'WARNING: Relative mode may behave unexpectedly with conformers '
                             'from multiple input molecules (default: %(default)s)')

    parser.add_argument('-l', '--limit',
                        dest='torsion_limit',
                        type=int,
                        default=0,
                        metavar='<limit>',
                        help='Maximum number of unusual torsions for a passing molecule '
                             '(default: %(default)s)')

    parser.add_argument('-d', '--local-density',
                        dest='local_density_threshold',
                        type=float,
                        default=10.0,
                        metavar='<threshold>',
                        help='Local density threshold for classifying a torsion as unusual '
                             '(default: %(default)s)')

    parser.add_argument('--incl-organometallics',
                        dest='incl_organometallics',
                        action='store_true',
                        help='Include organometallic compounds in the search '
                             '(default: organic compounds only)')

    parser.add_argument('--generalisation',
                        action='store_true',
                        help='Turn on generalisation for searches')

    parser.add_argument('--successfn',
                        default='successes.mol',
                        metavar='<file>',
                        help='Output file for molecules that pass the filter '
                             '(default: %(default)s)')

    parser.add_argument('--failurefn',
                        default='failures.mol',
                        metavar='<file>',
                        help='Output file for molecules that fail the filter '
                             '(default: %(default)s)')

    parser.add_argument('-u', '--unusual-torsions',
                        dest='unusualtorsionsfn',
                        default='unusual_torsions.csv',
                        metavar='<file>',
                        help='Output CSV file for unusual torsion details '
                             '(default: %(default)s)')

    return parser.parse_args()


def create_mogul_engine(local_density_threshold, incl_organometallics, generalisation):
    """Create and configure a geometry analyser engine.

    Args:
        local_density_threshold: Threshold for classifying torsions as unusual
        incl_organometallics: Whether to include organometallic compounds
        generalisation: Whether to enable generalisation for searches

    Returns:
        Configured ccdc.conformer.GeometryAnalyser instance
    """
    engine = conformer.GeometryAnalyser()

    engine.settings.bond.analyse = False
    engine.settings.angle.analyse = False
    engine.settings.ring.analyse = False

    engine.settings.torsion.local_density_threshold = local_density_threshold
    engine.settings.generalisation = generalisation
    engine.settings.organometallic_filter = 'all' if incl_organometallics else 'organics_only'

    return engine


def analysis(torsion_limit, input_filename, mode, engine, success_file, failure_file, unusual_torsion_file):
    """Analyze molecules for unusual torsions and filter based on criteria.

    Args:
        torsion_limit: Maximum number of unusual torsions allowed
        input_filename: Path to input molecule file
        mode: 'absolute' or 'relative' filtering mode
        engine: Configured GeometryAnalyser instance
        success_file: Path to output file for passing molecules
        failure_file: Path to output file for failing molecules
        unusual_torsion_file: Path to CSV file for unusual torsion details
    """
    # Analyze all molecules and collect unusual torsion data
    molecules = []
    min_unusual_torsions = float('inf')

    with io.MoleculeReader(input_filename) as mol_reader:
        for molecule in mol_reader:
            molecule.standardise_aromatic_bonds()
            molecule.standardise_delocalised_bonds()

            geometry_analysed_molecule = engine.analyse_molecule(molecule)

            molecule.unusual_torsions = [
                t for t in geometry_analysed_molecule.analysed_torsions
                if t.unusual and t.enough_hits
            ]
            molecule.num_unusual_torsions = len(molecule.unusual_torsions)
            molecules.append(molecule)

            min_unusual_torsions = min(min_unusual_torsions, molecule.num_unusual_torsions)

    # Write results
    with io.MoleculeWriter(success_file) as passed_writer, \
         io.MoleculeWriter(failure_file) as failed_writer, \
         open(unusual_torsion_file, 'w', newline='') as csv_file:

        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['MoleculeIndex', 'Value', 'Zscore', 'LocalDensity', 'NumHits', 'Atoms'])

        for idx, molecule in enumerate(molecules):
            threshold = torsion_limit if mode == 'absolute' else min_unusual_torsions + torsion_limit
            failed = molecule.num_unusual_torsions > threshold

            if failed:
                failed_writer.write(molecule)
                for torsion in molecule.unusual_torsions:
                    csv_writer.writerow([
                        idx,
                        torsion.value,
                        torsion.z_score,
                        torsion.local_density,
                        torsion.nhits,
                        ' '.join(torsion.atom_labels)
                    ])
            else:
                passed_writer.write(molecule)


def run():
    """Main entry point for the script."""
    args = parse_args()

    if args.torsion_limit < 0:
        raise ValueError('Torsion limit must be >= 0')

    engine = create_mogul_engine(
        args.local_density_threshold,
        args.incl_organometallics,
        args.generalisation
    )

    analysis(
        args.torsion_limit,
        args.inmolfn,
        args.mode,
        engine,
        args.successfn,
        args.failurefn,
        args.unusualtorsionsfn,
    )


if __name__ == '__main__':
    run()
