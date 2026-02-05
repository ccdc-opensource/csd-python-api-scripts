# Conformer Filter Density

Filter conformers using a variety of metrics described below.

A csv file will be produced with various points of analysis for unusual torsions. Optionally, files containing the
 conformers which pass and fail within the limits can be written to files.

CCDC Python API Licence required, minimum version: 3.0.15

Script can be run with any multimolecule file e.g. sdf. 

## Instructions on Running

positional arguments: input molecule file

Input file (single- or multi-molecule file)

options:
* -h, --help; Show this help message and exit
* -m {absolute,relative}, --mode; Limit mode: absolute (fixed threshold) or relative (threshold based on 
molecule with fewest unusual torsions). WARNING: Relative mode may behave unexpectedly with conformers from 
multiple input molecules (default: absolute)
* -l, --limit; Maximum number of unusual torsions for a passing molecule (default: 0)
* -d, --local-density; Local density threshold for classifying a torsion as unusual (default: 10.0)
* --incl-organometallics; Include organometallic compounds in the search (default: organic compounds only)
* --generalisation; Turn on generalisation for searches
* --successfn; Output file for molecules that pass the filter (default: successes.mol)
* --failurefn; Output file for molecules that fail the filter (default: failures.mol)
* -u, --unusual-torsions; Output CSV file for unusual torsion details (default: unusual_torsions.csv)


Originally created by Paul Sanschagrin
Updated by Chris Ringrose

    For feedback or to report any issues please contact support@ccdc.cam.ac.uk

