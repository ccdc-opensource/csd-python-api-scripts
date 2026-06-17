#! /usr/bin/env python
########################################################################################################################
#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2020-05-07: created by the Cambridge Crystallographic Data Centre
#
########################################################################################################################

import logging
import sys
from argparse import ArgumentParser
from platform import platform
from pathlib import Path
from os import mkdir, chdir
from shutil import copy
from time import time
from dataclasses import dataclass, field
from multiprocessing import Pool

import ccdc
from ccdc.io import EntryReader
from ccdc.docking import Docker

########################################################################################################################
#
# Program parameters...
#

# Default GOLD conf file to use...

CONF_FILE = 'gold.conf'

# Default number of parallel processes to use...

N_PROCESSES = 6

# Note that, although the number of batches is currently the same as the number of processes, it could in theory
# be greater. One  reason to do this might be to ensure that a contiguous block of large, flexible molecules
# in the input file don't make one batch run very much slower than the others. There is a cost to starting up
# new instances of GOLD, however, so care should be taken with this.
#
# The chunking is done such that the sizes are as even as possible, in that the last one will not be much shorter than
# the others. For example, for [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] split into 3 batches, we want
# [[1, 2, 3, 4], [4, 5, 6], [6, 7, 8, 9]] and not [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10]].

########################################################################################################################

# Record type to hold the parameters defining a batch...

@dataclass
class Batch:

    n: int                        # Batch number
    start: int                    # Index of first molecule in batch
    finish: int                   # Index of last molecule in batch
    conf_file: Path               # GOLD configuration file
    output_dir: Path              # Output dir, in which the batch sub-directory will be created
    dir: Path = field(init=False) # Sub-directory for batch, see __post_init__ below

    def __post_init__(self):

        self.dir = self.output_dir / f'chunk_{self.n:02d}'  # Directory will be created in do_chunk function

########################################################################################################################

# A summary of information about the script and where it is running, useful for debugging etc...

SCRIPT_INFO = f"""
Script:          {sys.argv[0]}
Platform:        {platform()}
Python exe:      {sys.executable}
Python version:  {'.'.join(str(x) for x in sys.version_info[:3])}
CSD API version: {ccdc.__version__}
"""

########################################################################################################################

def get_logger(name=__name__):
    logger = logging.getLogger(name)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('[%(asctime)s | %(levelname)s | %(name)s] %(message)s', datefmt='%y-%m-%d %H:%M:%S'))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    return logger


def do_batch(batch):

    """
    Dock a batch of the input file.

    :param batch: a record holding the parameters defining the batch

    As we can't return a GOLD results object from a pool process (it can't be pickled as it wraps C++ objects),
    we simply return a boolean recording whether GOLD exited with a 'success' status code.
    """

    logger = get_logger(f"Batch {batch.n}")

    ######

    # As Settings objects cannot be pickled they cannot be passed to pool processes, so we create a fresh copy...

    settings = Docker.Settings().from_file(str(batch.conf_file))

    # Create and enter the sub-directory for this batch...

    mkdir(batch.dir)

    chdir(batch.dir)

    settings.output_directory = '.'  # Ensure GOLD writes output to the batch sub-directory

    # Specify the batch of molecules to dock...

    ligand_file = settings.ligand_files[0]  # The ligand file info will be overwritten, so store for reference below

    settings.clear_ligand_files()

    settings.add_ligand_file(ligand_file.file_name, ndocks=ligand_file.ndocks, start=batch.start, finish=batch.finish)

    # Run docking...

    logger.info(f"Starting (indices {batch.start} - {batch.finish})...")

    docker = Docker(settings=settings)

    results = docker.dock()

    logger.info(f"...done")

    # As we can't return the results (as they are not picklable) and the poses have already been written to disk, we just return the status code

    return results.return_code

########################################################################################################################

def main():

    """
    Dock the molecules from the supplied input file in parallel.
    """

    logger = get_logger('Main')

    parser = ArgumentParser()

    parser.add_argument('conf_file', nargs='?', default=CONF_FILE, type=str, help=f"GOLD configuration file (default='{CONF_FILE}')")
    parser.add_argument('--n_processes', default=N_PROCESSES, type=int, help=f"No. of processes (default={N_PROCESSES})")

    config = parser.parse_args()

    conf_file = Path(config.conf_file)

    if not conf_file.exists():
        logger.error(f"Error! Configuration file '{conf_file}' not found!", file=sys.stderr)
        sys.exit(1)

    n_processes = config.n_processes

    if not n_processes > 0:
        logger.error(f"Error! Number of processes must be an integer greater than zero.")
        sys.exit(1)

    logger.info(SCRIPT_INFO)

    t0 = time()

    ######

    # Load setting from GOLD conf file...

    settings = Docker.Settings().from_file(str(conf_file))

    # Ensure the output directory exists (a sub-directory for each batch is created within it)...

    output_dir = Path(settings.output_directory)

    if not str(output_dir) == '.':  # Skip directory (re)creation if output dir is current directory

        if output_dir.exists():

            logger.error(f"Error! Output dir '{output_dir}' already exists.")

            sys.exit(1)

        mkdir(output_dir)

    # Count the molecules to dock in the input file...

    input_file = Path(settings.ligand_files[0].file_name)

    with EntryReader(str(input_file)) as reader:

        n_molecules = len(reader)

    logger.info(f"There are {n_molecules} molecules to dock on {n_processes} processes...")

    ######

    # Here we determine the sets of parameters defining the batches; recall that:
    # 1) We choose the number of batches to be the same as the number of processes.
    # 2) GOLD uses 1-based indexing for molecules.

    basic_size, remainder = n_molecules // n_processes, n_molecules % n_processes

    batches, start = [], 1

    for chunk_n in range(1, n_processes+1):

        finish = start + basic_size + (1 if chunk_n <= remainder else 0) - 1

        batches.append(Batch(n=chunk_n, start=start, finish=finish, conf_file=conf_file, output_dir=output_dir))

        start = finish + 1

    ######

    # Dock the batches in parallel...

    with Pool(n_processes) as pool:

        _ = pool.map(do_batch, batches)  # We are not currently checking the return codes

    ######

    # Combine output from batches into output directory and write combined 'bestranking.lst' file...

    bestranking, preamble_and_header = [], None  # Records and preamble plus data header from individual batch 'bestranking.lst' files

    for batch in batches:

        # Copy solution files...

        for soln_file in batch.dir.glob('gold_soln_*'):

            copy(str(soln_file), output_dir)

        # Load records from 'bestranking.lst' file, fixing file paths as we go...

        with (batch.dir / 'bestranking.lst').open('r') as file:

            lines = [x.replace(str(batch.dir), str(output_dir)).replace('\\.\\', '\\') for x in file.read().split('\n')]

        if batch.n == 1:  # Take preamble and data header from first batch only

            preamble_and_header = lines[:7]

        bestranking.extend(lines[7:-1])  # Records

    # Write 'bestranking.lst' file...

    with (output_dir / 'bestranking.lst').open('w') as file:

        file.write('\n'.join(preamble_and_header) + '\n')  # Preamble and data header

        file.write('\n'.join(bestranking) + '\n')  # Records

    # Copy over any other required files from first batch directory...

    chunk_dir = batches[0].dir

    for file_name in ['gold_protein.mol2']:

        copy(str(chunk_dir / file_name), output_dir)

    ######

    # All done.

    logger.info(f"Finished in {time() - t0:.1f} seconds.")

########################################################################################################################

if __name__ == '__main__':

    main()

########################################################################################################################
# End
########################################################################################################################
