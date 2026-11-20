"""
This script can be used for any purpose without limitation subject to the
conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
This permission notice and the following statement of attribution must be
included in all copies or substantial portions of this script.

"17/11/2026": created by Chris Ringrose for the Cambridge Crystallographic Data Centre
"""

import os
import random
import tarfile
import urllib.request

from ccdc import conformer, io
from ccdc.molecule import Molecule

HERE = os.path.dirname(os.path.abspath(__file__))

ARCHIVE_URL = 'https://drugdesign.unistra.fr/LIT-PCBA/Files/AVE_unbiased.tgz'
ARCHIVE_FILE = os.path.join(HERE, 'AVE_unbiased.tgz')

TARGET = 'MAPK1'          # ERK2: an ATP-competitive site
TEMPLATE = '4zzn'         # one of the PDB templates LIT-PCBA provides

N_ACTIVES = 20
N_INACTIVES = 200
RANDOM_SEED = 1           # fixed so the example is reproducible

LIBRARY_FILE = os.path.join(HERE, 'screening_library.sdf')


def read_smiles(archive, member):
    """Read a LIT-PCBA .smi file: one 'SMILES ID' pair per line."""
    with archive.extractfile(member) as handle:
        for line in handle:
            parts = line.decode().split()
            if len(parts) >= 2:
                yield parts[0], parts[1]


def main():
    if not os.path.exists(ARCHIVE_FILE):
        print('Downloading LIT-PCBA (about 57 MB)...')
        urllib.request.urlretrieve(ARCHIVE_URL, ARCHIVE_FILE)

    archive = tarfile.open(ARCHIVE_FILE)

    # The protein and its native ligand come ready for docking.
    for suffix in ('protein', 'ligand'):
        member = f'{TARGET}/{TEMPLATE}_{suffix}.mol2'
        with archive.extractfile(member) as source:
            with open(os.path.join(HERE, os.path.basename(member)), 'wb') as target:
                target.write(source.read())
        print(f'Extracted {os.path.basename(member)}')

    # Sample from the validation split.
    rng = random.Random(RANDOM_SEED)
    actives = rng.sample(list(read_smiles(archive, f'{TARGET}/active_V.smi')), N_ACTIVES)
    inactives = rng.sample(list(read_smiles(archive, f'{TARGET}/inactive_V.smi')), N_INACTIVES)

    # The identifier carries the known activity, so the poses GOLD writes out
    # can be scored without having to match them back to this file.
    generator = conformer.ConformerGenerator()
    written = 0
    with io.MoleculeWriter(LIBRARY_FILE) as writer:
        for label, entries in (('active', actives), ('inactive', inactives)):
            for smiles, identifier in entries:
                try:
                    molecule = Molecule.from_string(smiles, format='smiles')
                    molecule = generator.generate(molecule)[0].molecule
                except RuntimeError as error:
                    # A few benchmark SMILES fail to embed; skip rather than stop.
                    print(f'Skipping {identifier}: {error}')
                    continue
                molecule.identifier = f'{label}_{identifier}'
                writer.write(molecule)
                written += 1

    print(f'Wrote {written} molecules to {LIBRARY_FILE}')


if __name__ == '__main__':
    main()
