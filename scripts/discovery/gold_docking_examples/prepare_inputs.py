"""
This script can be used for any purpose without limitation subject to the
conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
This permission notice and the following statement of attribution must be
included in all copies or substantial portions of this script.

"17/11/2026": created by Chris Ringrose for the Cambridge Crystallographic Data Centre
"""

import math
import os
import sys
import urllib.request

from ccdc import io
from ccdc.protein import Protein

PDB_DOWNLOAD_URL = 'https://files.rcsb.org/download/{code}.cif'

# PDB entry -> (output directory, het code of the ligand of interest)
#
# PDB entry -> (output directory, het code, keep_metals)
#
# 1Y6B  VEGFR2 + oxazole inhibitor  Astex Diverse Set member, plain redocking
# 1KE7  CDK2 + LS3                  pharmacophore constraint
# 2VTQ  CDK2 + difluoro analogue    scaffold template for the constraint
# 2VU3  CDK2 + AT7519               the elaborated ligand to be docked
# 6OIM  KRAS G12C + sotorasib       covalent docking reference
# 3HS4  carbonic anhydrase II       metalloprotein docking, catalytic Zn kept

SYSTEMS = {
    '1Y6B': ('normal_docking', 'AAX', False),
    '1KE7': ('pharmacophore_constraint', 'LS3', False),
    '2VTQ': ('scaffold_constraint', 'LZA', False),
    '2VU3': ('scaffold_constraint', 'LZE', False),
    '6OIM': ('covalent_docking', 'MOV', False),
    '3HS4': ('metalloprotein_docking', 'AZM', True),
}


def fetch_mmcif(code, directory):
    """Download a PDB entry as mmCIF, unless it is already present."""
    cif_file = os.path.join(directory, f'{code}.cif')
    if not os.path.exists(cif_file):
        print(f'Fetching {code} from the RCSB PDB...')
        urllib.request.urlretrieve(PDB_DOWNLOAD_URL.format(code=code), cif_file)
    return cif_file


def pick_ligand(protein, het_code, keep_metals):
    """Choose which copy of the ligand to treat as the native one.

    A PDB entry often holds several copies of the same het group, only one of
    which is in the site of interest. For a metalloprotein the right copy is the
    one coordinating the metal; 3HS4, for instance, has three acetazolamides and
    only one of them is in the active site.
    """
    candidates = [l for l in protein.ligands if het_code in l.identifier]
    if not candidates:
        raise RuntimeError(f'No ligand matching {het_code} found')
    if not (keep_metals and protein.metals):
        return candidates[0]

    def distance_to_nearest_metal(ligand):
        return min(
            math.dist(atom.coordinates, metal.coordinates)
            for atom in ligand.atoms
            for metal in protein.metals
        )

    return min(candidates, key=distance_to_nearest_metal)


def prepare(code, directory, het_code, keep_metals=False, force=False):
    """Split a PDB entry into a docking-ready protein and its native ligand."""
    os.makedirs(directory, exist_ok=True)

    protein_file = os.path.join(directory, f'{code}_protein.mol2')
    reference_file = os.path.join(directory, f'{code}_reference.mol2')

    # Never clobber a file someone has prepared by hand: curated protonation or
    # a hand-fixed binding site is exactly the sort of work that is painful to
    # lose and easy to overwrite by rerunning a script.
    existing = [f for f in (protein_file, reference_file) if os.path.exists(f)]
    if existing and not force:
        print(f'  {code}: skipping, these already exist (use --force to replace):')
        for path in existing:
            print(f'    {path}')
        return protein_file, reference_file

    cif_file = fetch_mmcif(code, directory)

    protein = Protein.from_file(cif_file)
    protein.remove_all_waters()
    protein.remove_unknown_atoms()

    # Bond orders for the HET groups are not in the mmCIF file, so let the API
    # perceive them before we pull the ligand out.
    protein.detect_ligand_bonds()

    if keep_metals:
        # detect_ligand_bonds() bonds the metal to whatever it coordinates, so
        # in 3HS4 the zinc and the inhibitor it binds become a single "ligand"
        # and the zinc drops out of protein.metals. Deleting the metal bonds
        # separates them again. GOLD requires this anyway: a metal in the input
        # file must carry no bonds to its coordinating atoms, or it cannot work
        # out the coordination geometry.
        protein.remove_metal_bonds()

    native = pick_ligand(protein, het_code, keep_metals)
    native.add_hydrogens(mode='missing')

    with io.MoleculeWriter(reference_file) as writer:
        writer.write(native)

    # GOLD docks into an empty site, so every ligand comes out of the protein.
    for ligand in protein.ligands:
        protein.remove_ligand(ligand.identifier)
    protein.add_hydrogens()

    with io.MoleculeWriter(protein_file) as writer:
        writer.write(protein)

    print(f'  {code} -> {protein_file}, {reference_file}')
    return protein_file, reference_file


def main():
    arguments = sys.argv[1:]
    force = '--force' in arguments
    wanted = [c.upper() for c in arguments if c != '--force'] or list(SYSTEMS)
    for code in wanted:
        if code not in SYSTEMS:
            print(f'Unknown entry {code}; known entries are {", ".join(SYSTEMS)}')
            continue
        directory, het_code, keep_metals = SYSTEMS[code]
        prepare(code, directory, het_code, keep_metals=keep_metals, force=force)


if __name__ == '__main__':
    main()
