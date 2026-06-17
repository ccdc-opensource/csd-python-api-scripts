#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2022-02-28: Submitted by Abhik Mukhopadhyay , The Cambridge Crystallographic Data Centre
#

import requests
import urllib.request, urllib.parse, urllib.error
import os
import argparse

from ccdc import io
from ccdc import conformer
from ccdc.molecule import Molecule
from ccdc.descriptors import MolecularDescriptors

__doc__ = """
Starting from a list of PDB-codes, the script generates idealized conformers
for ligands and evaluates their RMSD to the conformation in the PDB.
"""

__author__ = 'Brandl, Giangreco, Higueruelo, Schaerfer and Sykes'


excluded_hetids = set(['NAP', 'MA4', 'EOH', 'AGC', 'EOM', 'BGC', 'HEM', 'CPS', 'CPT',
                       'TAU', 'TAR', 'TAS', 'SPK', '1FH', 'SPM', 'VA3', 'SPD',
                       'XPE', 'GD', 'GA', 'EHN', 'CO5', 'LAK', 'V70', 'PE4',
                       'MAC', 'LAT', 'CP2', 'CP3', 'AR', 'MAN', 'CMO', 'TSM',
                       'DTN', 'TSD', 'TSE', 'YBT', 'EDO', 'PCL', 'CB5', 'BEQ',
                       'F09', 'DTV', 'NO2', 'NO3', 'DZZ', 'IUM', 'UVW', 'NME',
                       'ZN', 'NFC', 'NFB', 'NFO', 'K', 'OH', 'NFS', 'NFR',
                       'MQD', 'MG', 'CBG', 'NOE', 'MO', 'MN', 'PC4', 'CBM',
                       'SE', 'CBX', 'MXE', 'BE7', 'BCT', 'DCE', 'W', 'GLO',
                       'MM4', 'GLC', 'DOD', 'SMO', 'GLY', 'DOX', 'FE', 'OF2',
                       'DET', 'C2C', 'FCY', 'MDD', 'AU3', 'FCL', 'FCO', 'HTG',
                       'AUC', 'VER', 'F', 'NFE', 'FMT', 'FMS', 'JEF', 'V',
                       'TOU', 'SX', 'OF3', 'SBT', 'SR', 'SBY', 'LMT', 'LMU',
                       'T3A', 'SM', 'SBE', 'SB', 'SBO', 'MMC', 'YH', 'DPR',
                       'YB', 'CAT', 'DPG', 'DPO', 'CAC', 'TF4', 'CAD', 'NHE',
                       'PGE', 'LA', 'PGO', 'PGL', 'HG2', 'LI', 'LU', 'PGR',
                       'PGQ', 'NHY', '144', 'NMU', 'PG6', 'NH4', 'WCC',
                       'PG5', 'ZRC', 'PG0', 'B7G', 'HGC', 'CFN', 'CFO', 'CFM',
                       'TFH', 'NCP', 'I42', 'TFA', 'PNL', 'NCO', 'Y1', 'CFT',
                       'KO4', 'SF4', 'ZNH', 'ZNO', 'SF3', 'FNE', 'EGL', '7PE',
                       'RU', 'LBT', 'SAL', 'RE', 'RB', 'BF4', 'URE', 'MRY',
                       'CHM', 'ACM', 'MRD', 'IPA', 'IPH', 'ZN3', 'ZN2', 'ZN1',
                       'LI1', '4OA', 'RHD', 'OCM', 'OCL', 'OCO', 'OCN', 'NGN',
                       'HDZ', 'P2K', 'V40', 'HDD', 'HDN', 'MG8', 'BMA', '2FU',
                       'BME', 'WO5', 'WO4', 'OX', 'YT3', 'WO2', 'CEQ', '2FH',
                       'AF3', 'ENC', 'GCD', '2HP', 'B3P', 'CE1', 'EU', 'XCC',
                       'MGF', 'TBR', 'TBU', 'PG4', 'OC8', 'BTC', 'BTB', 'OC1',
                       'OC3', 'OC2', 'OC5', 'OC4', 'OC7', 'OC6', 'NH3', '3OF',
                       'ATO', 'HGI', 'DMF', 'AST', 'MTG', 'SEA', 'SEK', 'GBL',
                       'EU3', 'CLN', 'CLO', 'DUM', 'UNX', 'BEF', 'CLF', 'PDT',
                       'XE', 'PDO', 'TLA', 'CLX', 'DMN', 'FRU', 'BBX', 'BBU',
                       'CLP', 'OMO', 'H2S', 'ZEM', 'KR', 'SRM', 'SE4', 'P6G',
                       '2ME', '2MO', 'CYS', '3PO', 'POR', 'POP', 'PON', '2BM',
                       'CYN', 'GOL', '202', 'F50', 'DHE', 'DHD', 'PO2', 'PO3',
                       'MLP', 'OTE', 'TFP', 'E1H', 'MCR', 'NRU', '3CO', '3CN',
                       '2PA', '2PE', 'NI2', 'NI3', 'NI1', '2PO', '2PN', 'TCN',
                       'HSG', 'HSH', 'TZZ', '12P', '543', 'OF0', 'OF1', 'HSW',
                       'C2O', 'MEE', 'N2O', 'PD', 'VXA', 'GPX', 'THJ', 'BNG',
                       'NIK', 'F2O', 'SYL', 'CCH', 'FDE', 'FDD', 'EMC', 'MO3',
                       'MO2', 'MO1', 'CRY', 'MO7', 'MO6', 'MO5', 'MO4', 'DMS',
                       'DMU', 'DMX', 'RGI', 'AEM', 'PS5', 'PR', 'LDA', 'PT',
                       'PB', 'ZO3', 'TCH', 'PI', 'MH2', 'MH3', 'MHM', 'BU2',
                       'BU3', 'DDQ', 'BU1', 'MHA', 'S', 'DDH', 'T1A', 'LCP',
                       'MOS', 'O4M', 'MOH', 'MOO', 'LCO', 'BCN', 'UMQ', 'NH2',
                       'HNI', 'PEU', 'PER', 'FSO', 'TMA', 'ETA', 'ETF', 'CD1',
                       'PEJ', 'CD3', 'FSX', 'DVT', 'PEG', 'COM', 'RU7', 'CON',
                       'HE5', 'HE6', 'BF2', '16D', '16A', 'P33', 'ND4', 'HEV',
                       'CO', 'CN', 'CM', 'HES', 'CA', 'HEZ', 'CD', 'HEG',
                       'HEA', 'HEB', 'HEC', 'CS', 'CR', 'CU', 'IOD', 'PE8',
                       'PSL', 'OES', 'PE5', 'PE7', 'PE3', 'FS1', 'FS2', 'FS3',
                       'FS4', 'PC3', 'OEC', 'NMO', 'NML', 'DIO', '6MO', 'CIT',
                       'DIA', 'C8E', '1PS', 'C10', 'C15', 'DIS', 'PPF', 'SOG',
                       'PPK', 'SOH', 'PPI', 'IR3', 'SOM', 'TBA', 'PPM', 'SOR',
                       'PPV', 'GAI', 'XL1', 'MBR', 'IR', '3NI', 'IRI', 'SO3',
                       'SO2', 'SO4', 'MTL', 'IN', 'VX', 'PP9', 'MTF', 'CM5',
                       'MTD', 'ARS', 'ART', 'SDS', 'ARF', '4MO', 'BA', 'NAO',
                       'CCN', 'NAW', 'BR', 'EPE', 'O2', 'BOG', 'PIN', 'MYQ',
                       'MYR', 'POL', 'PIG', 'N8E', 'F3S', 'IOH', 'PIS', 'CD5',
                       'BO3', 'BO4', 'DR6', 'NA2', 'PEO', 'NA6', 'CXE', 'NA5',
                       'OS', 'REO', 'U1', 'VO4', 'FLH', 'HO', 'FLO', 'COH',
                       'EEE', 'HG', 'PTN', 'ETI', 'SFO', 'ETN', 'MNR', 'VO3',
                       'MNH', 'OSM', 'ROP', 'SFN', 'MNC', 'PO4', '3BR', 'FMI',
                       'MN3', 'TRE', 'MN5', 'O', 'PBM', 'TRS', 'PT4', 'TRT',
                       'MW3', 'MW2', 'MPD', 'HTO', 'MPO', 'MTO', '6WO', 'MPR',
                       'MPS', 'AU', 'OXY', 'ACA', 'GD3', 'SUC', 'NEH', 'IDT',
                       'SUL', 'P4C', 'ACN', '2OF', 'S0H', 'IDO', 'HFM', 'ACU',
                       'ACT', '2OP', 'ACY', '2OS', 'NI', 'CU1', 'NO', 'NA',
                       'MW1', 'TEE', 'FE2', 'TEP', '1BO', 'FEC', 'FEA', 'FEO',
                       'FEL', 'CUZ', 'FES', 'CUA', 'CUO', 'CUN', 'CUM', 'CUL',
                       'PTL', 'MN6', 'BRP', 'BRO', 'BRM', 'BRJ', 'MES', 'OUT',
                       'AZI', 'MP1', '1PG', '1PE', 'MSM', 'HLT', 'MSF', 'CN1',
                       'CL', 'TL', 'SGM', 'HO3', 'TE', 'TB', 'ICI', 'AG',
                       'FPO', '1CP', 'AL', 'HOH', 'CNF', 'SCN', 'CNB', 'CE',
                       'CNN', 'CNM', 'ICT', '15P', 'ZH3', 'RTC', '2NO', 'DXG',
                       'DXE', 'VSO', 'MSE', 'FAR', 'FII', 'FPP', 'HFP', 'A4P', 'NAD',
                       'ADP', 'AMP', 'APC', 'ATP', 'IMP', 'RVP', 'XMP', 'SAH', 'SAM',
                       'FAD', 'NAP', 'NDP', 'UDP', 'ANP', 'COA', 'CME', 'UMP', 'NAI'])



class SearchAPI:
    def __init__(self, search_url):
        self.search_url = search_url
        self.search_options = '&wt=json&rows=100000'


    def url_response(self, url):
        """
        Getting JSON response from URL
        :param url: String
        :return: JSON
        """
        r = requests.get(url=url)
        # Status code 200 means 'OK'
        if r.status_code == 200:
            json_result = r.json()
            return json_result
        else:
            print(r.status_code, r.reason)
            return None

    def run_search(self):
        """
        Running search with terms
        Check pdbe search api documentation for more details
        :param pdbe_search_term: String
        :return: JSON
        """
        full_query = self.search_url
        response = self.url_response(full_query)
        return response


def get_sdf_hetid(pdbid, hetid, pdb_lig_filename):

    """
    downloads ligand sdf file using RCSB ligand expo API
    """

    hetid = str(hetid)
    sflag = False
    url = "http://ligand-expo.rcsb.org/files/%s/%s/isdf/" % (
        hetid[0], hetid)

    try:
        response = urllib.request.urlopen(url)
    except Exception:
        print(" can't find the ligand for " + str(pdbid))
    else:
        lines = response.readlines()
        sdf_file = ''
        for line in lines:
            try:
                sdf_filename = line.decode().split('<li><a href="')[1].split('"')[0]
                if sdf_filename.startswith(pdbid.lower()) and sdf_filename[-7] != 'D':
                    sdf_url = url + sdf_filename
                    sdf_file = urllib.request.urlopen(sdf_url).read()
                    break
            except Exception:
                continue
        if sdf_file:
            with open(pdb_lig_filename, 'wb') as f:
                f.write(sdf_file)
            sflag = True
        else:
            print(" can't find a ligand for " + str(pdbid) + " or it is disordered")
    return sflag


def get_smi_from_hetid(hetid):
    """
    obtain SMILES string for hetid using PDBe REST API
    """
    smi_query = SearchAPI(
        search_url='https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{}'.format(hetid))
    try:
        smi_call = smi_query.run_search()
    except Exception:
        print("problem accessing the compound API service")
    else:
        for x, y in smi_call.items():
            smi = (y[0]['smiles'][0]['name'])

    return smi


def get_3d_sdf_from_smi(smi, sdf_filename):

    """
    Generate 3D coordinate from SMILES string
    """

    sflag = False
    conformer_generator = conformer.ConformerGenerator()
    conformer_generator.settings.max_conformers = 1
    m = Molecule.from_string(smi)
    conformers = conformer_generator.generate(m)
    with io.EntryWriter(sdf_filename) as writer:
        writer.write(conformers[0].molecule)
    sflag = True

    return sflag


def get_hetids_from_pdbid(pdbid):
    """
    Find hetids in PDB entry using PDBe REST API and then return unique hetid that
    are not present in the excluded_hetids list above
    """
    chem_comp_query = SearchAPI(
        search_url='https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{}'.format(pdbid))
    try:
        cc_call = chem_comp_query.run_search()
    except Exception:
        print("problem accesing the entry API")
    else:
        all_cc_list = []
        for x, y in cc_call.items():
            for z in y:
                cc_id = (z['chem_comp_id'])
                all_cc_list.append(cc_id)
        cc_list = list(set(all_cc_list))

        hetids = [x for x in cc_list if x not in excluded_hetids]

    return hetids


def write_conformers(mol_input, mol_reference, conf_file):
    """
    Generate conformer and calculate rmsd wrt reference
    """
    top_rmsd = -1
    conformer_generator = conformer.ConformerGenerator()
    conformer_generator.settings.superimpose_conformers_onto_reference = True
    conformers = conformer_generator.generate(mol_input)
    if conformers:
        best_rmsd = [
            MolecularDescriptors.rmsd(mol_reference, c.molecule, None,
                                             True, True, True)
            for c in conformers
        ]
        best_rmsd.sort()
        top_rmsd = best_rmsd[0]
        with io.MoleculeWriter(conf_file) as mol_writer:
            for c in conformers:
                mol_writer.write(c.molecule)
            pass
        rmsd, norm_score = round(conformers[0].rmsd(), 2), conformers[
            0].normalised_score
    else:
        rmsd, norm_score = '', ''
    return rmsd, norm_score, top_rmsd


def generate_conformers_for_all_hetids(pdbid, folder_path='conformers_'):
    """
    generate conformer for all hetid
    """
    working_dir = folder_path + pdbid
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)

    list_hetids = get_hetids_from_pdbid(pdbid)
    results = []
    for hetid in list_hetids:
        print('doing hetid ' + str(hetid))
        reference_flag, input_flag = False, False
        pdb_lig_filename = os.path.join(working_dir, str(hetid) + '_reference.sdf')
        model_lig_filename = os.path.join(working_dir, str(hetid) + '_model.sdf')
        reference_flag = get_sdf_hetid(pdbid, hetid, pdb_lig_filename)
        input_smi = get_smi_from_hetid(hetid)
        if input_smi:
            input_flag = get_3d_sdf_from_smi(input_smi, model_lig_filename)
        if reference_flag:
            mol_reader_ref = io.MoleculeReader(pdb_lig_filename)
            try:
                mol_reference = mol_reader_ref[0]
            except Exception:
                continue
        if input_flag:
            mol_reader_input = io.MoleculeReader(model_lig_filename)
            try:
                mol_input = mol_reader_input[0]
            except Exception:
                continue
        count = 0
        if reference_flag and input_flag:
            if len(mol_input.heavy_atoms) == len(mol_reference.heavy_atoms):
                rb = len([b for b in mol_input.bonds if b.is_rotatable])
                if len(mol_input.heavy_atoms) > 6 and rb > 0:
                    # molecule needs hydrogens, set atom types,
                    # and aromatic bonds to match mogul libraries
                    mol_input.add_hydrogens()
                    mol_reference.add_hydrogens()
                    mol_input.standardise_delocalised_bonds()
                    mol_reference.standardise_delocalised_bonds()
                    mol_input.assign_bond_types(which='unknown')
                    mol_reference.assign_bond_types(which='unknown')
                    mol_input.standardise_aromatic_bonds()
                    mol_reference.standardise_aromatic_bonds()
                    mw = int(round(mol_input.molecular_weight, 0))
                    conf_file = model_lig_filename.replace('.sdf', '_conformers.sdf')
                    rmsd, norm_score, top_rmsd = write_conformers(mol_input,
                                                                  mol_reference,
                                                                  conf_file)
                    output_list = ','.join(map(str,
                                               [pdbid, hetid, mw, rmsd,
                                                norm_score, top_rmsd]))
                    results.append(output_list)
    return results


###############################################################################
# main
###############################################################################


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("targets", help="list of PDB-codes", type=str)

    args = parser.parse_args()
    input_filename = args.targets
    base_name = os.path.splitext(input_filename)[0]
    output_filename = base_name + "_CSD_Conformer_Results.csv"

    with open(input_filename) as input_file:
        pdbs_raw = input_file.readlines()

    pdbs = []
    for pdb_raw in pdbs_raw:
        pdb = pdb_raw.strip()
        if pdb and len(pdb) == 4:
            pdbs.append(pdb)
        else:
            print('this line is not a pdbid ' + str(pdb_raw))

    total = str(len(pdbs))
    i = 1

    with open(output_filename, 'w') as output_file:
        output_file.write('pdbid,hetid,mw,rmsd_lowest_energy_conformer,norm_score,lowest_rmsd\n')
        for pdbid in pdbs:
            print('processing ' + str(i) + ' of ' + total + ' PDB: ' + pdbid)
            i += 1
            results = generate_conformers_for_all_hetids(pdbid.strip())
            for result in results:
                output_file.write(result + '\n')


if __name__ == '__main__':
    main()



