# This file contains imports, configuration and utility functions that are used across multiple notebooks in the CSD Python API Notebooks collection.
# It is used by the notebooks via the '%run' cell magic, and is not designed to be run independently of the notebooks.

import warnings

from platform import platform
import sys
import os
import shutil
from pathlib import Path
import logging
import time

import pandas as pd
import altair

from IPython.display import HTML, SVG, IFrame

# For PyMOL
import subprocess
import xmlrpc.client as xmlrpclib
from time import sleep

# RDKit
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools
from rdkit.Chem.Draw import IPythonConsole

PandasTools.RenderImagesInAllDataFrames()
IPythonConsole.ipython_useSVG = True

# CCDC imports (must be last because of Qt issues)
import ccdc

from ccdc.io import EntryReader, EntryWriter, MoleculeReader, MoleculeWriter
from ccdc.entry import Entry
from ccdc.molecule import Molecule
from ccdc.diagram import DiagramGenerator

########################################################################################################################
# 
# Configuration
#

# Ensure environment variable CSDHOME is set...

assert 'CSDHOME' in os.environ, "Error! Environment variable 'CSDHOME' not set."

##########

# Info useful for debugging...

script_info = f"""
Platform:                     {platform()}

Python exe:                   {sys.executable}
Python version:               {'.'.join(str(x) for x in sys.version_info[:3])}

CSD version:                  {ccdc.io.csd_version()}
CSD directory:                {ccdc.io.csd_directory()}
API version:                  {ccdc.__version__}

CSDHOME:                      {os.environ.get('CSDHOME', 'Not set')}
CCDC_LICENSING_CONFIGURATION: {os.environ.get('CCDC_LICENSING_CONFIGURATION', 'Not set')}

RDKit version:                {rdkit.__version__}
"""

############

# Important! The PyMOL executable named in `pymol_exe` below must be in your path. If you installed PyMOL using `conda`
# this should be the case. If not, you will need to set `pymol_exe` to the name of your PyMOL executable, and to either
# ensure it is in your path or include the full path to the executable in `pymol_exe`.

pymol_exe = 'pymol'

# Utility to start PyMOL (somewhat) robustly...

def start_pymol():
    
    process = subprocess.Popen([pymol_exe, '-R'])

    pymol = xmlrpclib.ServerProxy('http://localhost:9123')

    for n_try in range(10):

        try:
            
            pymol.do('')  # No-op
            
            return pymol
        
        except ConnectionRefusedError as error:
        
            sleep(1)

    return None  # Failed to start PyMOL


############

# Some utility functions to help with display in Jupyter-Lab...


def show_dataframe(df):
    
    return HTML(df.to_html(escape=False).replace(r'\n', ''))


def mol2image(entity, atom_labels=[]):

    molecule = entity.molecule if type(entity) == ccdc.entry.Entry else entity
    
    atom_selection = [x for x in molecule.atoms if x.label in atom_labels]
            
    return diagram_generator.image(molecule, highlight_atoms=atom_selection)
    
    
def mol2html(entity, atom_labels=[]):

    molecule = entity.molecule if type(entity) == ccdc.entry.Entry else entity
    
    atom_selection = [x for x in molecule.atoms if x.label in atom_labels]
            
    return HTML(diagram_generator.image(molecule, highlight_atoms=atom_selection))


##########

# The API does not yet read stereochemical information, so requires a 3D structure as input. Thus, if the structures
# intended as input to the Conformer API are available as SMILES (or, say, a 2D SDF file), we recommend the use of
# [RDKit](http://rdkit.org/) to generate an initial 3D structure which is then used as input to the conformer generator.

# Utility to generate a 3D structure from an input 2D structure using RDKit...

def make_3d(mol, minimization_attempts=0): 
    
    """
    Generate a 3D structure from an input 2D structure using RDKit.
    
    mol is an RDKit molecule object.
    
    Set minimization_attempts to a positive number to enable MMFF minimization. Note that minimisation is no longer
    recommended by RDKit: https://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
    """
    
    # mol = Chem.Mol(mol)  # Copy mol
    
    mol = Chem.AddHs(mol)  # Hs are required for 3D structure generation (also copies mol)
    
    if AllChem.EmbedMolecule(mol) == -1:  # Generate 3D coordinates
    
        print(f"Error! Embedding failed!", file=sys.stderr)
        
        return None
    
    if minimization_attempts:  # Optional MMFF minimization

        for n in range(minimization_attempts):

            if AllChem.MMFFOptimizeMolecule(mol) == 0: break

        else:

            print(f"Warning! Minimisation did not finish after maximum of {minimization_attempts} attempts for '{mol.GetProp('_Name')}'", file=sys.stderr)

    return mol
    
########################################################################################################################   
# 
# Set platform dependent behaviour...
# 

csd_dir = Path(os.environ['CSDHOME']).parent  # CSD System directory

version = csd_dir.name.replace('CSD_', '')  # Release version

discovery_dir = csd_dir.parent / f'Discovery_{version}'  # Corresponding Discovery directory
    
if platform().startswith('Windows'):

    conqest_exe = str(csd_dir / 'ConQuest' / 'exe' / 'csds.exe')
    
    hermes_exe  = str(discovery_dir / 'Hermes' / 'hermes.exe')
    
else:  # Linux and MacOS

    conqest_exe = str(csd_dir / 'bin' / 'cq')
    
    hermes_exe = str(discovery_dir / 'bin' / 'hermes')
    
    os.environ['GOLD_DIR'] = str(discovery_dir / 'GOLD')
    
########################################################################################################################
# 
# Initialization
# 

# Set up a CCDC Diagram Generator...

diagram_generator = DiagramGenerator()

diagram_generator.settings.return_type = 'SVG'
diagram_generator.settings.explicit_polar_hydrogens = False
diagram_generator.settings.shrink_symbols = False

######

# Get logger and configure if necessary...

logger = logging.getLogger(__name__)

if not logger.hasHandlers():
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('[%(asctime)s %(levelname)-7s] %(message)s', datefmt='%y-%m-%d %H:%M:%S'))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

########################################################################################################################
# End
########################################################################################################################