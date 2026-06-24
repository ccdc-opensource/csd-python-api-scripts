# Working with morphologies and CSD-Particle

This is a collection of scripts from the November 2023 CCDC Webinar on crystal morphologies, presented by Dr Andrew G.
P. Maloney.

## Contents

This folder contains the following scripts:

- `calculate_morphologies_tabulate_input.py`
  Runs a calculation of BFDH morphology and VisualHabit morphology for `IBPRAC18` (ibuprofen), printing a table of
  facet properties.
- `exploring_surface_properties.py`
  Calculates the VisualHabit morphology of `IBPRAC18` (ibuprofen) and descriptors for the (100) surface. Also calculates
  particle rugosity.
- `morphology_plot.py`
  Contains a function to generate a 3D plot of a crystal morphology using Matplotlib. Running the script will plot the
  BFDH morphology of `IBPRAC18` (ibuprofen).

## Requirements

All of the modules required to run these scripts are available in the miniconda environment installed with the CSD
Python API.

## Licensing requirements

- CSD-Core, CSD-Particle

## Usage

Each script can be run independently from the command line:

`python script_name.py`

### Author

Pietro Sacchi

