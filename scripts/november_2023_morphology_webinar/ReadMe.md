# Scripts

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
  Contains a function to generate a 3D plot of a crystal morphology using Matplotlib.
- `morphology_shape_classification`
  Contains a definition of the `ShapeClassifier` class, which is used to calculate particle shape. This script also
  contains functions to generate Zingg plots.
- `shape_classification.py`
  Runs an example of the usage of the `ShapeClassifier` class using the structure `SUCACB02`.

## Requirements

All of the modules required to run these scripts are available in the miniconda environment installed with the CSD
Python API.

## Licensing requirements

- CSD-Core, CSD-Particle

### Author

Pietro Sacchi

