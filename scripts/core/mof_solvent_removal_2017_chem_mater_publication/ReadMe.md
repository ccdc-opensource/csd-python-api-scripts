# MOF solvent removal

## Summary

Scripts included in the supporting information of the article "Development of a Cambridge Structural Database Subset:
A Collection of Metal–Organic Frameworks for Past, Present, and Future", Peyman Z. Moghadam, Aurelia Li,
Seth B. Wiggin, Andi Tao, Andrew G. P. Maloney, Peter A. Wood, Suzanna C. Ward, and David Fairen-Jimenez
*Chem. Mater.* **2017**, 29, 7, 2618–2625, DOI: <https://doi.org/10.1021/acs.chemmater.7b00441>

Scripts are essentially equivalent: one is designed to be run through the Mercury CSD Python API menu to
remove solvent from a single structure present in the visualiser, the second runs from the command line
and takes a list of CSD entries (a .gcd file) to run through the solvent removal process in bulk.

## Requirements

Tested with CSD Python API 3.9.18

## Licensing Requirements

CSD-Core

## Instructions on running

For the script Mercury_MOF_solvent_removal.py:

- In Mercury, pick **CSD Python API** in the top-level menu, then **Options…** in the resulting pull-down menu.
- The Mercury Scripting Configuration control window will be displayed; from the  *Additional Mercury Script Locations*
section, use the **Add Location** button to navigate to a folder location containing the script
- It will then be possible to run the script directly from the CSD Python API menu, with the script running on the structure
shown in the visualiser

For the script Command_prompt_MOF_solvent_removal.py

```cmd
python Command_prompt_MOF_solvent_removal.py <search_results>.gcd
```

```cmd
positional arguments:
  input_file        CSD .gcd file from which to read MOF structures

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_DIRECTORY, --output-directory OUTPUT_DIRECTORY
                         Directory into which to write stripped structures
  -m, --monodentate
                        Whether or not to strip all unidenate (or monodentate) ligands from the structure
  -s SOLVENT_FILE, --solvent-file SOLVENT_FILE
                        The location of a solvent file
```

## Author

*S.B.Wiggin* (2016)

> For feedback or to report any issues please contact [support@ccdc.cam.ac.uk](mailto:support@ccdc.cam.ac.uk)
