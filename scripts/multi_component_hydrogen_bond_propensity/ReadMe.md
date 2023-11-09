# Multi-Component Hydrogen Bond Propensity 

## Summary

Performs a multi-component Hydrogen bond propensity calculation for a given library of co-formers
## Example 

Individual reports are generated for each coformer stored in indvidual folders. A summary multicomponent report with rankings is also generated. 

Individual reports include: 

- Predicted intermolecular hydrogen bond propensities 
- Hydrogen bond coordination likelihood

Summary report includes: 
- Chemical Diagram
- Table of ranked components 
- Multi-component hydrogen-bond propensity chart

## Requirements 
- docxtpl
- multi_component_hydrogen_bond_propensity_report.docx
- multi_component_pair_hbp_report.docx
  
## Licensing Requirements 

- CSD-Core

## Instructions on running

```cmd
python multi_component_hydrogen_bond_propensity_report.py HXACAN28
```

```cmd
positional arguments:
  input_structure       Refcode or mol2 file of the component to be screened

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        the working directory for the calculation
  -c COFORMER_LIBRARY, --coformer_library COFORMER_LIBRARY
                        the directory of the desired coformer library
  -f FAILURE_DIRECTORY, --failure_directory FAILURE_DIRECTORY
                        The location where the failures file should be generated
```

The default coformer library is the one supplied with your Mercury install
- for 2023.1 or later, in ```<CSD Install Location>\ccdc-software\mercury\molecular_libraries\ccdc_coformers```
- for 2022.3 or earlier, in ```<CSD Install Location>\Mercury\molecular_libraries\ccdc_coformers```

Ensure the input structure and coformers have the correct bond typing and any charges before running the script

If the calculation fails for any API-coformer combinations, they will be recorded with N/A in the summary report table

## Author
_Andrew Maloney_ (CCDC) 2017

> For feedback or to report any issues please contact [support@ccdc.cam.ac.uk](mailto:support@ccdc.cam.ac.uk)