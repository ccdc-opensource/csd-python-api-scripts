# Pharmacophore Query Generator

This tool provides the user with the ability to create pharmacophore queries from the results of a ligand overlay.

The pharmacophore queries are produced to be used with CrossMiner.

## Requirements

- [CSD Python API](https://downloads.ccdc.cam.ac.uk/documentation/API/) installed.
- Access to CSD CrossMiner and the feature definitions (`.cpf`) files.
- Access to the CSD Ligand Overlay Tool.

## Licensing Requirements

CSD-Discovery, CSD-Enterprise and Research Partner suites would all be sufficient.

## Instructions on Running

### Feature Definitions

The CrossMiner feature definition (`.cpf`) files are **not** shipped with this repo.
Supply the location of the feature definitions from your CrossMiner installation with
`-f`/`--feature_definitions`; this should be the directory containing the definition files
(either directly, or in `any`/`protein`/`small_molecule` subdirectories).
Usually, the location is `C:\users\<username>\CCDC\ccdc-software\csd-crossminer\feature_definitions`

```
python main.py -i <overlay_folder> -o <output_folder> -f <feature_definitions_folder>
```

The output folder (`-o`/`--output_folder`) is optional; if it is not supplied, the queries are
written to a `queries` folder created in the current directory.

### Options for Ligand Overlay Output

* `cluster`: Cluster the similar pharmacophore features based on proximity
* `projected`: Treat pharmacophore features as projected when appropriate e.g. acceptors and donors

There is also the option to specify a specific Ligand Overlay from all the results. If this is not specified, all
overlays are used.
When all overlays are used, the pharmacophore query will be a union of all the features from all the overlays.
These features are then clustered based on proximity AND prevalence.

### Using the Queries Generated

If you would like to use the queries generated with this tool, they can be opened in CrossMiner to run a search.
A file `crossminer_search.py` has also been provided which contains a Python function for the most simply kind of
CrossMiner search.
