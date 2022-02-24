# Concat Mol2

----

## Summary

Opens a set of mol2 files in the working directory and creates one concatenated multi-mol2 file.
Optionally delete the individual mol2 files other than the new concat.mol2

## Requirements

CSD Python API not required.
Mol2 files must be in the same directory.

## Licensing Requirements

No licence required

## Instructions on running

```cmd
> python concat_mol2.py
```

Help output:
```cmd
> python concat_mol2.py -h

usage: concat_mol2.py [-h] [-d]

optional arguments:
  -h, --help            show this help message and exit
  -d, --delete_contributors
                        Remove contributing individual mol2 files after
                        concatenation

```

## Author

_

> For feedback or to report any issues please contact [support@ccdc.cam.ac.uk](support@ccdc.cam.ac.uk)