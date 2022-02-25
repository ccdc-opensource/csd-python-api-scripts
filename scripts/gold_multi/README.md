# GOLD and multiprocessing

## Introduction

This repo contains a script, `gold_multi.py`, which is designed to illustrate how to use the [CSD Docking API](https://downloads.ccdc.cam.ac.uk/documentation/API/descriptive_docs/docking.html) and the standard Python [multiprocessing](https://docs.python.org/3.7/library/multiprocessing.html) module to parallelize GOLD docking. Also included is a simple example system to demonstrate the operation of the script.

Note that you will need to have both [GOLD](https://www.ccdc.cam.ac.uk/solutions/csd-discovery/components/gold/) and the [CSD Python API](https://downloads.ccdc.cam.ac.uk/documentation/API/) installed and licenced in order to use the script.

On a multi-core workstation, this approach should be suitable for docking some hundreds or thousands of ligands depending on the rigour of the docking protocol used; please consult the GOLD USer Guide for information about speed/accuracy tradeoffs in GOLD. Noet that the script is not useful for running GOLD on an HPC compute cluster or on the Cloud: the CCDC provides the GOLD Cluster and GOLD Cloud tools for those use-cases. For further details, please contact [support@ccdc.cam.ac.uk](mailto:support@ccdc.cam.ac.uk).

As ever when using multiprocessing techniques, increasing the number processes will at some point begin to degrade performance as available cores are saturated. At what point this happens will depend on the machine and the workload and thus can only really be determined by experimentation. A default of six was selected as the script was developed on an eight-core workstation and this seemed to give decent performance while leaving cores for other processes.

The script is designed to be as simple as possible in order to not obscure the mechanisms of parallelization. Thus, for example, configuration of the docking is taken entirely from the GOLD conf file. There is also the limitation that only a single input file of ligands is accepted. In addition, the implementation of error-handling and logging is rather lightweight. If a proper application was required then these matters could be addressed.

The script writes output to the directory specified in the GOLD configuration file, and the results can be inspected by loading the GOLD conf file in Hermes as normal (see the Hermes User Guide for details). A `bestranking.lst` file is also written, which records the best-scoring pose for each molecule. Other output normally written by GOLD is not created, although this could be implemented if necessary.

The script partitions the input ligand file into chunks and uses the Docking API and multiprocessing to dock these chunks in parallel using named subdirectories for their output. The solution files for the chunks are then copied to the main output directory and the full `bestranking.lst` file compiled from the partial chunk versions. The intermediate subdirectories are currently kept, but the script could easily be modified to delete them or use anonymous temporary directories if disk usage was to be an issue.

---

## Running the script

To run the script, an environment with the CCDC Python API installed must be active. Further information is available in
the [API installation notes](https://downloads.ccdc.cam.ac.uk/documentation/API/installation_notes.html).

The script is designed to be run from the command line only (and not, for example, from within Hermes).

On Windows, the command would be (in the folder where this archive was unzipped)...

```
> python.exe .\gold_multi.py
```

On Linux or MacOS, an equivalent would be (first making the script executable)...

```
$ chmod u+x ./gold_multi.py

$ ./gold_multi.py
```

In either case, add the option `--help` to show more information.

---

## Note on the input files provided

The example target provided here is SYK tyrosine kinase ([5LMA](https://www.ebi.ac.uk/pdbe/entry/pdb/5lma)).

The ligands in `input.sdf` were built from SMILES. If the name is a PDB code, it means the SMILES corresponded to the crystallographic ligand from that structure (with conventional ionization states assigned). If the name has a suffix, the SMILES is a manually-generated analogue. Note that not all these ligands can be correctly cross-docked into 5LMA, as there are induced-fit effects in SYK that GOLD cannot reproduce.