# Find Binding Conformation

## Summary

We know that most pharmaceutically relevant compounds bind to their targets in a relaxed conformation. The challenge in discovery is to figure out rapidly which conformations are readily accessible for the molecules we are considering. There is now a new solution to address this based on statistical, rather than just energetic approaches.

Driven by the wealth and diversity of bond, angle and torsion information in the Cambridge Structural Database (CSD), the CSD Conformer Generator produces realistic ensembles of low energy ligand structures. These are ready to be exploited for drug design in the presence and also in the absence of detailed knowledge about the three-dimensional structure of the protein active site.

Starting from a list of PDB-codes, this script generates idealized conformers
for ligands and evaluates their RMSD to the conformation in the PDB.

The output  are subdirectories for each PDB entry with the conformers generated for each ligand, and a spreadsheet (.csv) with the results of the comparison.
## Requirements
- Tested with CSD Python API 3.0.9 
- This script uses PDBe's and RCSB's API to obtain PDB related information.

## Licensing Requirements 
- CSD-Core

## Instructions on Running

```cmd
python find_binding_conformation.py pdb_example.txt
```

## Author
_'Brandl, Giangreco, Higueruelo, Schaerfer and Sykes'_



> For feedback or to report any issues please contact [support@ccdc.cam.ac.uk](mailto:support@ccdc.cam.ac.uk)