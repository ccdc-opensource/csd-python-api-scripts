# Normal docking — VEGFR2

Redocking a ligand into its own crystal structure.

System: VEGFR2 in complex with a 2-anilino-5-aryl-oxazole inhibitor (PDB 1Y6B),
a member of the Astex Diverse Set, the standard high-quality test set for
validating docking performance (Hartshorn *et al.*, *J. Med. Chem.* **2007**,
*50*, 726).

## Running it

```
python ../prepare_inputs.py 1Y6B
python normal_docking.py
```

## Files

| File | Description |
|---|---|
| `normal_docking.py` | The example. |
| `1Y6B_protein.mol2` | Prepared VEGFR2 coordinates, written by `prepare_inputs.py`. |
| `1Y6B_reference.mol2` | The native ligand in its crystallographic pose; both the docking input and the RMSD reference. |
| `output/` | GOLD's working directory: docked poses, configuration and logs. |
