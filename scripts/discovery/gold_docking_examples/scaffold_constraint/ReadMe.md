# Scaffold constraint docking — CDK2, AT7519 series

Holding a known core in place while GOLD searches the rest.

System: the Astex CDK2 series that produced the clinical candidate AT7519
(Wyatt *et al.*, *J. Med. Chem.* **2008**, *51*, 4986).

| Entry | Ligand | Role |
|---|---|---|
| PDB 2VTQ | the 2,6-difluoro analogue (het code `LZA`) | the template |
| PDB 2VU3 | AT7519 itself, the 2,6-dichloro compound (het code `LZE`) | docked here |

In this example, one member of a congeneric
series has been solved crystallographically, and you want the binding mode of a
close analogue. Docking the analogue freely can slide the shared core out of
the known pose; a scaffold constraint pins the core and lets GOLD search only
the part you are uncertain about.

The constraint substructure-matches the template onto the ligand and adds a
distance constraint per matched atom pair, so it is a spring rather than a weld.
Raise `weight` to hold the core more tightly.

## Running it

```
python ../prepare_inputs.py 2VTQ 2VU3
python scaffold_constraint.py
```

## Two practical points this example has to deal with

* **The template must be the shared core, not the whole template ligand.** The
  maximum common substructure is computed in the script to get the shared core.
* **Bond orders perceived for the two pyrazoles differ between the two PDB
  entries**, so the MCS is run with `check_bond_type` switched off. Leaving it
  on truncates the core at the amide and loses the pyrazole entirely.

2VTQ and 2VU3 are isomorphous, so their coordinates are already in a common
frame (the two ligand centroids sit 0.12 Å apart). For structures that are not,
superpose the proteins first with `ccdc.protein.Protein.ChainSuperposition`.

## Files

| File | Description |
|---|---|
| `scaffold_constraint.py` | The example. |
| `2VU3_protein.mol2` | Prepared CDK2 coordinates for the docking target. |
| `2VTQ_reference.mol2` | The difluoro analogue in its crystal pose; the scaffold template. |
| `2VU3_reference.mol2` | AT7519 in its crystal pose; the docking input and RMSD reference. |
| `output/` | GOLD's working directory: docked poses, configuration and logs. |
