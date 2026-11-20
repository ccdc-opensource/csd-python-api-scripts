# Metalloprotein docking — carbonic anhydrase II

Coordinating a catalytic zinc.

System: human carbonic anhydrase II with acetazolamide (PDB 3HS4, 1.1 Å).
Acetazolamide (marketed as Diamox) with carbonic anhydrase is a zinc metalloenzyme: the Zn is held by His94, His96 and His119, and the
inhibitor's sulfonamide nitrogen takes the fourth coordination position.

GOLD works out metal coordination geometry itself. For a zinc it tries
tetrahedral, trigonal bipyramidal and octahedral templates, superimposes them
on the coordinating protein atoms, keeps whichever fits best, and turns the
unoccupied positions into fitting points the ligand can bind to. On this system
it picks **tetrahedral** (RMSD 0.166, against 0.316 and 0.453 for the
alternatives), which is correct. The chosen geometry is reported in
`output/gold_protein.log` and `output/gold_protein.mol2`
contains dummy atoms at the idealised coordination positions.

## Running it

```
python ../prepare_inputs.py 3HS4
python metalloprotein_docking.py
```

## Files

| File | Description |
|---|---|
| `metalloprotein_docking.py` | The example. |
| `3HS4_protein.mol2` | Prepared carbonic anhydrase II coordinates, catalytic Zn retained and unbonded. |
| `3HS4_reference.mol2` | Acetazolamide in its crystal pose, neutral; the RMSD reference. |
| `acetazolamide_anion.mol2` | The deprotonated ligand actually docked, written by the script. |
| `output/` | GOLD's working directory, including `gold_protein.log` and `gold_protein.mol2`. |

RMSD is measured against the crystallographic pose, which still has its original
protonation. GOLD compares heavy atoms, so the missing hydrogen does not affect
the comparison.
