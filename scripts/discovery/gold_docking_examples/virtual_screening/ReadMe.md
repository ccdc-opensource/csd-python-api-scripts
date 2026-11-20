# Virtual screening — MAPK1 / ERK2

Docking a library and measuring enrichment.

System: MAPK1 (ERK2), PDB template 4ZZN, with a library sampled from LIT-PCBA.

Redocking tells you whether a protocol reproduces a known pose. Screening asks
a different question: does the score rank actives above inactives? The answer is
summarised by a ROC curve and its area (AUC), where 0.5 is random.

Because actives are rare in practice, the early part of the curve matters more
than the whole of it. The enrichment factor asks how many more actives you find
at the top of the ranked list than chance would give. It is taken at **5%** here
only because the sampled library is small: at 1% the cut would be a couple of
molecules and the number would be mostly noise. On a full screen, report EF at
1%.

## Running it

```
python prepare_screening_library.py
python virtual_screening.py
```

No `prepare_inputs.py` step — this example gets its protein and native ligand
from the LIT-PCBA download instead.

## The library

LIT-PCBA (Tran-Nguyen, Jacquemard & Rognan, *J. Chem. Inf. Model.* **2020**,
*60*, 4263) is a virtual screening benchmark drawn from confirmatory PubChem
dose-response assays. Unlike older sets built from property-matched decoys, its
inactives are molecules that were actually tested and found inactive, which
avoids rewarding a method for simply telling actives and decoys apart on gross
properties.

`prepare_screening_library.py` downloads the AVE-unbiased release (about 57 MB),
takes the MAPK1 target, samples 20 actives and 200 inactives from its validation
split with a fixed random seed, and generates 3D conformers. It also extracts
the GOLD-ready protein and ligand files that LIT-PCBA ships for the 4ZZN
template. The identifier of each molecule carries its known activity
(`active_…` / `inactive_…`), so the poses GOLD writes out can be scored without
matching them back to the library file.

A real screen uses every molecule. The sample here is sized to finish on a
laptop, which means the resulting ROC curve is noisy.

## Files

| File | Description |
|---|---|
| `prepare_screening_library.py` | Downloads LIT-PCBA and builds the library. |
| `virtual_screening.py` | The example: docks the library, reports AUC and enrichment, plots the ROC curve. |
| `AVE_unbiased.tgz` | The LIT-PCBA download. |
| `4zzn_protein.mol2`, `4zzn_ligand.mol2` | Extracted from LIT-PCBA; the target and the ligand defining its site. |
| `screening_library.sdf` | The sampled library with 3D conformers. |
| `roc_curve.png` | Written by `virtual_screening.py`. |
| `output/` | GOLD's working directory: docked poses, configuration and logs. |
