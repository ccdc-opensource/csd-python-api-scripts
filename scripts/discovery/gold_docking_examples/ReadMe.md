# Docking with the CSD Python API

Short, self-contained examples of the main GOLD docking modes, using the
CSD Python API. Each one uses a literature-backed system, runs end to end, and
reports the RMSD of its poses to the crystal structure.

Written using **CSD Python API 3.9.0**. Docking requires a CSD-Discovery or
CSD-Enterprise licence.

## Getting started

Every example reads prepared input files that are fetched from the RCSB PDB and
cleaned up by a shared script. Run it first:

```
python prepare_inputs.py          # all systems
python prepare_inputs.py 1KE7     # or just one
```

That downloads each entry, removes waters, takes the native ligand out of the
binding site, adds hydrogens, and writes a `<code>_protein.mol2` plus a
`<code>_reference.mol2` into the relevant example directory. Each example needs
the same three things, and they come from those two files:

* a protein with waters removed, the native ligand taken out and hydrogens added
* the native ligand on its own, used to define the binding site
* a reference ligand in its crystallographic pose, used to measure RMSD

You can then run whichever example you are interested in, from its own directory:

```
cd normal_docking
python normal_docking.py
```

## The examples

Each directory has its own `ReadMe.md` covering that example in full; the
sections below are the short version.

| Directory | Mode | System |
|---|---|---|
| `normal_docking` | Plain redocking | VEGFR2 (PDB 1Y6B) |
| `pharmacophore_constraint` | Pharmacophore constraint | CDK2 + LS3 (PDB 1KE7) |
| `scaffold_constraint` | Scaffold match constraint | CDK2, AT7519 series (PDB 2VTQ, 2VU3) |
| `covalent_docking` | Covalent, with warhead transformation | KRAS G12C + sotorasib (PDB 6OIM) |
| `metalloprotein_docking` | Coordinating a catalytic metal | Carbonic anhydrase II + acetazolamide (PDB 3HS4) |
| `virtual_screening` | Library screening and enrichment | MAPK1 / ERK2 (PDB 4ZZN) |

### Normal docking

Redocking a ligand into the structure it was solved in. 1Y6B belongs to the Astex Diverse
Set (Hartshorn *et al.*, *J. Med. Chem.* **2007**, *50*, 726), the standard
high-quality test set for docking validation. Success is conventionally an RMSD
of 2 Å or better.

### Pharmacophore constraint

A pharmacophore constraint is a scoring reward: poses
that put an atom of the requested type near a chosen point score better, and
poses that miss it are less competitive. Used here to require the hinge
hydrogen bond that ATP-competitive kinase inhibitors almost always make. In
1KE7, that's Leu83 NH donating to an LS3 oxygen.

### Scaffold constraint

For a congeneric series where one member has been solved and you want the
binding mode of a close analogue. The constraint substructure-matches a template
onto the ligand and adds a distance constraint per matched atom pair, so the
shared core is held on a spring while GOLD searches the rest.

2VTQ and 2VU3 are isomorphous, so no superposition is needed. For structures
that are not, superpose first with `Protein.ChainSuperposition`.

### Metalloprotein docking

Carbonic anhydrase II is a zinc enzyme: the metal is held by His94,
His96 and His119, and acetazolamide's sulfonamide takes the fourth coordination
position. GOLD determines the coordination geometry itself. It superimposes
tetrahedral, trigonal bipyramidal and octahedral templates on the coordinating
protein atoms and keeps the best fit, then turns the unoccupied positions into
fitting points.

### Virtual screening

Redocking asks whether a protocol reproduces a known pose. Screening asks a
different question: does the score rank actives above inactives? The example
docks a small library and reports a ROC curve, its AUC, and an enrichment
factor.

The library is sampled from **LIT-PCBA** (Tran-Nguyen, Jacquemard & Rognan,
*J. Chem. Inf. Model.* **2020**, *60*, 4263), whose inactives are molecules that
were assayed and found to be inactive, rather than property-matched decoys
assumed to be inactive. `prepare_screening_library.py` downloads the AVE-unbiased
release (about 57 MB), takes 20 actives and 200 inactives from the MAPK1
validation split, and generates 3D conformers.

A real screen would typically use the whole set. This sample however is sized to finish quickly on a laptop,
which means the ROC curve is noisy and the enrichment factor is quoted at 5%
rather than the usual 1%.

## Reading the results

Every script prints the fitness and the RMSD to the reference for each pose.
You can set `settings.reference_ligand_file`
and read the `Gold.Reference.RMSD` attribute back off each docked pose:

```python
for pose in docker.results.ligands:
    print(pose.identifier, pose.fitness(), pose.attributes['Gold.Reference.RMSD'])
```

## Plotting alongside the CSD API

There is an additional optional package, MatPlotLib, used to plot the ROC curve in the virtual screening example.
It is imported in the virtual screening example but can be skipped if not installed.
The other examples do not require it.
