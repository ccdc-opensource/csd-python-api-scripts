# Pharmacophore constraint docking — CDK2

Docking the inhibitor LS3 into CDK2 (PDB 1KE7) with a pharmacophore constraint
on the hinge hydrogen bond.

A pharmacophore constraint is a scoring *reward* rather than a hard restraint:
poses placing an atom of the requested type near the given point score higher.
Use it when you know, from a series of actives or from the biology, that a
particular interaction must be made. Here, that's the hinge hydrogen bond that almost
all ATP-competitive kinase inhibitors form.

Note that it does not forbid poses that miss the point; it makes them less
competitive. With `function_shape='block'` an atom outside the radius simply
contributes nothing, while `'gaussian'` lets the reward decay with distance.

## Running it

```
python ../prepare_inputs.py 1KE7
python pharmacophore_constraint.py
```

## Files

| File | Description |
|---|---|
| `pharmacophore_constraint.py` | The example. |
| `1KE7_protein.mol2` | Prepared CDK2 coordinates, written by `prepare_inputs.py`. |
| `1KE7_reference.mol2` | LS3 in its crystallographic pose; the docking input and the RMSD reference. |
| `1ke7.pdb` | Original protein structure as downloaded from the PDB. |
| `LS3_ideal.sdf` | LS3 coordinates from the PDB chemical component dictionary, in an idealised geometry rather than the bound one. |
| `reference_LS3.mol2` | Earlier copy of the LS3 bound coordinates, kept from the original version of this example. |

## The pharmacophore point

The constraint places an acceptor point at **(-10.910, 52.354, 11.763)**, the
position of the LS3 oxygen that accepts a hydrogen bond from the Leu83 backbone
NH (N···O 2.68 Å), with a 1.5 Å radius and a gaussian function shape.

Place the point where the *ligand atom* should sit, not on the protein partner,
and keep the radius wide enough for a good pose to reach it. A point set a few
ångströms off a real acceptor position with a tight radius will quietly push
docking away from the correct answer rather than towards it.

`use_as_fitting_point=True` also biases where GOLD tries to place the
ligand, rather than only rewarding the pose after the fact.

Valid types: `'acceptor'`, `'donor_acceptor'`, `'donor'`, `'aromatic_center'`
(or `'aromatic_centre'`), `'ring_center'` (or `'ring_centre'`).
