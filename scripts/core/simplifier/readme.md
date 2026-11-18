# Simplifier

A molecular simplifier, akin to a "ribbon"-type protein representation, designed to run inside the Mercury CSD Python API menu or from the command line.

`Simplifier` takes a set of user-selected atoms in a molecule and creates a reusable SMARTS simplification rule and connectivity skeleton. This rule can then be applied to any crystal structure
containing the same molecular substructure (such as polymorphs, solvates/hydrates, salts, or co-crystals) to replace the full molecule with a simplified skeletal representation.

Information in the publication and guidelines in its supporting information:
"Simplifying crystal structures: a new tool for visualising packing motifs in
molecular crystals as skeletal frameworks"
(DOI: <https://doi.org/10.1039/d6ce00813e>).

## Setup in Mercury

1. Place the `simplifier` folder in a location where you have read/write permissions.
2. Open the latest version of **Mercury**, select **CSD Python API** in the top toolbar, and click **Options...**.
3. Under **Additional Mercury Script Locations**, click **Add Location**, select the parent folder that contains the `simplifier` directory, and click **Save**.
4. In the **CSD Python API** menu in Mercury, you will now see **simplifier** as a submenu containing `define_simplifier.py`, `use_simplifier.py`, and `use_simplifier_remove_unmatched.py`. We
   recommend starting with `define_simplifier.py`.

This has been tested and validated against the 2026.2 version of Mercury and adjacent CSD Python API 3.0.8

## Mercury Workflow

### 1. Defining a Simplification Rule (`define_simplifier.py`)

1. Load a crystal structure in Mercury.
2. Select **between 2 and 10 atoms** belonging to the molecule you want to simplify:
   - **Important**: You must select **more than 1 atom** (at least 2 atoms, up to a maximum of 10 atoms). If 0 or 1 atom is selected when running `define_simplifier.py`, this README will be displayed
     in the Mercury HTML report instead.
   - For **2 selected atoms**, a single connecting edge (`0-1`) is created between them.
   - For **3 to 10 selected atoms**, connectivity is automatically determined using the shortest topological path distances across the molecule in the order the atoms were clicked.
3. Run **CSD Python API > simplifier > define_simplifier.py**.
4. Mercury will:
   - Generate an atom-mapped SMARTS string and connectivity string and store them in the cache (`simplifier_data.json`).
   - Replace the matched molecules in the current view with their simplified skeletal representation.
   - Open an HTML report displaying the SMARTS string, the connectivity pairs, a 2D chemical diagram highlighting the mapped atoms (`:0`, `:1`, ...), and the path to the cache file.

### 2. Using a Simplification Rule (`use_simplifier.py` / `use_simplifier_remove_unmatched.py`)

Once a molecule's simplification rule has been defined and cached, open another structure containing the same substructure (e.g., another polymorph, hydrate, salt, or co-crystal):

- **`use_simplifier.py`**: Replaces matched molecules with their simplified skeleton while **keeping unmatched components** (such as solvent molecules, counterions, or co-formers) intact.
- **`use_simplifier_remove_unmatched.py`**: Replaces matched molecules with their simplified skeleton and **removes all unmatched components** from the output structure.

## Command-Line Interface (CLI) Usage

`simplifier.py` can also be run directly from the command line using the CSD Python API environment:

```bash
python simplifier.py <idents...> [options]
```

### CLI Arguments

- `idents`: One or more CSD refcodes (e.g., `HXACAN28`), structure files
  (`.cif`, `.mol2`), a `.gcd` refcode list file, or `cwd` to process all
  `.cif` and `.mol2` files in the current working directory.
- `-s`, `--smarts`: Override the SMARTS simplification string. Defaults to the
  first matching entry in `simplifier_data.json`.
- `-c`, `--connectivity`: Override the connectivity string (e.g., `011223`).
  Defaults to the cached connectivity or sequential chain connectivity.
- `-od`, `--output_directory`: Directory for simplified output structures
  (`<ident>_s.<ext>` and `combined.mol2`). Defaults to the current directory.
- `-ff`, `--file_format`: Output format for individual structures (default:
  `mol2`).
- `-p`, `--parameters`: Path to a JSON file containing CLI parameter overrides.
- `-v`, `--verbose`: Print progress for each structure processed.
- `-rm`, `--remove_unmatched`: Remove unmatched molecular components from the
  simplified crystal structure.
- `-ma`, `--match_all`: Match each molecule independently rather than using the
  same SMARTS for all.
- `-as`, `--add_smarts`: Add a SMARTS rule (`ident` + `-s <smarts>` and optional
  `-c <connectivity>`) to `simplifier_data.json`, then exit.
- `-sd`, `--simplifier_data`: Path to the `simplifier_data.json` cache file.
  Defaults to that file in the `simplifier` script directory.

### CLI Examples

Simplify a CSD entry using existing cached rules and write outputs to `./out`:

```bash
python simplifier.py HXACAN28 -od ./out -v
```

Simplify all `.cif` and `.mol2` files in the current directory and remove unmatched solvent molecules:

```bash
python simplifier.py cwd -rm
```

Add a custom SMARTS simplification rule and connectivity to the cache file:

```bash
python simplifier.py HYQUIN08 -as -s "[#8H:0]c1ccc([#8H:1])cc1" -c "01"
```

## Cache Storage (`simplifier_data.json`)

- **In Mercury**: `MercuryInterface` stores the cache in the user space `$PROFILE/simplifier_cache/simplifier_data.json` so that rules defined in one Mercury session remain accessible across
  different structures and runs.
- **In the CLI**: Defaults to `<script_location>/simplifier/simplifier_data.json` unless overridden with `-sd` / `--simplifier_data`.
- **Privacy Note**: `simplifier_data.json` contains SMARTS/SMILES strings for structures you have defined. Do not share this file between organisations if it contains proprietary chemical structures.

## Files

- `simplifier.py`: Core library (`Simplifier`, `SimplifierData`, `MercuryInterface`) and command-line interface.
- `define_simplifier.py`: Mercury entry point to define a new simplification rule from selected atoms (2 to 10 atoms).
- `use_simplifier.py`: Mercury entry point to apply a cached simplification rule while retaining unmatched molecules.
- `use_simplifier_remove_unmatched.py`: Mercury entry point to apply a cached simplification rule and delete unmatched molecules.
- `simplifier_data.json`: JSON cache storing `smiles` (SMARTS) and `connectivity` mappings keyed by structure identifier.
- `tests/test_simplifier.py`: Unit test suite covering `SimplifierData`, `Simplifier`, CLI options, and `MercuryInterface`.

## Running Unit Tests

Run the unit test suite using the CSD Python API Python environment:

```bash
python -m unittest discover -s simplifier/tests -v
```

## Author

Chris Kingsbury and Alex Moldovan, CCDC (2024)

For feedback or to report any issues please contact [support@ccdc.cam.ac.uk](mailto:support@ccdc.cam.ac.uk)
