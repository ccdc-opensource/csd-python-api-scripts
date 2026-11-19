#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2024-10-30 Created by Alex Moldovan, the Cambridge Crystallographic Data Centre
# ccdc-licence-features not-applicable
#

from pathlib import Path

import ccdc.crystal
import ccdc.io
import ccdc.molecule
import ccdc.search
import json
import numpy as np
import os
import sys
import tempfile
import unittest
import warnings
from unittest.mock import MagicMock, patch

SIMPLIFIER_DIR = Path(__file__).resolve().parent.parent
if str(SIMPLIFIER_DIR) not in sys.path:
    sys.path.insert(0, str(SIMPLIFIER_DIR))

import simplifier as simplifier_mod  # noqa: E402
from simplifier import (  # noqa: E402
    MercuryInterface,
    NothingToSmileAbout,
    Simplifier,
    SimplifierData,
    get_args,
    markdown_to_html,
    simplifier_cl,
    swap_two,
)

warnings.filterwarnings("ignore", category=DeprecationWarning)

ABABEL_MOL2_PATH = Path(__file__).resolve().parent / "ABABEL.mol2"


def load_ababel_crystal(identifier: str = None) -> ccdc.crystal.Crystal:
    """Load the static ABABEL.mol2 crystal fixture."""
    crystal = ccdc.io.CrystalReader(str(ABABEL_MOL2_PATH))[0]
    if identifier is not None:
        crystal.identifier = identifier
    return crystal


def create_multicomponent_ababel_crystal(
        identifier: str = "ABABEL_HYDRATE",
) -> ccdc.crystal.Crystal:
    """Create a multi-component Crystal from ABABEL.mol2 plus a water molecule."""
    crystal = load_ababel_crystal(identifier=identifier)
    mol = crystal.molecule
    water = ccdc.molecule.Molecule.from_string("O", format="smiles")
    for idx, atom in enumerate(water.atoms):
        atom.coordinates = (25.0 + float(idx) * 1.0, 25.0, 25.0)
        atom.label = f"W{atom.atomic_symbol}{idx + 1}"
    mol.add_molecule(water)
    crystal.molecule = mol
    return crystal


class TestSwapTwo(unittest.TestCase):
    """Tests for the swap_two matrix utility function."""

    def test_swap_rows_and_columns(self):
        arr = np.array(
            [
                [0, 1, 5],
                [1, 0, 2],
                [5, 2, 0],
            ]
        )
        swapped = swap_two(arr, 0, 2)
        expected = np.array(
            [
                [0, 2, 5],
                [2, 0, 1],
                [5, 1, 0],
            ]
        )
        np.testing.assert_array_equal(swapped, expected)

    def test_swap_does_not_mutate_original_array(self):
        arr = np.arange(9).reshape(3, 3)
        original_copy = arr.copy()
        _ = swap_two(arr, 0, 1)
        np.testing.assert_array_equal(arr, original_copy)

    def test_swap_same_index_is_identity(self):
        arr = np.arange(16).reshape(4, 4)
        swapped = swap_two(arr, 2, 2)
        np.testing.assert_array_equal(swapped, arr)


class TestSimplifierData(unittest.TestCase):
    """Tests for SimplifierData cache management."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_invalid_parent_directory_raises_not_a_directory(self):
        invalid_cache = self.temp_path / "non_existent_subdir" / "cache.json"
        with self.assertRaises(NotADirectoryError):
            SimplifierData(cache_storage=invalid_cache)

    def test_creates_new_cache_file_when_missing(self):
        cache_file = self.temp_path / "new_cache.json"
        self.assertFalse(cache_file.exists())

        data = SimplifierData(cache_storage=str(cache_file))
        self.assertIsInstance(data.cache_storage, Path)
        self.assertTrue(cache_file.is_file())
        self.assertEqual(data.smiles, {})
        self.assertEqual(data.connectivity, {})

        disk_content = json.loads(cache_file.read_text())
        self.assertEqual(disk_content, {"smiles": {}, "connectivity": {}})

    def test_loads_existing_cache_file_on_init(self):
        cache_file = self.temp_path / "existing_cache.json"
        initial_payload = {
            "smiles": {"HYQUIN08": "[#8H:0]c1ccc([#8H:1])cc1"},
            "connectivity": {"HYQUIN08": "01"},
        }
        cache_file.write_text(json.dumps(initial_payload))

        data = SimplifierData(cache_storage=cache_file)
        self.assertEqual(data.smiles, initial_payload["smiles"])
        self.assertEqual(data.connectivity, initial_payload["connectivity"])

    def test_load_from_file_raises_file_not_found_if_deleted(self):
        cache_file = self.temp_path / "temp_cache.json"
        data = SimplifierData(cache_storage=cache_file)
        cache_file.unlink()

        with self.assertRaises(FileNotFoundError):
            data.load_from_file()

    def test_update_data_persists_to_disk(self):
        cache_file = self.temp_path / "update_cache.json"
        data = SimplifierData(cache_storage=cache_file)

        data.update_data("ACEMID01", "CC([#7H2:1])=[#8:0]", "01")
        self.assertEqual(data.smiles["ACEMID01"], "CC([#7H2:1])=[#8:0]")
        self.assertEqual(data.connectivity["ACEMID01"], "01")

        reloaded = SimplifierData(cache_storage=cache_file)
        self.assertEqual(reloaded.smiles["ACEMID01"], "CC([#7H2:1])=[#8:0]")
        self.assertEqual(reloaded.connectivity["ACEMID01"], "01")


class TestSimplifier(unittest.TestCase):
    """Tests for the core Simplifier class methods."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        self.cache_file = self.temp_path / "simplifier_data.json"
        self.simplifier = Simplifier(data_path=self.cache_file)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_properties_getters_and_setters(self):
        self.assertIsNone(self.simplifier.selected_smiles)
        self.assertIsNone(self.simplifier.connectivity)
        self.assertIsNone(self.simplifier.selected_molecule)

        mol = ccdc.molecule.Molecule.from_string("CCO", format="smiles")
        self.simplifier.selected_smiles = "[#6:0][#6:1][#8:2]"
        self.simplifier.connectivity = "0112"
        self.simplifier.selected_molecule = mol

        self.assertEqual(self.simplifier.selected_smiles, "[#6:0][#6:1][#8:2]")
        self.assertEqual(self.simplifier.connectivity, "0112")
        self.assertEqual(self.simplifier.selected_molecule, mol)

    def test_connectivity_typog(self):
        self.assertIsNone(self.simplifier.connectivity_typog())

        self.simplifier.connectivity = "01"
        self.assertEqual(self.simplifier.connectivity_typog(), "0-1")

        self.simplifier.connectivity = "011223"
        self.assertEqual(self.simplifier.connectivity_typog(), "0-1, 1-2, 2-3")

    def test_optimise_association_matrix_branching_true(self):
        # 4 atoms: atom 2 is closest to atom 0, atom 3 is closest to atom 1
        distances = np.array(
            [
                [0, 2, 1, 4],
                [2, 0, 3, 1],
                [1, 3, 0, 5],
                [4, 1, 5, 0],
            ]
        )
        atoms = ["A0", "A1", "A2", "A3"]
        out_atoms, conn = Simplifier.optimise_association_matrix_with_connectivity(
            distances, atoms, branching=True
        )
        self.assertEqual(out_atoms, ["A0", "A1", "A2", "A3"])
        # row 0 & 1 -> '01'; row 2 min in [:2] is idx 0 -> '02'; row 3 min in [:3] is idx 1 -> '13'
        self.assertEqual(conn, "010213")

    def test_optimise_association_matrix_branching_false(self):
        # Linear chain A - B - C where initial order is [A, C, B] (distances: A-C=4, C-B=2, A-B=2)
        distances = np.array(
            [
                [0, 4, 2],
                [4, 0, 2],
                [2, 2, 0],
            ]
        )
        atoms = ["A", "C", "B"]
        out_atoms, conn = Simplifier.optimise_association_matrix_with_connectivity(
            distances, atoms, iterations=5, branching=False
        )
        self.assertEqual(conn, "0112")
        # Swapping index 1 and 2 reduces super-diagonal sum from (4 + 2 = 6) to (2 + 2 = 4)
        self.assertEqual(out_atoms, ["A", "B", "C"])

    def test_get_selected_atoms_and_connectivity_counts(self):
        crystal = load_ababel_crystal()
        self.simplifier.selected_molecule = crystal.molecule

        # 0 atoms -> ValueError
        with self.assertRaises(ValueError):
            self.simplifier.get_selected_atoms_and_connectivity([])

        # 1 atom -> ValueError (must select more than 1 atom)
        self.simplifier.selected_molecule = load_ababel_crystal().molecule
        with self.assertRaises(ValueError):
            self.simplifier.get_selected_atoms_and_connectivity(
                [self.simplifier.selected_molecule.atoms[0]]
            )

        # 2 atoms -> '01'
        self.simplifier.selected_molecule = load_ababel_crystal().molecule
        atoms_2, conn_2 = self.simplifier.get_selected_atoms_and_connectivity(
            self.simplifier.selected_molecule.atoms[:2]
        )
        self.assertEqual(len(atoms_2), 2)
        self.assertEqual(conn_2, "01")

        # 3 atoms -> 2 pairs (4 characters)
        self.simplifier.selected_molecule = load_ababel_crystal().molecule
        heavy_atoms = [a for a in self.simplifier.selected_molecule.atoms if a.atomic_number > 1]
        atoms_3, conn_3 = self.simplifier.get_selected_atoms_and_connectivity(
            heavy_atoms[:3]
        )
        self.assertEqual(len(atoms_3), 3)
        self.assertEqual(len(conn_3), 4)

    def test_alter_molecule_smiles(self):
        mol = load_ababel_crystal().molecule
        selected = [mol.atom("Cl1"), mol.atom("O1"), mol.atom("N1")]
        self.assertEqual(len(selected), 3)

        altered = Simplifier.alter_molecule_smiles(selected, mol)
        self.assertIn("[#17:0]", altered)
        self.assertIn("[#8:1]", altered)
        self.assertIn("[#7H:2]", altered)

        smarts = ccdc.search.SMARTSSubstructure(altered)
        fresh_mol = load_ababel_crystal().molecule
        self.assertTrue(smarts.match_molecule(fresh_mol))

    def test_define_simplification_updates_cache(self):
        crystal = load_ababel_crystal()
        selected_atoms = [
            crystal.molecule.atom("Cl1"),
            crystal.molecule.atom("O1"),
            crystal.molecule.atom("N1"),
        ]
        self.simplifier.define_simplification(crystal, selected_atoms)

        self.assertIn("ABABEL", self.simplifier.simplifier_data.smiles)
        self.assertIn("ABABEL", self.simplifier.simplifier_data.connectivity)
        self.assertEqual(
            self.simplifier.simplifier_data.smiles["ABABEL"],
            self.simplifier.selected_smiles,
        )

    def test_get_smarts_success_and_nothing_to_smile_about(self):
        crystal = load_ababel_crystal()

        # Empty cache raises NothingToSmileAbout
        with self.assertRaises(NothingToSmileAbout):
            self.simplifier.get_smarts(crystal)

        # Define simplification on a fresh copy so cache is populated for ABABEL
        define_crystal = load_ababel_crystal()
        self.simplifier.define_simplification(
            define_crystal,
            [define_crystal.molecule.atom("Cl1"), define_crystal.molecule.atom("O1")],
        )
        expected_smiles = self.simplifier.selected_smiles

        smarts = self.simplifier.get_smarts(crystal)
        self.assertIsInstance(smarts, ccdc.search.SMARTSSubstructure)
        self.assertEqual(self.simplifier.selected_smiles, expected_smiles)
        self.assertEqual(self.simplifier.connectivity, "01")

    def test_use_simplification_single_component(self):
        define_crystal = load_ababel_crystal()
        self.simplifier.define_simplification(
            define_crystal,
            [define_crystal.molecule.atom("Cl1"), define_crystal.molecule.atom("O1")],
        )

        target_crystal = load_ababel_crystal()
        self.simplifier.use_simplification(target_crystal)
        simplified_mol = self.simplifier.selected_molecule

        # Only the 2 mapped atoms (Cl1 and O1) should remain, connected by a single bond
        self.assertEqual(len(simplified_mol.atoms), 2)
        self.assertEqual({a.atomic_symbol for a in simplified_mol.atoms}, {"Cl", "O"})
        self.assertEqual(len(simplified_mol.bonds), 1)

    def test_use_simplification_with_explicit_smarts_and_no_connectivity(self):
        define_crystal = load_ababel_crystal()
        self.simplifier.define_simplification(
            define_crystal,
            [define_crystal.molecule.atom("Cl1"), define_crystal.molecule.atom("O1")],
        )
        ababel_smarts = self.simplifier.selected_smiles

        target_crystal = load_ababel_crystal()
        fresh_simplifier = Simplifier(data_path=self.cache_file)
        fresh_simplifier.use_simplification(
            target_crystal,
            smarts=ababel_smarts,
            connectivity="",
        )
        simplified_mol = fresh_simplifier.selected_molecule
        self.assertEqual(len(simplified_mol.atoms), 2)
        self.assertEqual(len(simplified_mol.bonds), 1)

    def test_use_simplification_leave_vs_remove_unmatched_molecules(self):
        define_crystal = load_ababel_crystal()
        self.simplifier.define_simplification(
            define_crystal,
            [define_crystal.molecule.atom("Cl1"), define_crystal.molecule.atom("O1")],
        )

        # Multi-component crystal: ABABEL + water (3 atoms)
        multi_crystal_keep = create_multicomponent_ababel_crystal()
        multi_crystal_remove = create_multicomponent_ababel_crystal()

        # Case 1: leave_unmatched_molecules=True (keeps water component = 3 atoms + 2 simplified atoms = 5)
        simplifier_keep = Simplifier(
            data_path=self.cache_file, leave_unmatched_molecules=True
        )
        simplifier_keep.use_simplification(multi_crystal_keep)
        self.assertEqual(len(simplifier_keep.selected_molecule.atoms), 5)

        # Case 2: leave_unmatched_molecules=False (removes water component, leaving only 2 simplified atoms)
        simplifier_remove = Simplifier(
            data_path=self.cache_file, leave_unmatched_molecules=False
        )
        simplifier_remove.use_simplification(multi_crystal_remove)
        self.assertEqual(len(simplifier_remove.selected_molecule.atoms), 2)

    def test_simplification_diagram_returns_base64_img_tag(self):
        define_crystal = load_ababel_crystal()
        self.simplifier.define_simplification(
            define_crystal,
            [define_crystal.molecule.atom("Cl1"), define_crystal.molecule.atom("O1")],
        )
        html_img = self.simplifier.simplification_diagram()
        self.assertTrue(
            html_img.startswith(
                '<img id = "simplification_ats" src="data:image/png;base64, '
            )
        )
        self.assertIn('alt="simplification_ats"', html_img)

    def test_write_output_and_write_output_molecule(self):
        define_crystal = load_ababel_crystal()
        self.simplifier.define_simplification(
            define_crystal,
            [define_crystal.molecule.atom("Cl1"), define_crystal.molecule.atom("O1")],
        )
        expected_smiles = self.simplifier.selected_smiles

        crystal = load_ababel_crystal()
        self.simplifier.use_simplification(crystal)

        # Test write_output_molecule
        out_mol2 = self.temp_path / "output_simplified.mol2"
        self.simplifier.write_output_molecule(str(out_mol2), crystal)
        self.assertTrue(out_mol2.is_file())
        reloaded = ccdc.io.CrystalReader(str(out_mol2))[0]
        self.assertEqual(len(reloaded.molecule.atoms), 2)

        # Test write_output with mocked ApplicationInterface
        mock_interface = MagicMock()
        mock_interface.identifier = "ABABEL"
        written_chunks = []
        mock_report = MagicMock()
        mock_report.write.side_effect = written_chunks.append
        mock_interface.html_report.return_value.__enter__.return_value = mock_report

        self.simplifier.write_output(mock_interface)
        combined_html = "".join(written_chunks)
        self.assertIn(f'<p id="smiles">{expected_smiles}</p>', combined_html)
        self.assertIn("Connectivity: 0-1", combined_html)
        self.assertIn("data:image/png;base64,", combined_html)


class TestSimplifierCLI(unittest.TestCase):
    """Tests for command-line parsing (get_args) and execution (simplifier_cl)."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        self.cache_file = self.temp_path / "simplifier_data.json"
        simplifier = Simplifier(data_path=self.cache_file)
        define_crystal = load_ababel_crystal()
        simplifier.define_simplification(
            define_crystal,
            [define_crystal.molecule.atom("Cl1"), define_crystal.molecule.atom("O1")],
        )
        self.ababel_smarts = simplifier.selected_smiles

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_get_args_parsing(self):
        test_argv = [
            "simplifier.py",
            "structure1.cif",
            "structure2.mol2",
            "-s",
            self.ababel_smarts,
            "-c",
            "01",
            "-od",
            str(self.temp_path),
            "-ff",
            "cif",
            "-v",
            "-rm",
            "-ma",
            "-as",
            "-sd",
            str(self.cache_file),
        ]
        with patch.object(sys, "argv", test_argv):
            args = get_args()

        self.assertEqual(args.idents, ["structure1.cif", "structure2.mol2"])
        self.assertEqual(args.smarts, self.ababel_smarts)
        self.assertEqual(args.connectivity, "01")
        self.assertEqual(args.output_directory, str(self.temp_path))
        self.assertEqual(args.file_format, "cif")
        self.assertTrue(args.verbose)
        self.assertTrue(args.remove_unmatched)
        self.assertTrue(args.match_all)
        self.assertTrue(args.add_smarts)
        self.assertEqual(args.simplifier_data, str(self.cache_file))

    def test_simplifier_cl_add_smarts(self):
        test_argv = [
            "simplifier.py",
            "NEWREF",
            "-as",
            "-s",
            "[#7H2:0]c1ccc([#8H:1])cc1[#6H:2]",
            "-sd",
            str(self.cache_file),
        ]
        with patch.object(sys, "argv", test_argv):
            simplifier_cl()

        data = SimplifierData(cache_storage=self.cache_file)
        self.assertEqual(data.smiles["NEWREF"], "[#7H2:0]c1ccc([#8H:1])cc1[#6H:2]")
        self.assertEqual(data.connectivity["NEWREF"], "0112")

    def test_simplifier_cl_local_mol2_file_and_parameters(self):
        params_path = self.temp_path / "params.json"
        params_path.write_text(json.dumps({"verbose": True}))

        out_dir = self.temp_path / "outputs"
        out_dir.mkdir()

        test_argv = [
            "simplifier.py",
            str(ABABEL_MOL2_PATH),
            "-od",
            str(out_dir),
            "-p",
            str(params_path),
            "-sd",
            str(self.cache_file),
        ]
        with patch.object(sys, "argv", test_argv):
            simplifier_cl()

        simplified_file = out_dir / "ABABEL_s.mol2"
        combined_file = out_dir / "combined.mol2"
        self.assertTrue(simplified_file.is_file())
        self.assertTrue(combined_file.is_file())
        self.assertGreater(combined_file.stat().st_size, 0)

        reloaded = ccdc.io.CrystalReader(str(simplified_file))[0]
        self.assertEqual(len(reloaded.molecule.atoms), 2)

    def test_simplifier_cl_cwd_and_smarts_override(self):
        work_dir = self.temp_path / "cwd_test"
        work_dir.mkdir()
        local_mol2 = work_dir / "sample.mol2"
        local_mol2.write_text(ABABEL_MOL2_PATH.read_text())

        old_cwd = os.getcwd()
        try:
            os.chdir(work_dir)
            test_argv = [
                "simplifier.py",
                "cwd",
                "-s",
                self.ababel_smarts,
                "-sd",
                str(self.cache_file),
            ]
            with patch.object(sys, "argv", test_argv):
                simplifier_cl()

            self.assertTrue((work_dir / "sample_s.mol2").is_file())
            self.assertTrue((work_dir / "combined.mol2").is_file())
        finally:
            os.chdir(old_cwd)

    def test_simplifier_cl_warns_on_unmatched_structure(self):
        empty_cache = self.temp_path / "empty_cache.json"
        SimplifierData(cache_storage=empty_cache)

        out_dir = self.temp_path / "unmatched_out"
        out_dir.mkdir()

        test_argv = [
            "simplifier.py",
            str(ABABEL_MOL2_PATH),
            "-od",
            str(out_dir),
            "-sd",
            str(empty_cache),
        ]
        with patch.object(sys, "argv", test_argv):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                simplifier_cl()
                warning_messages = [str(w.message) for w in caught]
                self.assertTrue(
                    any("has no smiles stored in cache" in m for m in warning_messages)
                )


class TestMercuryInterface(unittest.TestCase):
    """Tests for MercuryInterface with mocked ApplicationInterface."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        # MercuryInterface uses Path(working_directory_path).parents[2] / "simplifier_cache"
        self.working_dir = self.temp_path / "level1" / "level2" / "work_dir"
        self.working_dir.mkdir(parents=True)

    def tearDown(self):
        self.temp_dir.cleanup()

    @patch.object(simplifier_mod, "ApplicationInterface")
    def test_mercury_interface_init_and_write_readme(self, mock_app_interface_cls):
        mock_app = MagicMock()
        mock_app.working_directory_path = str(self.working_dir)
        mock_app.current_entry.crystal = load_ababel_crystal()
        output_html = self.working_dir / "output.html"
        mock_app.output_html_file = str(output_html)
        mock_app.output_mol2_file = str(self.working_dir / "output.mol2")
        mock_app.selected_atoms = []
        mock_app_interface_cls.return_value = mock_app

        mi = MercuryInterface(leave_unmatched_molecules=False)
        expected_cache_dir = self.temp_path / "simplifier_cache"
        self.assertTrue(expected_cache_dir.is_dir())
        self.assertFalse(mi.simplifier.leave_unmatched_molecules)

        # When no atoms are selected, write_readme_if_none_selected renders readme.md to HTML and returns True
        self.assertTrue(mi.write_readme_if_none_selected())
        self.assertTrue(output_html.is_file())
        rendered_html = output_html.read_text(encoding="utf-8")
        expected_html = markdown_to_html(
            (SIMPLIFIER_DIR / "readme.md").read_text(encoding="utf-8")
        )
        self.assertEqual(rendered_html, expected_html)
        self.assertIn("<h1", rendered_html)
        self.assertIn("Simplifier", rendered_html)
        self.assertIn("Setup in Mercury", rendered_html)

        # When 1 atom is selected, write_readme_if_none_selected also returns True
        mock_app.selected_atoms = [mock_app.current_entry.crystal.molecule.atoms[0]]
        self.assertTrue(mi.write_readme_if_none_selected())

        # When >1 atoms are selected, write_readme_if_none_selected returns False
        mock_app.selected_atoms = mock_app.current_entry.crystal.molecule.atoms[:2]
        self.assertFalse(mi.write_readme_if_none_selected())

    @patch.object(simplifier_mod, "ApplicationInterface")
    def test_mercury_interface_write_output(self, mock_app_interface_cls):
        mock_app = MagicMock()
        mock_app.working_directory_path = str(self.working_dir)
        mock_app.current_entry.crystal = load_ababel_crystal()
        mock_app.output_mol2_file = str(self.working_dir / "output.mol2")
        mock_app_interface_cls.return_value = mock_app

        mi = MercuryInterface(leave_unmatched_molecules=True)
        mi.simplifier.write_output = MagicMock()
        mi.simplifier.write_output_molecule = MagicMock()

        mi.write_output()
        mi.simplifier.write_output.assert_called_once_with(interface=mock_app)
        mi.simplifier.write_output_molecule.assert_called_once_with(
            file_output=mock_app.output_mol2_file,
            crystal=mi.crystal,
        )


if __name__ == "__main__":
    unittest.main()
