#!/usr/bin/env python
#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2026-06-29: created by Ming-Yu Guo
#

"""Demo: export full disorder-aware CIF data from the CSD Python API.

This script demonstrates a complete workflow for extracting crystal structures
from the Cambridge Structural Database (CSD) with atom-site occupancy,
disorder-group metadata, anisotropic displacement parameters, and bond
connectivity preserved in a manually assembled CIF file.

Why this script exists:
    The built-in CSD Python API CIF writers are convenient, but for some
    disordered entries they do not expose a simple path to preserve the full
    occupancy/disorder atom-site metadata needed for downstream enumeration or
    custom reconstruction workflows.

Pipeline:
    1. Read an entry from the CSD with ``EntryReader``.
    2. Access ``crystal.disordered_molecule`` when available.
    3. Collect atom-site occupancy, disorder assembly/group, ADP, and bonds.
    4. Write a CIF with explicit ``_atom_site_occupancy`` and disorder tags.

Examples:
    python demo_csd_disorder_workflow.py --refcodes ABACIR ABABUB --output-dir ./output
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any


LICENSE_REQUIREMENT = "CSD-Materials or another licence tier with CSD Python API database access"


def _escape_cif_string(value: str | None) -> str:
    return (value or "?").replace("'", "''")


def _format_float(value: Any, precision: int = 4) -> str:
    if value is None:
        return "?"
    return f"{float(value):.{precision}f}"


def _format_esd(value: Any, esd: Any, precision: int = 4) -> str:
    if value is None:
        return "?"
    base = f"{float(value):.{precision}f}"
    if esd is None:
        return base
    try:
        if float(esd) > 0:
            return f"{base}({esd})"
    except (TypeError, ValueError):
        pass
    return base


def _build_disorder_map(crystal: Any) -> dict[str, tuple[int | None, str | None]]:
    disorder_map: dict[str, tuple[int | None, str | None]] = {}
    disorder = getattr(crystal, "disorder", None)
    if disorder is None:
        return disorder_map

    try:
        for assembly in disorder.assemblies:
            assembly_id = str(getattr(assembly, "id", ""))
            for group in assembly.groups:
                group_id = int(getattr(group, "id", 0))
                for atom in group.atoms:
                    disorder_map[str(atom.label)] = (group_id, assembly_id)
    except Exception:
        return disorder_map
    return disorder_map


def _extract_occupancy(atom: Any) -> float:
    try:
        if atom.occupancy is not None:
            return float(atom.occupancy)
    except Exception:
        pass
    return 1.0


def _extract_displacement(atom: Any) -> tuple[str, float | None, dict[str, Any] | None]:
    displacement_type = "?"
    u_iso = None
    aniso = None
    try:
        displacement = atom.displacement_parameters
        if displacement is None:
            return displacement_type, u_iso, aniso
        u_iso = float(displacement.isotropic_equivalent)
        if displacement.type != "Anisotropic":
            return "Uiso", u_iso, aniso
        values = displacement.values
        uncertainties = displacement.uncertainties
        aniso = {
            "u11": values[0][0], "u22": values[1][1], "u33": values[2][2],
            "u12": values[0][1], "u13": values[0][2], "u23": values[1][2],
            "e11": uncertainties[0][0], "e22": uncertainties[1][1], "e33": uncertainties[2][2],
            "e12": uncertainties[0][1], "e13": uncertainties[0][2], "e23": uncertainties[1][2],
        }
        return "Uani", u_iso, aniso
    except Exception:
        return displacement_type, u_iso, aniso


def _extract_atom_data(molecule: Any, disorder_map: dict[str, tuple[int | None, str | None]]) -> list[dict[str, Any]]:
    atom_data: list[dict[str, Any]] = []
    for atom in molecule.atoms:
        label = str(atom.label) if atom.label else "?"
        symbol = str(atom.atomic_symbol) if atom.atomic_symbol else "?"
        fractional = atom.fractional_coordinates
        displacement_type, u_iso, aniso = _extract_displacement(atom)
        disorder_group, disorder_assembly = disorder_map.get(label, (None, None))
        atom_data.append({
            "label": label,
            "symbol": symbol,
            "fx": float(fractional.x) if fractional else None,
            "fy": float(fractional.y) if fractional else None,
            "fz": float(fractional.z) if fractional else None,
            "occupancy": _extract_occupancy(atom),
            "u_iso": u_iso,
            "displacement_type": displacement_type,
            "disorder_group": disorder_group,
            "disorder_assembly": disorder_assembly,
            "aniso": aniso,
        })
    return atom_data


def _extract_bonds(molecule: Any) -> list[tuple[str, str, float, str]]:
    bonds: list[tuple[str, str, float, str]] = []
    try:
        for bond in molecule.bonds:
            atom_1, atom_2 = bond.atoms
            bonds.append((str(atom_1.label), str(atom_2.label), float(bond.length), str(bond.bond_type)))
    except Exception:
        return bonds
    return bonds


def _build_header_lines(refcode: str, entry: Any, crystal: Any) -> list[str]:
    cell_a, cell_b, cell_c = crystal.cell_lengths
    alpha, beta, gamma = crystal.cell_angles
    return [
        f"data_{refcode}",
        "_audit_creation_method            'CSD Python API disorder workflow demo'",
        f"_chemical_name_common             '{_escape_cif_string(entry.chemical_name)}'",
        f"_chemical_formula_sum             '{entry.formula or '?'}'",
        f"_cell_length_a                    {_format_float(cell_a)}",
        f"_cell_length_b                    {_format_float(cell_b)}",
        f"_cell_length_c                    {_format_float(cell_c)}",
        f"_cell_angle_alpha                 {_format_float(alpha, 2)}",
        f"_cell_angle_beta                  {_format_float(beta, 2)}",
        f"_cell_angle_gamma                 {_format_float(gamma, 2)}",
        f"_cell_volume                      {_format_float(crystal.cell_volume, 2)}",
        f"_cell_formula_units_Z             {int(crystal.z_value) if crystal.z_value else '?'}",
        f"_symmetry_space_group_name_H-M    '{crystal.spacegroup_symbol or '?'}'",
        "",
    ]


def _build_symmetry_lines(crystal: Any) -> list[str]:
    if not crystal.symmetry_operators:
        return []
    lines = ["loop_", "_symmetry_equiv_pos_as_xyz"]
    for operator in crystal.symmetry_operators:
        lines.append(f"  '{operator}'")
    lines.append("")
    return lines


def _build_atom_site_lines(atom_data: list[dict[str, Any]]) -> list[str]:
    disorder_present = any(atom["disorder_group"] is not None for atom in atom_data)
    lines = [
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_occupancy",
        "_atom_site_U_iso_or_equiv",
        "_atom_site_thermal_displace_type",
    ]
    if disorder_present:
        lines.extend([
            "_atom_site_disorder_assembly",
            "_atom_site_disorder_group",
        ])

    for atom in atom_data:
        if atom["fx"] is None:
            continue
        row = [
            f"{atom['label']:<8s}",
            f"{atom['symbol']:<4s}",
            f"{atom['fx']:12.6f}",
            f"{atom['fy']:12.6f}",
            f"{atom['fz']:12.6f}",
            f"{atom['occupancy']:8.4f}",
            f"{_format_float(atom['u_iso']):>10s}",
            f"{atom['displacement_type']:>5s}",
        ]
        if disorder_present:
            row.extend([
                f"{str(atom['disorder_assembly']) if atom['disorder_assembly'] is not None else '.':>4s}",
                f"{str(atom['disorder_group']) if atom['disorder_group'] is not None else '.':>4s}",
            ])
        lines.append("  " + " ".join(row))
    lines.append("")
    return lines


def _build_anisotropic_lines(atom_data: list[dict[str, Any]]) -> list[str]:
    anisotropic_atoms = [atom for atom in atom_data if atom["aniso"]]
    if not anisotropic_atoms:
        return []
    lines = [
        "loop_",
        "_atom_site_aniso_label",
        "_atom_site_aniso_U_11",
        "_atom_site_aniso_U_22",
        "_atom_site_aniso_U_33",
        "_atom_site_aniso_U_12",
        "_atom_site_aniso_U_13",
        "_atom_site_aniso_U_23",
    ]
    for atom in anisotropic_atoms:
        u = atom["aniso"]
        lines.append(
            f"  {atom['label']:<8s}"
            f" {_format_esd(u['u11'], u['e11']):>12s}"
            f" {_format_esd(u['u22'], u['e22']):>12s}"
            f" {_format_esd(u['u33'], u['e33']):>12s}"
            f" {_format_esd(u['u12'], u['e12']):>12s}"
            f" {_format_esd(u['u13'], u['e13']):>12s}"
            f" {_format_esd(u['u23'], u['e23']):>12s}"
        )
    lines.append("")
    return lines


def _build_bond_lines(bonds: list[tuple[str, str, float, str]]) -> list[str]:
    if not bonds:
        return []
    lines = [
        "loop_",
        "_geom_bond_atom_site_label_1",
        "_geom_bond_atom_site_label_2",
        "_geom_bond_distance",
        "_ccdc_geom_bond_type",
    ]
    for atom_1, atom_2, distance, bond_type in bonds:
        lines.append(f"  {atom_1:<8s} {atom_2:<8s} {distance:8.4f} {bond_type}")
    lines.append("")
    return lines


def _write_cif_file(output_path: Path, lines: list[str]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")


_SAFE_REFCODE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]{1,19}$")


def _sanitize_refcode(refcode: str) -> str:
    """Ensure refcode contains only safe alphanumeric characters."""
    if not _SAFE_REFCODE_RE.match(refcode):
        raise ValueError(
            f"Invalid refcode {refcode!r}: must be 2-20 alphanumeric characters starting with a letter"
        )
    return refcode


def _safe_output_dir(raw_path: str) -> Path:
    """Resolve and validate output directory, rejecting path traversal."""
    # Only allow relative paths under the current working directory
    if raw_path != Path(raw_path).as_posix().replace("\\", "/"):
        cleaned = Path(raw_path).as_posix()
    else:
        cleaned = raw_path
    if ".." in cleaned.split("/"):
        raise ValueError(f"Path traversal not allowed in output directory: {raw_path!r}")
    resolved = Path.cwd() / cleaned
    resolved = resolved.resolve()
    cwd = Path.cwd().resolve()
    if not str(resolved).startswith(str(cwd)):
        raise ValueError(f"Output directory must be under working directory: {resolved}")
    return resolved


def export_full_cif_from_csd(refcode: str, output_dir: Path) -> dict[str, Any]:
    """Export a CSD structure to CIF with occupancy and disorder metadata."""
    from ccdc.io import EntryReader

    safe_name = _sanitize_refcode(refcode)
    output_path = output_dir / f"{safe_name}.cif"

    with EntryReader("CSD") as reader:
        entry = reader.entry(safe_name)
        crystal = reader.crystal(safe_name)

    molecule = crystal.disordered_molecule or crystal.molecule
    disorder_map = _build_disorder_map(crystal)
    atom_data = _extract_atom_data(molecule, disorder_map)
    bonds = _extract_bonds(molecule)
    lines = []
    lines.extend(_build_header_lines(refcode, entry, crystal))
    lines.extend(_build_symmetry_lines(crystal))
    lines.extend(_build_atom_site_lines(atom_data))
    lines.extend(_build_anisotropic_lines(atom_data))
    lines.extend(_build_bond_lines(bonds))
    lines.append("#END")
    _write_cif_file(output_path, lines)

    return {
        "refcode": refcode,
        "output": str(output_path),
        "n_atoms": len(atom_data),
        "n_partial_occupancy": sum(1 for atom in atom_data if atom["occupancy"] < 0.999),
        "n_bonds": len(bonds),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refcodes", nargs="+", required=True, help="CSD refcodes to export")
    parser.add_argument("--output-dir", default="output", help="Directory for exported CIF files")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = _safe_output_dir(args.output_dir)

    print(f"Licence requirement: {LICENSE_REQUIREMENT}")
    for refcode in args.refcodes:
        summary = export_full_cif_from_csd(refcode, output_dir)
        print(
            f"Exported {summary['refcode']} -> {summary['output']} "
            f"(atoms={summary['n_atoms']}, partial_occ={summary['n_partial_occupancy']}, bonds={summary['n_bonds']})"
        )


if __name__ == "__main__":
    main()
