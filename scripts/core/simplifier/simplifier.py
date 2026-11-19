#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2024-04-08 Created by Chris Kingsbury, the Cambridge Crystallographic Data Centre
# 2026-10-21 Update to public version (Chris Kingsbury, Alex Moldovan)
# ORCID 0000-0002-4694-5566
#
# ccdc-licence-features not applicable
# flake8: noqa
from pathlib import Path

import base64
import ccdc.diagram
import ccdc.io
import ccdc.molecule
import ccdc.search
import dataclasses
import io
import itertools
import json
import numpy as np
import os
import re
import warnings
from ccdc.utilities import ApplicationInterface
from copy import copy
from dataclasses import dataclass, field
from typing import TypeVar, List, Dict, Union

Molecule = TypeVar('Molecule')
Crystal = TypeVar('Crystal')
Atom = TypeVar('Atom')
Entry = TypeVar('Entry')
Connectivity = TypeVar('Connectivity')


class NothingToSmileAbout(Exception):
    pass


def markdown_to_html(md_text: str) -> str:
    """Convert Markdown text to HTML, using the markdown package if installed or a built-in fallback."""
    try:
        import markdown

        return markdown.markdown(md_text, extensions=["tables", "fenced_code"])
    except ImportError:
        import html

        lines = md_text.splitlines()
        out = []
        in_code = False
        in_list = None
        in_table = False

        def format_inline(text: str) -> str:
            escaped = html.escape(text)
            escaped = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', escaped)
            escaped = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", escaped)
            escaped = re.sub(r"`([^`]+)`", r"<code>\1</code>", escaped)
            return escaped

        def close_blocks():
            nonlocal in_list, in_table
            if in_list:
                out.append(f"</{in_list}>")
                in_list = None
            if in_table:
                out.append("</tbody></table>")
                in_table = False

        for line in lines:
            stripped = line.strip()
            if stripped.startswith("```"):
                close_blocks()
                if not in_code:
                    out.append("<pre><code>")
                    in_code = True
                else:
                    out.append("</code></pre>")
                    in_code = False
                continue

            if in_code:
                out.append(html.escape(line))
                continue

            if not stripped:
                close_blocks()
                continue

            header_match = re.match(r"^(#{1,6})\s+(.*)$", stripped)
            if header_match:
                close_blocks()
                level = len(header_match.group(1))
                title = header_match.group(2)
                slug = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
                out.append(f'<h{level} id="{slug}">{format_inline(title)}</h{level}>')
                continue

            if stripped.startswith("|") and stripped.endswith("|"):
                cells = [c.strip() for c in stripped.strip("|").split("|")]
                if all(re.match(r"^:?-+:?$", c) for c in cells):
                    continue
                if not in_table:
                    close_blocks()
                    out.append("<table><thead><tr>")
                    for cell in cells:
                        out.append(f"<th>{format_inline(cell)}</th>")
                    out.append("</tr></thead><tbody>")
                    in_table = True
                else:
                    out.append("<tr>")
                    for cell in cells:
                        out.append(f"<td>{format_inline(cell)}</td>")
                    out.append("</tr>")
                continue

            ul_match = re.match(r"^[-*]\s+(.*)$", stripped)
            if ul_match:
                if in_table:
                    out.append("</tbody></table>")
                    in_table = False
                if in_list != "ul":
                    if in_list:
                        out.append(f"</{in_list}>")
                    out.append("<ul>")
                    in_list = "ul"
                out.append(f"<li>{format_inline(ul_match.group(1))}</li>")
                continue

            ol_match = re.match(r"^\d+\.\s+(.*)$", stripped)
            if ol_match:
                if in_table:
                    out.append("</tbody></table>")
                    in_table = False
                if in_list != "ol":
                    if in_list:
                        out.append(f"</{in_list}>")
                    out.append("<ol>")
                    in_list = "ol"
                out.append(f"<li>{format_inline(ol_match.group(1))}</li>")
                continue

            close_blocks()
            out.append(f"<p>{format_inline(stripped)}</p>")

        if in_code:
            out.append("</code></pre>")
        close_blocks()
        return "\n".join(out)


class MercuryInterface:

    def __init__(self, leave_unmatched_molecules: bool = False):
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        self.interface = ApplicationInterface(parse_commandline=False)
        self.interface.parse_commandline()
        try:
            self.crystal = self.interface.current_entry.crystal
        except OSError:
            self.crystal = ccdc.io.CrystalReader(self.interface.input_mol2_file).__next__()
            self.crystal.molecule.assign_bonds()
        # This essentially puts the cache in a folder that lives in the user documents along with the ref code
        # folders so it's visible by all runs
        simplifier_cache_folder = Path(self.interface.working_directory_path).parents[2] / "simplifier_cache"
        if not os.path.isdir(simplifier_cache_folder):
            os.mkdir(simplifier_cache_folder)
        simplifier_cache_file = os.path.join(simplifier_cache_folder, "simplifier_data.json")
        self.simplifier = Simplifier(data_path=simplifier_cache_file,
                                     leave_unmatched_molecules=leave_unmatched_molecules)

    def write_output(self):
        self.simplifier.write_output(interface=self.interface)
        self.simplifier.write_output_molecule(file_output=self.interface.output_mol2_file,
                                              crystal=self.crystal)

    def write_readme_if_none_selected(self):
        if len(self.interface.selected_atoms) <= 1:
            readme_md = (Path(__file__).parent / "readme.md").read_text(encoding="utf-8")
            readme_html = markdown_to_html(readme_md)
            Path(self.interface.output_html_file).write_text(readme_html, encoding="utf-8")
            return True
        else:
            return False


@dataclass
class SimplifierData:
    cache_storage: Path
    smiles: Dict = field(default_factory=dict)
    connectivity: Dict = field(default_factory=dict)

    def __post_init__(self):
        if isinstance(self.cache_storage, str):
            self.cache_storage = Path(self.cache_storage)

        self._loaded_from_file = False
        if not self.cache_storage.parent.is_dir():
            raise NotADirectoryError(f'{self.cache_storage} Does not have a valid directory')
        if len(self.smiles.keys()) == 0:
            if self.cache_storage.is_file() and (len(self.cache_storage.read_text()) > 1):
                self.load_from_file()
            else:
                self.write_cache_storage()

    def load_from_file(self):
        if not self.cache_storage.is_file():
            raise FileNotFoundError(f"Could not find the cache file for the simplifier data in {self.cache_storage}")
        with open(self.cache_storage, "r") as f:
            raw_data = json.load(f)
        self.smiles = raw_data['smiles']
        self.connectivity = raw_data['connectivity']

    def write_cache_storage(self):
        with open(self.cache_storage, "w") as file:
            safe_dict = {k: v for k, v in dataclasses.asdict(self).items() if k not in ("cache_storage")}
            file.write(json.dumps(safe_dict, separators=(", \n", ": ")))

    def update_data(self, identifier: str, smiles: str, connectivity: str) -> None:
        self.smiles.update({identifier: smiles})
        self.connectivity.update({identifier: connectivity})
        self.write_cache_storage()


class Simplifier:
    def __init__(self, data_path: Union[str, Path] = None, leave_unmatched_molecules: bool = True):
        self._selected_molecule = None
        self._selected_smiles = None
        self._connectivity = None
        self.simplifier_data = SimplifierData(cache_storage=data_path)
        self.leave_unmatched_molecules = leave_unmatched_molecules

    @property
    def selected_smiles(self) -> str:
        if self._selected_smiles:
            return self._selected_smiles

    @selected_smiles.setter
    def selected_smiles(self, value: str):
        self._selected_smiles = value

    @property
    def connectivity(self) -> str:
        if self._connectivity:
            return self._connectivity

    @connectivity.setter
    def connectivity(self, value: str):
        self._connectivity = value

    @property
    def selected_molecule(self) -> Molecule:
        if self._selected_molecule:
            return self._selected_molecule

    @selected_molecule.setter
    def selected_molecule(self, value: Molecule):
        self._selected_molecule = value

    def get_selected_atoms_and_connectivity(self, user_selected_atoms: List[object]) -> tuple[List[Atom], Connectivity]:
        if not 2 <= len(user_selected_atoms) <= 10:
            raise ValueError("Please select between 2 and 10 atoms")
        self._selected_molecule.assign_bond_types(which="unknown")
        self._selected_molecule = [
            x for x in self._selected_molecule.components if user_selected_atoms[0] in x.atoms
        ][0]
        selected_atoms_labels = [x.label for x in user_selected_atoms]
        selected_atoms = [self._selected_molecule.atom(x) for x in selected_atoms_labels]
        if len(selected_atoms) > 2:
            distances = np.array(
                [
                    [self._selected_molecule.shortest_path(at1, at2) for at2 in selected_atoms]
                    for at1 in selected_atoms
                ]
            )
            selected_atoms, connectivity = self.optimise_association_matrix_with_connectivity(distances, selected_atoms,
                                                                                              iterations=10)
        elif len(selected_atoms) == 2:
            connectivity = '01'
        else:
            raise ValueError("Please select more than one atom")
        return selected_atoms, connectivity

    @staticmethod
    def optimise_association_matrix_with_connectivity(distances: np.ndarray,
                                                      selected_atoms: List[Atom],
                                                      iterations: int = 10,
                                                      branching=True) -> tuple[List[Atom], Connectivity]:
        if not branching:
            for _ in range(iterations):
                total_distance = np.sum(np.diagonal(distances, 1))
                for i, j in itertools.combinations(range(len(selected_atoms)), 2):
                    if np.sum(np.diagonal(swap_two(distances, i, j), 1)) < total_distance:
                        selected_atoms[i], selected_atoms[j] = (
                            selected_atoms[j],
                            selected_atoms[i],
                        )
                        distances = swap_two(distances, i, j)
            connectivity = ''.join([f"{i}{i + 1}" for i in range(len(selected_atoms) - 1)])
        else:
            connectivity = '01' + ''.join(
                [f"{int(np.argmin(row[:row_ix]))}{row_ix}" for row_ix, row in enumerate(distances) if row_ix > 1])

        return selected_atoms, connectivity

    @staticmethod
    def alter_molecule_smiles(selected_atoms: List[Atom], molecule: Molecule) -> str:
        original_atomic_numbers = [x.atomic_number for x in selected_atoms]
        repl_list = ["Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn"]

        for atom, repl in zip(selected_atoms, repl_list):
            atom.atomic_symbol = repl

        smiles_altered = molecule.smiles
        pattern = re.compile("|".join([f"\\[{Xx}?.?.?.?.\\]" for Xx in repl_list]))

        keys = pattern.findall(smiles_altered)

        for key in keys:
            atom_sym = key[1:3]
            if "H" in key:
                extras = key[3:][:-1]
            else:
                extras = ""
            try:
                ind = repl_list.index(atom_sym)
            except ValueError:
                continue
            smiles_altered = smiles_altered.replace(
                key, f"[#{original_atomic_numbers[ind]}{extras}:{ind}]"
            )
        return smiles_altered

    def get_smarts(self, output_crystal: Crystal) -> str:
        for refcode, stored_smiles in self.simplifier_data.smiles.items():
            smarts = ccdc.search.SMARTSSubstructure(stored_smiles)
            self.selected_smiles = None
            if smarts.match_molecule(output_crystal.molecule):
                print(f"Found {refcode}")
                self._selected_smiles = stored_smiles
                self._connectivity = self.simplifier_data.connectivity.get(refcode, None)
                break
        if self.selected_smiles is None:
            raise NothingToSmileAbout("No smiles found in dictionary.")
        return smarts

    def use_simplification(self, crystal: Crystal, smarts: object = None, connectivity: object = None):
        molecule = crystal.molecule
        molecule.assign_bond_types(which="unknown")
        crystal.molecule = molecule

        output_crystal = copy(crystal)
        self._selected_molecule = output_crystal.molecule
        if smarts is None:
            smarts = self.get_smarts(output_crystal)
        elif isinstance(smarts, str):
            self._selected_smiles = smarts
            smarts = ccdc.search.SMARTSSubstructure(smarts)
        elif hasattr(smarts, "smarts"):
            self._selected_smiles = smarts.smarts

        if connectivity is None:
            connectivity = self._connectivity
        else:
            self._connectivity = connectivity
        sss = ccdc.search.SubstructureSearch()
        sss.add_substructure(smarts)

        slabs = [i for i in range(10) if f":{i}" in self._selected_smiles]
        indices = [smarts.label_to_atom_index(str(i)) for i in slabs]

        if not connectivity:
            for i, (a, b) in enumerate(zip(indices[1:], indices)):
                sss.add_vector(f"vec{i}", (0, a), (0, b))
                number_of_bonds = len(indices) - 1
        else:
            for i, (a, b) in enumerate(zip(connectivity[::2], connectivity[1::2])):
                sss.add_vector(f"vec{i}", (0, indices[int(a)]), (0, indices[int(b)]))
                number_of_bonds = len(connectivity[::2])

        hits = sss.search(output_crystal)
        simple_labels = set()

        for hit in hits:
            for i in range(number_of_bonds):
                at_a_l, at_b_l = [x.label for x in hit.vector_atoms(f"vec{i}")]
                at_a, at_b = [self._selected_molecule.atom(x) for x in (at_a_l, at_b_l)]
                simple_labels.update({at_a_l, at_b_l})
                if at_a not in at_b.neighbours:
                    self._selected_molecule.add_bond(1, at_a, at_b)

        if self.leave_unmatched_molecules:
            remove_labels = set()
            for comp in self._selected_molecule.components:
                if len([x for x in comp.atoms if (x.label in simple_labels)]) > 0:
                    remove_labels.update(
                        {x.label for x in comp.atoms if x.label not in simple_labels}
                    )
            self._selected_molecule.remove_atoms([x for x in self._selected_molecule.atoms if x.label in remove_labels])
        else:
            self._selected_molecule.remove_atoms(
                [x for x in self._selected_molecule.atoms if x.label not in simple_labels])

    def define_simplification(self, crystal: Crystal, user_selected_atoms: List[Atom]) -> None:
        identifier = crystal.identifier
        self._selected_molecule = crystal.molecule
        selected_atoms, self._connectivity = self.get_selected_atoms_and_connectivity(user_selected_atoms)

        self.selected_smiles = self.alter_molecule_smiles(selected_atoms, self._selected_molecule)
        print(self._selected_smiles)

        self.simplifier_data.update_data(identifier, self._selected_smiles, self._connectivity)

    def simplification_diagram(self):
        dg = ccdc.diagram.DiagramGenerator()

        smarts = ccdc.search.SMARTSSubstructure(self._selected_smiles)
        slabs = [i for i in range(10) if f":{i}" in self._selected_smiles]
        indices = [smarts.label_to_atom_index(str(i)) for i in slabs]

        cleaned_smiles = self._selected_smiles.replace('~', '')
        mol = ccdc.molecule.Molecule.from_string(cleaned_smiles)

        hats = [x for i, x in enumerate(mol.atoms) if i in indices]
        for i, x in zip(np.argsort(indices), hats):
            x.label = x.label + f":{i}"

        with io.BytesIO() as buf:
            dg.image(mol, hats, label_atoms=hats).save(buf, 'png')
            image_base64 = base64.b64encode(buf.getvalue()).decode('utf-8').replace('\n', '')

        self.simplification_image_html = f'<img id = "simplification_ats" src="data:image/png;base64, {image_base64}" alt="simplification_ats", width = "300">'
        return self.simplification_image_html

    def connectivity_typog(self):
        if not self.connectivity:
            return None
        return ', '.join([f"{a}-{b}" for a, b in zip(self.connectivity[::2], self.connectivity[1::2])])

    def write_output(self, interface: ApplicationInterface, smiles_altered: str = None) -> None:
        if smiles_altered is None:
            smiles_altered = self.selected_smiles
        with interface.html_report(title=f"smiles for {interface.identifier}") as report:
            report.write(
                f'<p id="smiles">{smiles_altered}</p>'
            )
            report.write(
                f'''<p> Connectivity: {self.connectivity_typog()} </p>'''
            )
            report.write(
                self.simplification_diagram()
            )
            report.write(
                f'''<p> Simplifier Cache Storage: <a href="{self.simplifier_data.cache_storage}">{self.simplifier_data.cache_storage}</a> </p>'''
            )

    def write_output_molecule(self, file_output: str, crystal: Crystal) -> None:
        with ccdc.io.CrystalWriter(file_output) as crystal_writer:
            crystal.molecule = self._selected_molecule
            crystal_writer.write(crystal)


# swaps the rows and columns of an array to a new array. Useful for
# substitution of two nodes in an association matrix
def swap_two(arr, i, j):
    A = arr.copy()
    A[[i, j], :] = A[[j, i], :]
    A[:, [i, j]] = A[:, [j, i]]
    return A


# Gets the relevant arguments for the simplifier to the command line or in Prototypes
def get_args():
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "idents",
        type=str,
        nargs="+",
        help="the files/refcodes that will be analysed (or .gcd)",
    )

    parser.add_argument(
        "-s",
        "--smarts",
        type=str,
        default=None,
        help="override SMARTS simplification string - default is first found in simplifier_data.json",
    )

    parser.add_argument(
        "-c",
        "--connectivity",
        type=str,
        default=None,
        help="override connectivity - default is first found in simplifier_data.json",
    )

    parser.add_argument(
        "-od",
        "--output_directory",
        type=str,
        default=None,
        help="put outputs in this directory",
    )
    parser.add_argument(
        "-ff",
        "--file_format",
        type=str,
        default="mol2",
        help="file format for output",
    )

    parser.add_argument(
        "-p",
        "--parameters",
        type=str,
        default=None,
        help="parameters from JSON",
    )

    parser.add_argument("-v", "--verbose", action="store_true", help="print as you go")

    parser.add_argument(
        "-rm",
        "--remove_unmatched",
        action="store_true",
        help="remove unmatched molecules",
    )
    parser.add_argument(
        "-ma",
        "--match_all",
        action="store_true",
        help="match each molecule independently (instead of using the same SMARTS for all)",
    )

    parser.add_argument(
        "-as",
        "--add_smarts",
        action="store_true",
        help="add a smarts string to the json (as ident:smarts) and then quit",
    )
    parser.add_argument(
        '-sd',
        '--simplifier_data',
        default=f'{os.path.join(os.path.abspath(os.path.dirname(__file__)), "simplifier_data.json")}',
        help="Provide the location of the simplifier data cache. Default with CMD is the script location"
    )
    args = parser.parse_args()
    return args


def simplifier_cl(settings=None):
    args = get_args()
    simplify = Simplifier(data_path=args.simplifier_data, leave_unmatched_molecules=not args.remove_unmatched)

    if args.parameters:
        import json

        with open(args.parameters, "rt") as f:
            jsondict = json.load(f)
            args.__dict__.update(jsondict)
    if settings:
        args.__dict__.update(settings)

    if args.add_smarts:
        smarts_val = getattr(args, "smiles", None) or args.smarts
        simplify.selected_smiles = smarts_val
        if args.connectivity:
            simplify.connectivity = args.connectivity
        else:
            simplify.connectivity = "".join([f"{ix}{ix + 1}" for ix in range(smarts_val.count(":") - 1)])
        simplify.simplifier_data.update_data(args.idents[0], simplify._selected_smiles, simplify._connectivity)
        print(f"added {args.idents[0]}: {smarts_val}")
        return None

    if "." not in "".join(args.idents):
        csd_reader = ccdc.io.CrystalReader("CSD")

    if args.idents[0] == "cwd":
        directory = Path(os.getcwd())
        args.idents = [
            str(directory / filenm)
            for filenm in list(directory.iterdir())
            if filenm.suffix in [".cif", ".mol2"]
        ]

    if args.idents[0].endswith(".gcd"):
        csd_reader = ccdc.io.CrystalReader("CSD")
        refcodes = open(args.idents[0], "r").readlines()
        args.idents = [
            x.rstrip("\n") for x in refcodes if len(x.rstrip("\n")) in [6, 8]
        ]

    for index, ident in enumerate(args.idents):
        if args.verbose:
            print(f"{ident} ({index + 1} / {len(args.idents)})")

        if "." in ident:
            source = ccdc.io.CrystalReader(ident)[0]
        else:
            source = csd_reader.crystal(ident.upper())
        if len(source.molecule.atoms) < 1:
            warnings.warn(f"{ident} has no atoms that are viable for this tool")
            continue
        output = copy(source)
        print(source.identifier)
        try:
            if args.smarts:
                simplify.selected_smiles = args.smarts
                smarts = ccdc.search.SMARTSSubstructure(args.smarts)
                if args.connectivity:
                    connectivity = args.connectivity
                else:
                    connectivity = "".join([f"{ix}{ix+1}" for ix in range(simplify.selected_smiles.count(":")-1)])
                simplify.use_simplification(crystal=source, smarts=smarts, connectivity=connectivity)
            else:
                simplify.use_simplification(crystal=source)
        except NothingToSmileAbout:
            warnings.warn(f"{ident} has no smiles stored in cache. Please use mercury and inspect your structure.")
            continue
        filenm = f"{Path(ident).stem}_s.{args.file_format}"
        output_molecule = str(Path(args.output_directory or os.getcwd()) / filenm)
        simplify.write_output_molecule(file_output=output_molecule, crystal=output)

    with open(
        str(Path(args.output_directory or os.getcwd()) / "combined.mol2"), "w"
    ) as outfile:
        for index, ident in enumerate(args.idents):
            filenm = f"{Path(ident).stem}_s.{args.file_format}"
            filepath = Path(args.output_directory or os.getcwd()) / filenm
            if os.path.isfile(filepath):
                with open(filepath, "r") as fob:
                    outfile.write(fob.read())
                    outfile.write("\n\n")


if __name__ == "__main__":
    simplifier_cl()
