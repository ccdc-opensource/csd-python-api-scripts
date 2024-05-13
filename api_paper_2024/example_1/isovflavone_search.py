#!/usr/bin/env python
#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2024-05-02: created by the Cambridge Crystallographic Data Centre
#

"""
Credit - this code was written by Jason C. Cole, Natalie Johnson and Alex Moldovan
"""
import argparse

from ccdc.io import EntryReader
from ccdc.search import SubstructureSearch, SMARTSSubstructure


class TrihydroxyIsoflavoneHitfinder(SubstructureSearch.HitProcessor):
    """Post process hits - pick out hits that are only tri-substituted
       and then retrieve information from them to tabulate.
    """

    def __init__(self, query_string, database='CSD'):

        self._query = SMARTSSubstructure(query_string)
        self._query_string = query_string
        self._hits = []
        self._extracted_data = []
        self._entry_reader = EntryReader(database)
        self._hits_without_filtering = []
        self._database = database

    @staticmethod
    def _hydroxyl_or_hydroxylate(atom):
        return (atom.atomic_number == 8) and \
            (len(atom.neighbours) == 1 and atom.formal_charge == -1.0) or \
            (len(atom.neighbours) == 2 and (
                    atom.neighbours[0].atomic_number == 1 or atom.neighbours[1].atomic_number == 1))

    def _find_hydroxyls(self, hit):
        return [atom for atom in hit.match_atoms() if self._hydroxyl_or_hydroxylate(atom)]

    def _count_bound_hydroxyls(self, hit):
        return len(self._find_hydroxyls(hit))

    def _get_entry_data_other(self, identifier):
        entry = self._entry_reader.entry(identifier)
        return entry.attributes

    def _get_entry_data_csd(self, identifier):
        entry = self._entry_reader.entry(identifier)

        synonyms = entry.synonyms
        chemical_name = entry.chemical_name
        return {"Chemical Name": chemical_name, "Synonyms": synonyms}

    def _substitution_pattern(self, hit):
        # Get the labelled atoms in the query

        pattern = []
        label_index_lookup = {i: self._query.label_to_atom_index(i) for i in self._query._matches.keys()}
        print(label_index_lookup)
        match_atoms = hit.match_atoms()
        for k in label_index_lookup.keys():

            atom = match_atoms[label_index_lookup[k]]

            if self._hydroxyl_or_hydroxylate(atom):
                pattern.append(str(k))

        return ",".join(sorted(pattern))

    def _get_entry_data(self, hit):
        if self.database == 'CSD':
            d = self._get_entry_data_csd(hit.identifier)
        else:
            d = self._get_entry_data_other(hit.identifier)
        return d | {"Substitution Pattern": self._substitution_pattern(hit)}

    def add_hit(self, hit):
        """
        This is the key method that gets called in the search
        """
        self._hits_without_filtering.append(hit)
        if self._count_bound_hydroxyls(hit) == 3:
            self._hits.append(hit)

    def tabulate(self):
        """
        Generate a dictionary of dictionaries of relevant data from the hits
        """
        data = {}
        for hit in self._hits:
            data[hit.identifier] = self._get_entry_data(hit)

        return data

    def run(self):
        searcher = SubstructureSearch()
        searcher.add_substructure(self._query)
        super().search(searcher, self._entry_reader)


if __name__ == "__main__":

    sub = "$([#1]),$([OH1]),$([OX1H0]),$(O[CH3]),$(Oc1ccccc1)"
    query_string = (f"c(!@[{sub}:1])1c(!@[{sub}:2])c(!@[{sub}:3])c(!@[{sub}:4])c(OC(!@[{sub}:5])"
                    f"=C(c2c(!@[{sub}:6])c(!@[{sub}:7])c(!@[{sub}:8])c(!@[{sub}:9])c(!@[{sub}:10])2)C(=O)c1)")  # noqa

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--database', default='CSD',
                        help='Path to the file to search or "CSD" to use the CSD')
    args = parser.parse_args()

    database = 'tiny.sdf'
    filtered_search = TrihydroxyIsoflavoneHitfinder(query_string, args.database)
    filtered_search.run()

    data = filtered_search.tabulate()
    sorted_ids = sorted(data.keys())
    info_keys = sorted(data[sorted_ids[0]].keys())

    print(f"Without filtering, we have: {len(filtered_search._hits_without_filtering)} hits")
    print(f"After filtering we have: {len(filtered_search._hits)} hits")

    print(",".join(["Identifier"] + info_keys))
    for key in sorted_ids:
        datum = data[key]

        print(",".join([key] + [str(datum[x]) for x in info_keys]))
