#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2016-03-22: created by Anthony Reilly, The Cambridge Crystallographic Data Centre
# 2016-12-06: updated by Anthony Reilly, The Cambridge Crystallographic Data Centre
#

"""Packing_Similarity_Dendrogram.py - Construct a dendrogram for an input set of structures based on packing-similarity
analysis
"""

import sys
import argparse
import matplotlib

matplotlib.use('Agg')
from ccdc.io import EntryReader, CrystalWriter, MoleculeWriter
from ccdc.crystal import PackingSimilarity
import numpy as np

import matplotlib.pyplot as plt
import os


def strip_terminal(name, reader):
    """
    Iteratively removes terminal atoms and chains from molecules, stopping at atoms in the allowed list or in rings
    Saves the resulting structures to "stripped.cif" for use in code
    """
    allowed_terminal_atoms = ["S", "O", "N"]
    xtal_writer = CrystalWriter(name + "_stripped.cif")

    for entry in reader:
        crystal = entry.crystal
        crystal.assign_bonds()
        molecule = crystal.molecule
        # Zealously strip atoms with only one bond, unless allowed by filter above (e.g. C=O)
        atom_diff = -1
        while atom_diff != 0:
            n_atoms_before = len(molecule.atoms)

            for atom in molecule.atoms:
                if len(atom.bonds) == 1 and atom.atomic_symbol not in allowed_terminal_atoms:
                    molecule.remove_atom(atom)
            n_atoms_after = len(molecule.atoms)
            atom_diff = n_atoms_after - n_atoms_before

        crystal.molecule = molecule
        xtal_writer.write(crystal)


# Cluster functions
def compare_clusters(c1, c2, relations, mode):
    # Set sensible bounds on the cluster levels
    if mode == 'complete':
        level = 10000
        for id1 in c1['identifiers']:
            for id2 in c2['identifiers']:
                l = relations[id1][id2]
                level = min(level, l)
        return level
    elif mode == 'single':
        level = -1
        for id1 in c1['identifiers']:
            for id2 in c2['identifiers']:
                l = relations[id1][id2]
                level = max(level, l)
        return level
    elif mode == 'average':
        l = []
        for id1 in c1['identifiers']:
            for id2 in c2['identifiers']:
                l.append(relations[id1][id2])
        level = sum(l) / len(l)
        return level


def merge_equal_levels(cluster):
    children = cluster['children']
    for child in children:
        merge_equal_levels(child)
    pos = 0
    while pos < len(children):
        child = children[pos]
        if child['level'] == cluster['level']:
            children.remove(child)
            children.extend(child['children'])
        else:
            pos += 1


def merge_clusters(c1, c2, level):
    new_cluster = {'level': level, 'identifiers': c1['identifiers'] + c2['identifiers'], 'children': [c1, c2]}
    return new_cluster


def plot_dendrogram(cluster_list, n_ps_mols, filename):
    """
    Function for producing a dendrogram from an input cluster hierarchy
    """

    def get_terminal_count(cluster, total):
        """
        Find the number of terminals in the current cluster
        """
        if not cluster['children']:
            total += 1
        else:
            for child in cluster['children']:
                total = get_terminal_count(child, total)
        return total

    def get_terminals(cluster, terminals, hs):
        """
        Find the heights of the cluster terminals
        """
        if not cluster['children']:
            terminals.append(float(hs[str(cluster['identifiers'])]))
        else:
            for child in cluster['children']:
                terminals = get_terminals(child, terminals, hs)
        return terminals

    def assign_y_positions(cluster, c, terminal_y_positions):
        """
        Assign heights to each of the terminal positions in the current cluster
        """
        if not cluster['children']:
            c += 1
            terminal_y_positions.update({str(cluster['identifiers']): c})
        else:
            for child in cluster['children']:
                terminal_y_positions, c = assign_y_positions(child, c, terminal_y_positions)
        return terminal_y_positions, c

    def get_midpoint(cluster, hs):
        """
        Find the midpoint of the current cluster
        """
        terminal_list = []
        terminal_list = get_terminals(cluster, terminal_list, hs)
        midpoint = sum(terminal_list) / len(terminal_list)
        return midpoint

    def plot_tree(cluster, x_start_positions, y_start_positions, hs):
        """
        Recursive function to plot a dendrogram
        """
        node = cluster['level']

        if cluster['children']:
            terminal_list = []
            terminal_list = get_terminals(cluster, terminal_list, hs)
            midpoint = sum(terminal_list) / len(terminal_list)

            xpos = [x_start_positions[1], node]
            ypos = [y_start_positions[1], midpoint]
        else:
            if len(cluster['identifiers']) != 1:
                # Starting groups that are merged start at the highest level to indicate this
                xpos = [x_start_positions[1], node]
            else:
                # Point joins at the next level
                xpos = [x_start_positions[1], x_start_positions[1] + 1]
            ypos = [y_start_positions[1], hs[str(cluster['identifiers'])]]

        if node != 0:
            # Plotting x across and y up or down
            plt.plot(xpos, [ypos[1], ypos[1]], marker=None, linestyle='-', linewidth=1.0,
                     color="Black", zorder=2)
            plt.plot([xpos[0], xpos[0]], ypos, marker=None, linestyle='-', linewidth=1.0,
                     color="Black", zorder=2)

        if cluster['children']:
            # Now repeat for each child
            for child in cluster['children']:
                plot_tree(child, xpos, ypos, hs)
        else:
            # Plot terminal point and add structure indices
            plt.scatter(xpos[1], ypos[1], color="Blue", zorder=3, s=10)
            plt.annotate(str(",".join(cluster['identifiers'])), xy=(xpos[1], ypos[1]),
                         xytext=(xpos[1] + 0.15, ypos[1]), verticalalignment='center', fontsize='3')

    # Setup tree plotting by getting each terminal's height
    heights, count = assign_y_positions(cluster_list[0], 0, {})

    # Set start of the tree - middle of the plot and 1,1
    xpositions = [1, 1]
    ypositions = [1, get_midpoint(cluster_list[0], heights)]

    # Plot tree
    plot_tree(cluster_list[0], xpositions, ypositions, heights)

    # Plot formatting
    ax = plt.axes()
    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    levels = range(n_ps_mols, -1, -1)
    highlighted_levels = range(n_ps_mols, -1, -1)
    #    highlighted_levels.append(0)

    for level in highlighted_levels:
        plt.plot([level, level], [0, count], "--", linewidth=0.5, color="Gray", zorder=1)

    # Pad the plot to have enough space for structure indices
    plt.xlim(-1, n_ps_mols + 5)
    plt.ylim(0, count + 1)
    ax.set_xticks(levels)
    ax.tick_params(axis='x', bottom='off', top='off')
    plt.xlabel('Packing Similarity / ' + str(n_ps_mols) + ' Molecules', fontsize='large')
    # Save output
    plt.savefig(filename + "_packing_similarity_tree.png", dpi=1000, bbox_inches='tight')
    print("Packing tree diagram saved to " + filename + "_packing_similarity_tree.png")


def main(input_file, matrix_file, n_ps_mols, output_ps_results, conf_threshold, ps_angles, ps_distances, strip,
         n_struct, allow_mol_diff, cluster_mode):
    # Initialise Packing Similarity
    ps = PackingSimilarity()
    ps.settings.ignore_hydrogen_positions = True
    ps.settings.ignore_bond_types = True
    ps.settings.match_entire_packing_shell = False
    # Deal with e.g. salt forms
    if allow_mol_diff:
        ps.settings.allow_molecular_differences = True
        ps.settings.ignore_hydrogen_counts = True
        ps.settings.ignore_bond_counts = True
    else:
        ps.settings.allow_molecular_differences = False
        ps.settings.ignore_hydrogen_counts = False
        ps.settings.ignore_bond_counts = False

    # Deal with solvates
    ps.settings.ignore_smallest_components = True
    ps.settings.packing_shell_size = n_ps_mols
    ps.settings.angle_tolerance = ps_angles
    ps.settings.distance_tolerance = ps_distances

    refcodes = []

    input_name = input_file.rsplit(".")[0]
    print("--------------------------------------------------------")

    if not matrix_file:
        # Read in the input file
        print("Reading input database/gcd/structures:", input_file)

        if not strip:
            structure_reader = EntryReader(input_file)
        else:
            structure_reader = EntryReader(input_file)
            strip_terminal(input_name, structure_reader)
            structure_reader.close()
            structure_reader = EntryReader(input_name + "_stripped.cif")

        if n_struct:
            # noinspection PyTypeChecker
            structure_size = min(len(structure_reader), n_struct)
        else:
            # noinspection PyTypeChecker
            structure_size = len(structure_reader)

        # Initialise matrix
        matrix = np.zeros((structure_size, structure_size))

        print("Generating matrix of packing similarities")

        overlay_folder = input_name + "_overlays"
        if output_ps_results:
            g = open("packing_similarity_results.txt", "w")
            g.write("Packing Similarity Analysis for: " + input_file + "\n")
            if not os.path.exists(overlay_folder):
                os.makedirs(overlay_folder)

        for i in range(0, structure_size):
            refcodes.append(str(structure_reader[i].identifier))

        for i in range(0, structure_size):
            entry_i = structure_reader[i]
            crystal_i = entry_i.crystal
            refcodes.append(str(structure_reader[i].identifier))

            for j in range(i, structure_size):
                if i == j:
                    matrix[i, i] = ps.settings.packing_shell_size
                    continue
                entry_j = structure_reader[j]
                crystal_j = entry_j.crystal

                result = ps.compare(crystal_i, crystal_j)

                if result is not None:
                    if result.nmatched_molecules != 1:
                        matrix[i, j] = result.nmatched_molecules
                        matrix[j, i] = result.nmatched_molecules
                    # For single-molecule matches enforce a threshold for conformation RMSDs
                    elif result.rmsd < conf_threshold:
                        matrix[i, j] = result.nmatched_molecules
                        matrix[j, i] = result.nmatched_molecules
                    else:
                        matrix[i, j] = 0
                        matrix[j, i] = 0
                    if output_ps_results:
                        g.write(entry_i.identifier + " " + entry_j.identifier + ": " + str(int(matrix[i, j])) +
                                " molecules \n")

                        overlay_writer = MoleculeWriter(os.path.join(overlay_folder,
                                                                     "overlay_{0}_{1}.mol2".format(refcodes[i],
                                                                                                   refcodes[j])))
                        mols = result.overlay_molecules()
                        for mol in mols:
                            overlay_writer.write(mol)

                else:
                    matrix[i, j] = 0
                    matrix[j, i] = 0
                    if output_ps_results:
                        g.write(entry_i.identifier + " " + entry_j.identifier + ": " + str(int(matrix[i, j])) +
                                " molecules \n")

        np.savetxt(input_name + '_similarity_matrix.txt', matrix, delimiter=',')
        print("Packing similarity matrix saved to similarity_matrix.txt")
        if output_ps_results:
            g.close()
    else:
        # Use the input matrix
        print("Reading input matrix:", matrix_file)
        matrix = np.loadtxt(matrix_file, delimiter=',')
        # Check the matrix is consistent with options specified and the input structures
        if int(matrix[1, 1]) != n_ps_mols:
            print("Error - input matrix and requested packing-shell size do not match")
            sys.exit(1)
        structure_reader = EntryReader(input_file)
        # noinspection PyTypeChecker
        structure_size = len(structure_reader)
        for i in range(0, structure_size):
            refcodes.append(str(structure_reader[i].identifier))

        structure_reader.close()
        if len(matrix) != structure_size:
            print("Error - input matrix does not contain the same number of structures as inputted")
            sys.exit(1)

    print("--------------------------------------------------------")

    # Set structure label based on location in input file
    labels = []
    for i in range(0, structure_size):
        labels.append(str(i + 1))

    # Generate a look-up table of results from matrix
    relations = {}
    for i in range(0, structure_size):
        relations.update(
            {refcodes[i]: {refcodes[x]: matrix[i, x].tolist() for x in range(0, structure_size)}})

    # Generate a cluster hierarchy - populate it initially with every structure
    cluster_list = []
    for i in range(0, structure_size):
        cluster_list.append({'level': n_ps_mols, 'identifiers': [refcodes[i]], 'children': []})

    # Merge the structures - getting best match for each structure
    while len(cluster_list) > 1:
        best_score = -1
        best_cluster1 = None
        best_cluster2 = None
        for i in cluster_list:
            for j in cluster_list:
                if i == j:
                    continue
                score = compare_clusters(i, j, relations, cluster_mode)
                if score > best_score:
                    best_score = score
                    best_cluster1 = i
                    best_cluster2 = j
        new_cluster = merge_clusters(best_cluster1, best_cluster2, best_score)
        cluster_list.remove(best_cluster1)
        cluster_list.remove(best_cluster2)
        cluster_list.append(new_cluster)

    # Tidy cluster hierarchy by merging groups with equal matches
    merge_equal_levels(cluster_list[0])

    # Plot a heat map of the PS matrix
    x = np.arange(0, structure_size + 1, 1)
    y = np.arange(0, structure_size + 1, 1)

    plot = plt.pcolor(x, y, matrix, cmap=plt.get_cmap('rainbow', (n_ps_mols - 1)), vmin=1, vmax=n_ps_mols)
    plt.xticks(np.arange(0, structure_size + 1, 5) - 0.5, np.arange(0, structure_size + 1, 5))
    plt.yticks(np.arange(0, structure_size + 1, 5) - 0.5, np.arange(0, structure_size + 1, 5))
    cb = plt.colorbar(plot, ticks=range(1, 16))

    plt.xlim(0, structure_size)
    plt.ylim(0, structure_size)
    ax = plt.gca()
    ax.set_xlabel("Structure Index", fontsize='x-large')
    ax.set_ylabel("Structure Index", fontsize='x-large')
    cb.set_label('Packing Similarity  /' + str(n_ps_mols) + ' Molecules', fontsize='x-large')
    plt.savefig(input_name + "_heat_map.png", dpi=300)
    print("Packing similarity heat map saved to " + input_name + "_heat_map.png")
    plt.close()

    # Plot a dendrogram
    plot_dendrogram(cluster_list, n_ps_mols, input_name)

    print("--------------------------------------------------------")

    sys.exit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('input_file', help='Set of structures to perform analysis on [.mol2/cif/res/ind]')
    parser.add_argument('-m', '--matrix', type=str, help='NumPy matrix containing existing packing similarity results.',
                        metavar='similarity_matrix.txt')
    parser.add_argument('-ns', '--n_structures', type=int, help='Number of structures to take from input set.',
                        metavar='25')
    parser.add_argument('-nm', '--n_molecules', type=int, default=15,
                        help='Size of molecular packing shell to use for analysis '
                             '(must be consistent with input matrix, if used).', metavar="20")
    parser.add_argument('-o', action="store_true",
                        help='Flag for whether to save packing similarity results (text file and mol2 overlays).')
    parser.add_argument('--allow_molecular_differences', action="store_true",
                        help='Flag for whether to allow for molecular differences between structures (e.g. for salts).')
    parser.add_argument('--clustering_type', choices=['complete', 'single', 'average'], default='single',
                        help='Type of clustering to employ')
    parser.add_argument('-s', '--strip', action="store_true",
                        help='Strip all terminal atoms and alkyl chains, up to any hetero atom (O, N, S) or cyclic '
                             'atom. This cuts the molecule down to core structural features and may reveal'
                             ' more general structural similarities, including those between molecules with different '
                             'conformations.')
    parser.add_argument('-ct', '--conf_tol', type=float, default=0.5,
                        help='RMSD threshold for considering two conformations to be the same (when merging at level 1)'
                             '.', metavar="0.5")
    parser.add_argument('-at', '--angle_tol', type=float, default=25, metavar="25",
                        help="Tolerance for angles (in degrees) used by packing similarity.")
    parser.add_argument('-dt', '--dist_tol', type=float, default=0.25, metavar="0.25",
                        help="Fractional tolerance for distances (0.0 - 1.0) used by packing similarity.")
    args = parser.parse_args()
    if not os.path.isfile(args.input_file):
        parser.error('%s not found.' % args.input_file)
    if args.matrix:
        if not os.path.isfile(args.matrix):
            parser.error('%s not found.' % args.matrix)

    main(args.input_file, args.matrix, args.n_molecules, args.o, args.conf_tol, args.angle_tol,
         args.dist_tol, args.strip, args.n_structures, args.allow_molecular_differences,
         args.clustering_type)
