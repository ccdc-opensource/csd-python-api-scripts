"""
This script can be used for any purpose without limitation subject to the
conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
This permission notice and the following statement of attribution must be
included in all copies or substantial portions of this script.

"17/11/2026": created by Chris Ringrose for the Cambridge Crystallographic Data Centre
"""

import os

try:
    from matplotlib import pyplot as plt
except ImportError as _error:
    plt = None
    print(f'Plotting disabled (matplotlib unavailable: {_error})')

from ccdc import docking, io

HERE = os.path.dirname(os.path.abspath(__file__))

PROTEIN_FILE = os.path.join(HERE, '4zzn_protein.mol2')
NATIVE_LIGAND_FILE = os.path.join(HERE, '4zzn_ligand.mol2')
LIBRARY_FILE = os.path.join(HERE, 'screening_library.sdf')
OUTPUT_DIRECTORY = os.path.join(HERE, 'output')


def dock_library():
    """Dock every molecule in the library once and return the docked poses."""
    settings = docking.Docker.Settings()
    settings.add_protein_file(PROTEIN_FILE)

    protein = settings.proteins[0]
    native = io.MoleculeReader(NATIVE_LIGAND_FILE)[0]
    native.identifier = 'native_ligand'

    settings.binding_site = settings.BindingSiteFromLigand(protein, native, 6.0)

    # One pose per molecule: in screening you need a score for each compound,
    # not a thorough exploration of any one of them.
    settings.add_ligand_file(LIBRARY_FILE, ndocks=1)

    settings.fitness_function = 'plp'
    # Screening trades search effort for throughput.
    settings.autoscale = 30.
    settings.output_directory = OUTPUT_DIRECTORY
    settings.output_format = 'mol2'
    settings.output_file = 'docked_library.mol2'

    docker = docking.Docker(settings=settings)
    docker.dock(file_name=os.path.join(OUTPUT_DIRECTORY, 'virtual_screening.conf'))
    return docker.results


def roc_curve(scored):
    """Return ROC points and the AUC from (score, is_active) pairs."""
    ranked = sorted(scored, key=lambda pair: pair[0], reverse=True)
    n_actives = sum(1 for _, active in ranked if active)
    n_inactives = len(ranked) - n_actives

    # Walk down the ranked list, stepping up for an active and right otherwise.
    # The area accumulated under those steps is the AUC.
    true_positives = false_positives = 0
    points = [(0.0, 0.0)]
    auc = 0.0
    for _, active in ranked:
        if active:
            true_positives += 1
        else:
            false_positives += 1
            auc += true_positives / n_actives   # area of the strip just swept
        points.append((false_positives / n_inactives, true_positives / n_actives))

    return points, auc / n_inactives


def enrichment_factor(scored, fraction=0.01):
    """How many times more actives appear in the top fraction than by chance."""
    ranked = sorted(scored, key=lambda pair: pair[0], reverse=True)
    cut = max(1, int(round(len(ranked) * fraction)))
    hit_rate = sum(1 for _, active in ranked[:cut] if active) / cut
    overall_rate = sum(1 for _, active in ranked if active) / len(ranked)
    return hit_rate / overall_rate


def main():
    results = dock_library()

    # The identifier prefix records the known activity (see the prep script).
    scored = [
        (pose.fitness(), pose.identifier.startswith('active'))
        for pose in results.ligands
    ]

    points, auc = roc_curve(scored)
    n_actives = sum(1 for _, active in scored if active)
    print(f'Docked {len(scored)} molecules, {n_actives} of them actives')
    print(f'AUC                  {auc:.3f}')
    print(f'Enrichment factor 5%  {enrichment_factor(scored, 0.05):.1f}')

    plot_roc(points, auc)


def plot_roc(points, auc):
    """Save the ROC curve, if matplotlib is available."""
    if plt is None:
        print('Skipping the ROC plot: matplotlib did not import.')
        return

    plt.figure(figsize=(5, 5))
    plt.plot([x for x, _ in points], [y for _, y in points],
             label=f'GOLD PLP (AUC = {auc:.3f})')
    plt.plot([0, 1], [0, 1], '--', color='grey', label='random')
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('MAPK1 virtual screen (LIT-PCBA sample)')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(os.path.join(HERE, 'roc_curve.png'), dpi=150)
    print('Wrote roc_curve.png')


if __name__ == '__main__':
    main()
