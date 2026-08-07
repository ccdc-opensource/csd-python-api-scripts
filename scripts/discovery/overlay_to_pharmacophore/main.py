"""
This script can be used for any purpose without limitation subject to the
conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
This permission notice and the following statement of attribution must be
included in all copies or substantial portions of this script.

"07/08/2026": created by the Cambridge Crystallographic Data Centre
"""

import argparse
from pathlib import Path

from cluster import cluster_features
from datastructures import OverlayData
from overlay import OverlayToPharmFeatures
from write_query import FeaturesToCrossMinerQuery


def str_to_bool(value: str) -> bool:
    return value.lower() in {'t', 'true', '1', 'yes', 'y'}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create Pharmacophore Features from a Ligand Overlay"
    )
    parser.add_argument('-i', '--input_folder', type=str, required=True,
                        help='Input file(s) path.')
    parser.add_argument('-o', '--output_folder', type=str, default=None,
                        help="Output folder path. Defaults to a 'queries' folder in the current directory.")
    parser.add_argument('-f', '--feature_definitions', type=str, required=True,
                        help='Path to the directory containing the CrossMiner feature definition (.cpf) files.')
    parser.add_argument('-c', '--cluster', type=str_to_bool, default=False,
                        help='Cluster features if they are close together or common across multiple inputs.')
    parser.add_argument('-p', '--projected', type=str_to_bool, default=False,
                        help='Use projected acceptor features or point features.')
    parser.add_argument('-id', '--overlay_id', type=int, default=0,
                        help='Overlay ID to process. If 0 or not specified, all overlays will be processed.')

    return parser.parse_args()


def main():
    args = parse_args()

    input_folder = Path(args.input_folder)
    if not input_folder.exists():
        raise FileNotFoundError(f"Input folder {input_folder} does not exist.")

    feature_definitions = Path(args.feature_definitions)
    if not feature_definitions.is_dir():
        raise FileNotFoundError(f"Feature definitions folder {feature_definitions} does not exist.")

    output_folder = Path(args.output_folder) if args.output_folder else Path('queries')
    output_folder.mkdir(parents=True, exist_ok=True)

    if (args.overlay_id == 0) or (args.overlay_id is None):
        overlay_files = sorted(input_folder.glob('solution_*.mol2'))
        pharm_files = sorted(input_folder.glob('pharmacophores/solution_pharm_*.mol2'))
    else:
        overlay_files = [input_folder / f'solution_{args.overlay_id:02}.mol2']
        pharm_files = [input_folder / f'pharmacophores/solution_pharm_{args.overlay_id:02}.mol2']
    feature_sets = []
    for pharm_file, overlay_file in zip(pharm_files, overlay_files):
        overlay_data = OverlayData(
            input_folder=input_folder,
            output_folder=output_folder,
            pharm_file=pharm_file,
            overlay_file=overlay_file
        )
        feature_sets.append(OverlayToPharmFeatures(overlay_data, projected=args.projected).features)

    for i, feature_set in enumerate(feature_sets, 1):
        if args.cluster:
            feature_set = cluster_features(feature_set)

        query = FeaturesToCrossMinerQuery(
            pharm_feature_points=feature_set,
            feature_definitions=feature_definitions,
            output_file=output_folder / f'features_{i}.cm',
        )
        query.write_feature_file()


if __name__ == '__main__':
    main()
