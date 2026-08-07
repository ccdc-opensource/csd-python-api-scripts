"""
This script can be used for any purpose without limitation subject to the
conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
This permission notice and the following statement of attribution must be
included in all copies or substantial portions of this script.

"07/08/2026": created by the Cambridge Crystallographic Data Centre
"""

from pathlib import Path

import pytest

from cluster import cluster_features
from datastructures import OverlayData, PharmFeaturePoint
from overlay import OverlayToPharmFeatures
from write_query import FeaturesToCrossMinerQuery

TESTS_DIR = Path(__file__).parent
REPO_ROOT = TESTS_DIR.parent
JOB_FOLDER = TESTS_DIR / "job1"


@pytest.fixture(autouse=True)
def run_from_repo_root(monkeypatch):
    """Run every test from the repository root so ``feature_definitions/`` resolves."""
    monkeypatch.chdir(REPO_ROOT)


@pytest.fixture
def overlay_data(tmp_path) -> OverlayData:
    """Overlay/pharmacophore file locations for the job1 example."""
    return OverlayData(
        input_folder=JOB_FOLDER,
        output_folder=tmp_path,
        pharm_file=JOB_FOLDER / "pharmacophores/solution_pharm_01.mol2",
        overlay_file=JOB_FOLDER / "solution_01.mol2",
    )


@pytest.fixture
def features(overlay_data) -> list[PharmFeaturePoint]:
    """Raw (unclustered) features extracted from the example pharmacophore."""
    return OverlayToPharmFeatures(overlay_data).features


@pytest.fixture
def feature_definitions(features, tmp_path) -> Path:
    """
    Stub CrossMiner feature definition directory.

    The real ``.cpf`` files ship with CrossMiner rather than this repo, so minimal
    placeholders are generated for each label present in the example pharmacophore.
    """
    definitions = tmp_path / "feature_definitions" / "any"
    definitions.mkdir(parents=True)
    for label in {f.label for f in features}:
        (definitions / f"features_{label}.cpf").write_text(f"SUBSTRUCTURE {label}\n")
    return definitions.parent


def test_features_extracted_from_pharmacophore(features):
    """Parsing the pharmacophore mol2 yields labelled feature points."""
    assert len(features) > 0, "No features extracted from pharmacophore."
    assert all(feature.label is not None for feature in features)


def test_projected_features_have_virtual_points(features):
    """Every feature labelled as projected carries a virtual point."""
    projected = [f for f in features if f.label.endswith("projected")]
    assert projected, "No projected features extracted from pharmacophore."
    assert all(f.virtual_point is not None for f in projected)


def test_unprojected_acceptors_have_no_virtual_points(overlay_data):
    """With `projected=False` acceptors are kept as single points."""
    features = OverlayToPharmFeatures(overlay_data, projected=False).features
    acceptors = [f for f in features if f.label.startswith("acceptor")]
    assert acceptors, "No acceptor features extracted from pharmacophore."
    assert all(f.label == "acceptor" and f.virtual_point is None for f in acceptors)


def test_clustering_reduces_features(features):
    """Clustering merges nearby features without dropping any feature type."""
    clustered = cluster_features(features)
    assert 0 < len(clustered) <= len(features)
    assert {f.label for f in clustered} == {f.label for f in features}


def test_query_file_written(features, feature_definitions, tmp_path):
    """A CrossMiner query file is written containing the clustered features."""
    output_file = tmp_path / "query.cm"
    clustered = cluster_features(features)

    FeaturesToCrossMinerQuery(
        clustered, feature_definitions=feature_definitions, output_file=output_file
    ).write_feature_file()

    assert output_file.exists()
    contents = output_file.read_text()
    assert "FEATURE_LIBRARY_START" in contents
    assert "FEATURE_LIBRARY_END" in contents
    assert contents.count("PHARMACOPHORE_FEATURE ") == len(clustered)


def test_missing_feature_definitions_raises(features, tmp_path):
    """An error is raised when the supplied feature definitions directory is absent."""
    with pytest.raises(FileNotFoundError):
        FeaturesToCrossMinerQuery(
            features,
            feature_definitions=tmp_path / "does_not_exist",
            output_file=tmp_path / "query.cm",
        )
