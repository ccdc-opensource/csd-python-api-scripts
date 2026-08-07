from pathlib import Path

from ccdc.pharmacophore import Pharmacophore


def search(query_file: Path, database_file: Path):
    """
    This is here as an example for users to perform a simple CrossMiner
    pharmacophore search using the CCDC Python API.
    It is not used in the main workflow.
    """
    settings = Pharmacophore.Search.Settings()
    settings.max_hit_structures = 20
    settings.max_hits_per_structure = 1
    settings.max_hit_rmsd = 1.0
    searcher = Pharmacophore.Search(settings)
    feature_db = Pharmacophore.FeatureDatabase.from_file(database_file)
    query = Pharmacophore.Query.from_file(str(query_file))
    hits = searcher.search(
        model=query,
        database=feature_db,
        verbose=True,
    )
    return hits
