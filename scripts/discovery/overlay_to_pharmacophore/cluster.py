from collections import defaultdict

import numpy as np

from datastructures import PharmFeaturePoint


def cluster_features(features: list[PharmFeaturePoint]) -> list[PharmFeaturePoint]:
    """Cluster all features where relevant."""
    groups = defaultdict(list)
    for feature in features:
        groups[feature.label].append(feature)

    clustered_features = []
    for label, group in groups.items():
        if label == 'hydrophobe':
            clustered_features.extend(_cluster_similar_features(group))
        elif label in ('acceptor_projected', 'donor_projected'):
            clustered_features.extend(_cluster_similar_features(group, check_vp=True))
        else:
            clustered_features.extend(group)
    return clustered_features


def _cluster_similar_features(
        features: list[PharmFeaturePoint], cluster_radius: float = 2.0, check_vp: bool = False
) -> list[PharmFeaturePoint]:
    """
    Cluster specific features that are close to each other.
    For projected features, also check virtual point distances.

    Args:
        features: List of features to cluster
        cluster_radius: Distance between features to be included in clustering (Ang)
        check_vp: Whether to also require virtual points to be within the radius
    """
    if len(features) < 2:
        return features

    clusters = _connected_components(
        features,
        lambda a, b: _is_close(a, b, cluster_radius, check_vp),
    )
    return [_merge_cluster(cluster, check_vp) for cluster in clusters]


def _is_close(
        a: PharmFeaturePoint, b: PharmFeaturePoint, cluster_radius: float, check_vp: bool
) -> bool:
    """Whether two features are within ``cluster_radius`` (and, if ``check_vp``, their virtual points too)."""
    if np.linalg.norm(a - b) >= cluster_radius:
        return False
    if check_vp and np.linalg.norm(a.virtual_point - b.virtual_point) >= cluster_radius:
        return False
    return True


def _connected_components(
        features: list[PharmFeaturePoint], is_linked
) -> list[list[PharmFeaturePoint]]:
    """
    Single-linkage grouping: features are placed in the same cluster if a chain of
    ``is_linked`` neighbours connects them.
    """
    unclustered = features.copy()
    clusters = []

    while unclustered:
        cluster = [unclustered.pop(0)]
        # Grow the cluster breadth-first: any unclustered feature linked to a member joins it.
        i = 0
        while i < len(cluster):
            member = cluster[i]
            remaining = []
            for feature in unclustered:
                if is_linked(feature, member):
                    cluster.append(feature)
                else:
                    remaining.append(feature)
            unclustered = remaining
            i += 1
        clusters.append(cluster)

    return clusters


def _merge_cluster(cluster: list[PharmFeaturePoint], check_vp: bool) -> PharmFeaturePoint:
    """Collapse a cluster into a single feature at its centroid (singletons are returned unchanged)."""
    if len(cluster) == 1:
        return cluster[0]

    centroid = np.mean([c.coordinates for c in cluster], axis=0).round(4)
    vp_centroid = (
        np.mean([c.virtual_point for c in cluster], axis=0).round(4) if check_vp else None
    )
    return PharmFeaturePoint(centroid, label=cluster[0].label, virtual_point=vp_centroid)
