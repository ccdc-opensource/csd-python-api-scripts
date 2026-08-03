from functools import cached_property

import numpy as np

from datastructures import OverlayData, PharmFeaturePoint


class OverlayToPharmFeatures:
    def __init__(self, overlay_data: OverlayData, projected: bool = True):
        """
        Args:
            overlay_data: Pharmacophore and overlay data for a single ligand overlay.
            projected: Whether to use projected acceptor features
                (since acceptors are flexible and vary across overlays,
                it may be better to simply use point acceptors).
        """
        self.overlay_data = overlay_data
        self.projected = projected

    @cached_property
    def features(self) -> list[PharmFeaturePoint]:
        return self._feature_centres_from_pharm()

    def _read_atom_lines_from_pharm(self):
        """
        Generator that yields lines from the ATOM block of a pharmacophore file.
        This block is where the pharmacophore feature points are defined.
        """
        reading_atoms = False
        with open(self.overlay_data.pharm_file) as pharm_file:
            for line in pharm_file:
                line = line.strip()

                if line == "@<TRIPOS>ATOM":
                    reading_atoms = True
                    continue
                elif line.startswith("@<TRIPOS>") and reading_atoms:
                    break  # End of ATOM block

                if reading_atoms and line:
                    yield line.split()

    def _feature_centres_from_pharm(self) -> list[PharmFeaturePoint]:
        """
        Parse the pharmacophore mol2 file and extract feature points centres.

        returns:  List of PharmFeaturePoints, each labelled with the feature name
            (e.g. 'donor_projected', 'acceptor', 'hydrophobe'), and with associated virtual points if present.
        """
        feature_centres = []
        last_feature_point = None

        for parts in self._read_atom_lines_from_pharm():
            _feature_no, feature_name, x, y, z = parts[:5]
            # Extract first 2/3 chars which describe the feature type
            prefix = ''.join(c for c in feature_name if c.isalpha())

            if prefix in ('DON', 'ACC', 'HY'):
                label = self._prefix_to_label(prefix)
                last_feature_point = PharmFeaturePoint(
                    (float(x), float(y), float(z)),
                    label=label,
                )
                feature_centres.append(last_feature_point)
            elif prefix == 'VP':
                if last_feature_point is None:
                    raise ValueError("Virtual point found before its parent feature point.")
                # Acceptor feature with a virtual point is an acceptor projected feature.
                if last_feature_point.label == "acceptor":
                    # If preferring to use point acceptors, ignore the virtual point.
                    if not self.projected:
                        continue
                    last_feature_point.label = "acceptor_projected"
                vp_coords = np.array([x, y, z], dtype=float)
                # If the last feature has no virtual points, assign the virtual point to it.
                if last_feature_point.virtual_point is None:
                    last_feature_point.virtual_point = vp_coords
                # If the last feature does have virtual points, create a new feature.
                else:
                    new_feature = last_feature_point.copy()
                    new_feature.virtual_point = vp_coords
                    feature_centres.append(new_feature)
            elif prefix == 'NL':
                if last_feature_point is None:
                    raise ValueError("Virtual point found before its parent feature point.")
                last_feature_point.virtual_point = np.array([x, y, z], dtype=float)
                last_feature_point.label = 'ring_planar_projected'
            else:
                raise ValueError(f"Unexpected feature type: {feature_name}")
        return feature_centres

    @staticmethod
    def _prefix_to_label(prefix: str) -> str:
        """
        Convert a feature prefix to its corresponding label.
        Acceptor may later be converted to acceptor_projected if it has a virtual point.
        Args:
            prefix: Prefix to convert to label
        returns: Label for the feature type, e.g. 'acceptor', 'donor_projected', 'hydrophobe'
        """
        return {
            'ACC': 'acceptor',
            'DON': 'donor_projected',
            'HY': 'hydrophobe',
        }[prefix]
