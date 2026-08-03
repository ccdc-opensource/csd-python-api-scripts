from pathlib import Path

from datastructures import PharmFeaturePoint


class FeaturesToCrossMinerQuery:
    """Main class for converting a set of features to a CrossMiner feature file."""

    TEMPLATE_FEATURES_FROM_FILE = """
    FEATURE_SUBSTRUCTURE_START

        {}
    FEATURE_SUBSTRUCTURE_END
    """
    def __init__(
            self,
            pharm_feature_points: list[PharmFeaturePoint],
            feature_definitions: Path,
            output_file: Path = Path('test.cm'),
    ):
        """
        Args:
            pharm_feature_points: List of PharmFeaturePoints to be used in the CrossMiner query.
            feature_definitions: Directory containing the CrossMiner feature definition (``.cpf``) files.
                These are not shipped with this repo; supply the location of your CrossMiner installation's
                feature definitions (the directory containing the ``any``, ``protein`` and
                ``small_molecule`` subdirectories).
            output_file: Output file path where the features will be saved.
        """
        self.pharm_feature_points: list[PharmFeaturePoint] = pharm_feature_points
        self.feature_definitions = Path(feature_definitions)
        if not self.feature_definitions.is_dir():
            raise FileNotFoundError(
                f"Feature definitions directory {self.feature_definitions} does not exist."
            )
        self.output_file = output_file

    def _feature_definition_file(self, feature: PharmFeaturePoint) -> Path:
        """
        Locate the ``.cpf`` definition file for a feature within the supplied feature definitions directory.

        Args:
            feature: PharmFeaturePoint whose label determines the definition file name.
        returns: Path to the matching ``.cpf`` file.
        """
        if feature.label is None:
            raise LookupError("Feature label must be set to locate its feature definition file.")

        filename = f'features_{feature.label}.cpf'
        matches = sorted(self.feature_definitions.glob(f'*/{filename}'))
        if not matches:
            # Also allow a flat directory of definition files.
            flat = self.feature_definitions / filename
            if flat.is_file():
                return flat
            raise FileNotFoundError(
                f"No feature definition file '{filename}' found in {self.feature_definitions}."
            )
        return matches[0]


    @staticmethod
    def _create_feature_vector_text(feature: PharmFeaturePoint) -> str:
        """
        Create a formatted string for a feature vector based on the PharmFeaturePoint object.
        Args:
            feature: PharmFeaturePoint containing coordinates, label, and virtual point if applicable.
        """

        result = (
            f"PHARMACOPHORE_FEATURE {feature.label}\n"
            f"PHARMACOPHORE_SPHERE {' '.join(str(i) for i in feature.coordinates)} {feature.weight_parent}\n"
        )
        if feature.virtual_point is not None:
            result += f"PHARMACOPHORE_SPHERE {' '.join(str(i) for i in feature.virtual_point)} {feature.weight_vp}\n"

        result += (
            f"PHARMACOPHORE_FEATURE_SMALL_MOLECULE\n"
            f"PHARMACOPHORE_FEATURE_DESCRIPTION {feature.label}\n\n"
        )

        return result

    def write_feature_file(self):
        """
        Create a feature file with proper whitespace formatting.
        """
        output_data = 'FEATURE_LIBRARY_START\n'

        # Collect unique feature types and their corresponding files
        feature_files = dict()
        for feature in self.pharm_feature_points:
            if feature.label in feature_files:
                continue
            feature_file = self._feature_definition_file(feature)
            feature_files[feature.label] = feature_file

            # Append the feature data to the output stream
            with open(feature_file) as ff:
                output_data += self.TEMPLATE_FEATURES_FROM_FILE.format(
                    "        ".join(ff.readlines())
                )

        output_data += '\nFEATURE_LIBRARY_END\n\n'

        for feature in self.pharm_feature_points:
            output_data += self._create_feature_vector_text(feature)

        with open(self.output_file, 'w') as output_file:
            output_file.write(output_data)
