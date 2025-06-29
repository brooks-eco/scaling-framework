"""
Murray-Darling Basin Scaling-up Framework - Core Spatial Analysis Functions

=============================================================================
OVERVIEW
=============================================================================
This framework implements hierarchical spatial aggregation methods for wetland
and landscape ecology research. It enables scaling-up of fine-resolution ecological
data (e.g., individual wetlands, vegetation patches) to broader management units
(e.g., river reaches, catchments, bioregions).

The framework addresses the fundamental challenge in landscape ecology of
transferring information across spatial scales while preserving ecological
relevance and statistical validity.

=============================================================================
CORE CAPABILITIES
=============================================================================
1. Multi-scale Spatial Data Management:
   - Automated geometry validation and CRS standardization
   - Support for points (monitoring sites), lines (rivers), polygons (wetlands)
   - Hierarchical scale relationships (wetlands → reaches → catchments)

2. Ecologically-Informed Aggregation Methods:
   - Area-weighted means for habitat quality metrics
   - Length-weighted aggregation for riparian corridor analysis
   - Frequency-weighted analysis for species occurrence data
   - Custom weighting schemes for ecological importance

3. Functional Group Classification:
   - Substring-based reclassification (e.g., 'Permanent Wetland' → 'Wetland')
   - Wetland type grouping (ephemeral, seasonal, permanent)
   - Vegetation community aggregation

4. Statistical Validation and Quality Control:
   - Automated detection of invalid geometries
   - Duplicate record identification
   - Missing data handling with ecological context

=============================================================================
ECOLOGICAL APPLICATIONS
=============================================================================
- Wetland condition assessment across river basins
- Vegetation health monitoring from plot to landscape scale
- Species habitat suitability modeling at multiple scales
- Water quality indicator aggregation
- Biodiversity pattern analysis across management boundaries

=============================================================================
TECHNICAL IMPLEMENTATION
=============================================================================
Built on GeoPandas and Pandas for robust spatial data handling, with
optimized algorithms for large-scale ecological datasets common in
landscape ecology research.

Author: Shane Brooks
Version: 2.0
License: Open Source  - Creative Commons v4 with attribution CCby4
"""

import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import re
import math
import copy
from matplotlib.lines import Line2D
from typing import List, Optional, Union
from pathlib import Path
from geopandas import GeoDataFrame
from shapely.validation import make_valid, explain_validity
from pyproj import CRS


class SpatialScale:
    """
    Represents a single spatial scale in hierarchical ecological analysis.

    This class encapsulates spatial ecological data at one organizational level
    (e.g., individual wetlands, river reaches, catchments) and provides methods
    for scale-appropriate aggregation and analysis.

    In landscape ecology, spatial scales represent different levels of ecological
    organization. For example:
    - Fine scale: Individual wetland polygons with vegetation condition scores
    - Intermediate scale: River reaches containing multiple wetlands
    - Broad scale: Catchments encompassing multiple river reaches

    Key Ecological Concepts:
    - Scale dependency: Ecological patterns vary with observation scale
    - Hierarchical organization: Smaller units nest within larger units
    - Emergent properties: Higher-level patterns emerge from lower-level processes

    Technical Features:
    - Automatic geometry type detection and validation
    - Ecologically-appropriate weighting (area for habitats, length for corridors)
    - Functional group management for ecological classification
    - Spatial join optimization for large datasets
    - Results caching for computational efficiency

    Examples:
        >>> # Create wetland base scale
        >>> wetlands = SpatialScale(
        ...     name="wetlands",
        ...     source="wetlands.shp",
        ...     unique_id_fields="WETLAND_ID",
        ...     type_field="WETLAND_TYPE",
        ...     is_base_scale=True
        ... )

        >>> # Create catchment aggregation scale
        >>> catchments = SpatialScale(
        ...     name="catchments",
        ...     source="catchments.shp",
        ...     unique_id_fields="CATCH_ID"
        ... )

    Attributes:
        name (str): Human-readable identifier for this scale
        data (GeoDataFrame): Spatial data with computed weights
        geometry_type (str): 'polygon', 'line', or 'point'
        weighting_field (str): Field used for ecological weighting
        type_field (str): Field containing functional groups
        results (dict): Cached aggregation results
    """

    def __init__(
        self,
        name: str,
        source: Union[str, Path],
        unique_id_fields: Union[str, List[str]],
        weighting_field: Optional[str] = None,
        metric_fields: Optional[List[str]] = None,
        measure_multiplier: Optional[float] = None,
        type_field: Optional[str] = None,
        default_crs: str = None,
        is_base_scale: bool = False,
    ):
        """
        Initialize a spatial scale for hierarchical ecological analysis.

        This constructor loads and validates spatial ecological data, automatically
        detecting geometry types and computing ecologically-appropriate weights.

        The initialization process follows landscape ecology best practices:
        1. Spatial data validation and CRS standardization
        2. Geometry type detection (critical for appropriate weighting)
        3. Automatic weight calculation based on ecological relevance:
           - Polygons (habitats): Area in hectares
           - Lines (corridors): Length in meters
           - Points (sites): Count-based weighting
        4. Data quality validation and error reporting

        Args:
            name: Identifier for this spatial scale (e.g., 'wetlands', 'catchments')
            source: Path to shapefile containing spatial features
            unique_id_fields: Column(s) that uniquely identify each feature.
                For compound keys use list: ['REGION_ID', 'SITE_ID']
            weighting_field: Column containing ecological weights. If None,
                automatically computed based on geometry type
            metric_fields: Columns containing ecological metrics to aggregate.
                Examples: ['NDVI_mean', 'species_richness', 'water_quality']
            measure_multiplier: Scaling factor for weights (e.g., 0.0001 to
                convert m² to hectares)
            type_field: Column containing functional groups for ecological
                classification (e.g., 'WETLAND_TYPE', 'VEG_COMMUNITY')
            default_crs: Coordinate reference system. Defaults to EPSG:3577
                (Australian Albers) for MDB region
            is_base_scale: True for finest resolution data, False for
                aggregation targets

        Raises:
            FileNotFoundError: If shapefile doesn't exist
            ValueError: If required columns missing or data invalid

        Examples:
            >>> # Wetland base scale with vegetation metrics
            >>> wetlands = SpatialScale(
            ...     name="wetlands",
            ...     source="MDB_wetlands.shp",
            ...     unique_id_fields="WETLAND_ID",
            ...     type_field="WETLAND_TYPE",
            ...     metric_fields=['NDVI_2023', 'condition_score'],
            ...     is_base_scale=True
            ... )

            >>> # River reach aggregation scale
            >>> reaches = SpatialScale(
            ...     name="river_reaches",
            ...     source="MDB_reaches.shp",
            ...     unique_id_fields="REACH_ID"
            ... )
        """
        self.name = name
        self.source = Path(source)
        self.is_base_scale = is_base_scale

        # Check that the shapefile actually exists
        if not self.source.exists():
            raise FileNotFoundError(
                f"Cannot find shapefile: {source}\n"
                f"Please check that the file exists and the path is correct."
            )

        # Store unique_id as LIST (convert from string if needed)
        self.unique_id_fields = (
            [unique_id_fields]
            if isinstance(unique_id_fields, str)
            else unique_id_fields
        )

        # =====================================================================
        # STEP 1: DETERMINE REQUIRED COLUMNS FOR EFFICIENT DATA LOADING
        # =====================================================================
        # In landscape ecology, datasets can be very large (thousands of
        # wetlands, millions of vegetation plots). We optimize performance by
        # loading only the columns needed for analysis, reducing memory usage
        # and processing time for large ecological datasets.
        required_cols = set()
        required_cols.update(self.unique_id_fields)
        if weighting_field:
            required_cols.add(weighting_field)
        if type_field:
            required_cols.add(type_field)
        if metric_fields:
            if isinstance(metric_fields, str):
                metric_fields = [metric_fields]
            required_cols.update(metric_fields)

        # Ensure geometry is included (required for GeoPandas)
        required_cols = list(required_cols)
        required_cols.append("geometry")  # geometry must be present

        print(f"Loading {name} data from {self.source.name}...")
        print(f"   Requested: {required_cols}")

        # =====================================================================
        # STEP 2: COORDINATE REFERENCE SYSTEM STANDARDIZATION
        # =====================================================================
        # Consistent CRS is critical for accurate spatial analysis in ecology.
        # EPSG:3577 (Australian Albers) is optimal for area calculations across
        # the Murray-Darling Basin, minimizing distortion for ecological metrics.
        try:
            self.default_crs = (
                CRS.from_user_input(default_crs) if default_crs else CRS.from_epsg(3577)
            )
        except:
            print(
                f"   Invalid CRS '{default_crs}', using EPSG:3577 (Australian Albers)"
            )
            self.default_crs = CRS.from_epsg(3577)
        print(f"   Using CRS: {self.default_crs.name}.")

        # =====================================================================
        # STEP 3: LOAD AND VALIDATE SPATIAL DATA
        # =====================================================================
        # Load spatial data with geometry validation to ensure robust analysis.
        # Invalid geometries can cause errors in spatial joins and aggregation.
        self.data = self._validate_geometries(
            gpd.read_file(self.source).to_crs(self.default_crs), fix_invalid=True
        )

        # # Convert columns to strings (if not already)
        # self.data.columns = self.data.columns.map(str) # now in validate_metric_columns

        # # Convert required_cols to strings to match DataFrame columns
        # required_cols = [str(c) for c in required_cols] # now in validate_metric_columns

        # Select columns based on whether metric_fields is specified
        if metric_fields:
            # Only keep required columns when metric_fields is specified
            cols_to_keep = list(set(required_cols) | {"geometry"})
            self.data = self.data.loc[:, cols_to_keep]
        print(f"   No specific metric fields chosen so loading all columns.")
        cols = [col for col in self.data.columns if col != "geometry"]
        print(f"   Loaded: {cols}")
        # If metric_fields is None, retain all columns (no filtering)

        multi_geom_count = self.data.geometry.geom_type.isin(
            ["MultiPolygon", "MultiLineString", "MultiPoint"]
        ).sum()
        if multi_geom_count:
            print(
                f"   Found {multi_geom_count} multi-part geometries in {self.name} data.  Calculations will proceed OK. Just letting you know!  "
            )
        # commented out explode logic but retained for future reference
        # # Handle complex geometries by breaking them into simple parts
        # # E.g., a MultiPolygon becomes multiple separate Polygons
        # # This ensures each feature is treated as a separate unit for analysis
        # if multi_geom_count > 0:
        #     print(
        #         f"   Breaking {multi_geom_count} complex geometries into simple parts..."
        #     )
        #     self.data = self.data.explode(index_parts=False).reset_index(drop=True)

        # # Check that all required columns exist in the shapefile
        # now in validate_metric_columns
        # missing_cols = [c for c in required_cols if c not in self.data.columns]
        # if missing_cols:
        #     raise ValueError(
        #         f"Required columns {missing_cols} not found in {source}\n"
        #         f"Available columns: {list(self.data.columns)}\n"
        #         f"Please check your column names in config.py"
        #     )

        # =====================================================================
        # STEP 4: GEOMETRY TYPE DETECTION FOR ECOLOGICAL WEIGHTING
        # =====================================================================
        # Different ecological features require different weighting approaches:
        # - Wetland polygons: Area-based (habitat extent)
        # - River lines: Length-based (corridor connectivity)
        # - Monitoring points: Count-based (site representation)
        # This detection ensures ecologically-appropriate aggregation methods.
        all_types = set(self.data.geometry.geom_type.unique())
        valid_polygon = all_types <= {"Polygon", "MultiPolygon"}
        valid_line = all_types <= {"LineString", "MultiLineString"}
        valid_point = all_types <= {"Point", "MultiPoint"}

        if valid_polygon:
            geometry_type = "polygon"
        elif valid_line:
            geometry_type = "line"
        elif valid_point:
            geometry_type = "point"
        else:
            raise ValueError(f"Inconsistent or unsupported geometry types: {all_types}")

        print(f"   Detected geometry type: {geometry_type} ({len(self.data)} features)")

        # =====================================================================
        # STEP 5: COMPUTE ECOLOGICALLY-APPROPRIATE WEIGHTS
        # =====================================================================
        # Automatic weight calculation follows ecological principles:
        # - Larger habitats have greater ecological influence (area weighting)
        # - Longer corridors provide more connectivity (length weighting)
        # - Each monitoring site contributes equally (count weighting)
        if not weighting_field:
            if geometry_type == "polygon":
                self.data["Area_Ha"] = self.data.geometry.area / 10000
                weighting_field = "Area_Ha"
            elif geometry_type == "line":
                self.data["Length_m"] = self.data.geometry.length
                weighting_field = "Length_m"
            elif geometry_type == "point":
                self.data["Count"] = 1
                weighting_field = "Count"

        # Apply measure multiplier if specified (e.g., convert units)
        if measure_multiplier and weighting_field:
            print(f"   Applying multiplier {measure_multiplier} to {weighting_field}")
            self.data[weighting_field] *= measure_multiplier

        self.metric_fields = metric_fields or []  # Metrics available for aggregation
        self.weighting_field = (
            weighting_field  # Field used for weighting (area, length, etc.)
        )
        self.type_field = (
            type_field  # Optional field for type-based grouping/reclassification
        )

        # Centralized validation of metrics and UID uniqueness
        self.data = self._validate_metric_data(
            self.data,
            unique_id=self.unique_id_fields,
            metric_field_names=self.metric_fields,
            not_allowed=[self.weighting_field, self.type_field],
            set_index=False,
        )

        self.geometry_type = geometry_type
        self._join_cache = {}  # Cached spatial join results
        self.results = {}  # Stores aggregation results keyed by result name

    def __str__(self) -> str:
        """Create a human-readable description of this spatial scale."""
        scale_type = "BASE SCALE" if self.is_base_scale else "AGGREGATION SCALE"
        return f"{self.name} ({scale_type}) - {len(self.data)} features from {self.source.name}"

    def _validate_geometries(
        self, gdf: gpd.GeoDataFrame, fix_invalid: bool = True
    ) -> gpd.GeoDataFrame:
        """
        Validate and optionally fix invalid geometries.

        Args:
            gdf: GeoDataFrame to validate.
            fix_invalid: Whether to attempt fixing invalid geometries.

        Returns:
            GeoDataFrame with validated geometries.
        """
        invalid_mask = ~gdf.geometry.is_valid
        invalid_count = invalid_mask.sum()

        if invalid_count == 0:
            print("   All geometries are valid.")
            return gdf

        print(f"\u26a0 Found {invalid_count} invalid geometries.")

        if fix_invalid:
            print("   Attempting to fix invalid geometries...")
            gdf.loc[invalid_mask, "geometry"] = gdf.loc[invalid_mask, "geometry"].apply(
                make_valid
            )

            still_invalid = (~gdf.geometry.is_valid).sum()
            if still_invalid == 0:
                print("   All geometries successfully fixed.")
            else:
                print(
                    f"\u26a0 Warning: {still_invalid} geometries remain invalid after repair."
                )

        return gdf

    @staticmethod
    def validate_result_name(name: str, existing_keys: Optional[set] = None) -> str:
        if not isinstance(name, str) or not name.strip():
            raise ValueError("result_name must be a non-empty string.")
        if not re.match(r"^[\w\-]+$", name):
            raise ValueError(
                "result_name can only contain letters, numbers, underscores, or dashes."
            )

        if existing_keys and name in existing_keys:
            orig_name = name
            for i in range(1, 11):
                name = f"{orig_name}_{i}"
                if name not in existing_keys:
                    print(
                        f"\u26a0 Warning Result '{orig_name}' exists. Using '{name}' instead."
                    )
                    return name
            raise ValueError(
                f"Too many results named '{orig_name}'. Choose a different name."
            )

        return name

    def join_data(
        self,
        data_path: Union[str, Path],
        tabular_data_unique_id_fields: Optional[Union[str, List]] = None,
        pivot_row_id: Optional[List[str]] = None,
        pivot_columns: Optional[str] = None,
        pivot_values: Optional[str] = None,
        how: str = "inner",
    ) -> None:
        """
        Load external tabular data and join it to the SpatialScale's data.

        Supports:
        1. Direct join if data is in wide format (one row per unique_id).
        2. Pivoting long data (key-value) into wide before join.

        Args:
            data_path: Path to the CSV or tabular data file.
            unique_id_fields: Column in the CSV matching source unique ID. If None, defaults to self.unique_id.
            pivot_index: If pivoting, list of columns to use as index.
            pivot_columns: If pivoting, column name whose values become new columns.
            pivot_values: If pivoting, column name containing values to fill.
            how: Type of join to perform (default "left").

        Raises:
            ValueError if required columns not found.
        """

        # Load CSV data
        tabular_data = pd.read_csv(data_path)

        # Determine the ID field to use fall back to unique_id_fields

        tabular_data_unique_id_fields = (
            [tabular_data_unique_id_fields]
            if isinstance(tabular_data_unique_id_fields, str)
            else tabular_data_unique_id_fields or self.unique_id_fields
        )

        # if not tabular_data_unique_id_fields:
        #     tabular_data_unique_id_fields = self.unique_id_fields

        pivot_row_id = [pivot_row_id] if isinstance(pivot_row_id, str) else pivot_row_id

        # Check all unique_id fields exist in the DataFrame columns
        if not all(
            field in tabular_data.columns for field in tabular_data_unique_id_fields
        ):
            raise ValueError(
                f"unique_id_fields '{tabular_data_unique_id_fields}' not found in data columns: {list(tabular_data.columns)}"
            )

        # look for duplicate rows with same tabular_data_unique_id_fields (which may be a list) as indictor
        # data is in long format and needs to be pivoted
        duplicated_mask = tabular_data.duplicated(
            subset=tabular_data_unique_id_fields, keep=False
        )
        if not pivot_row_id and duplicated_mask.any():
            print(
                "\u26a0 Warning: Duplicate IDs found in tabular data. Consider pivoting?"
            )

        # Pivot if requested (assuming long format data)
        if pivot_row_id and pivot_columns and pivot_values:
            missing_cols = [
                c
                for c in [*pivot_row_id, pivot_columns, pivot_values]
                if c not in tabular_data.columns
            ]
            if missing_cols:
                raise ValueError(f"Missing columns required for pivot: {missing_cols}")

            tabular_data = tabular_data.pivot_table(
                index=pivot_row_id,
                columns=pivot_columns,
                values=pivot_values,
                aggfunc="first",
            ).reset_index()

        # ensure all column headers are strings
        tabular_data_clean = self._validate_metric_data(
            tabular_data,
            unique_id=self.unique_id_fields,
            not_allowed=[self.weighting_field, self.type_field],
        )

        # Single merge operation handles both same and different column names
        self.data = self.data.merge(
            tabular_data_clean,
            how=how,
            left_on=self.unique_id_fields,
            right_on=tabular_data_unique_id_fields,
            suffixes=("", "_external"),
            validate="one_to_one",
        )

        print(
            f"Joined external data from {data_path} on '{tabular_data_unique_id_fields}'"
        )

    def spatial_join(
        self, target_scale: "SpatialScale", how: str = "intersects"
    ) -> GeoDataFrame:
        """
        Perform and cache a spatial join between this scale and a target scale.

        Args:
            target_scale: The target SpatialScale to join with.
            how: The spatial predicate (e.g., "intersects", "contains").

        Returns:
            GeoDataFrame: Spatially joined data.
        """
        cache_key = f"{target_scale.name}_{how}"

        if cache_key not in self._join_cache:
            print(
                f"Performing spatial join: {self.name} to {target_scale.name} using '{how}'"
            )

            # Resolve column name conflicts
            source_cols = set(self.data.columns) - {"geometry"}
            target_cols = set(target_scale.data.columns) - {"geometry"}
            conflicting = source_cols & target_cols
            rename_map = {col: f"{col}_target" for col in conflicting}

            if rename_map:
                print(
                    f"\u26a0 Warning: Renaming {len(rename_map)} conflicting columns in target scale: {rename_map}"
                )
                target_scale.data = target_scale.data.rename(columns=rename_map)

            # Update target's unique_id fields (noting uid is a list) if renamed
            # Handle compound key (list of fields)
            updated_uids = [
                rename_map.get(uid, uid) for uid in target_scale.unique_id_fields
            ]

            # If any field was renamed, update the unique_id_field
            if updated_uids != target_scale.unique_id_fields:
                target_scale.unique_id_fields = updated_uids
                print(
                    f"\u26a0 Warning: Updated target scale unique_id_field to '{updated_uids}'"
                )

            # Run the spatial join (inner join only, to maintain overlaps)
            joined = gpd.sjoin(self.data, target_scale.data, how="inner", predicate=how)
            self._join_cache[cache_key] = joined

        else:
            print(f"Using cached spatial join: {cache_key}")

        return self._join_cache[cache_key]

    def _apply_reclassification(
        self,
        df: pd.DataFrame,
        reclass_map: dict,
        keep_unmatched_types: bool = True,
        new_class_field: str = "regrouped",
    ) -> tuple[pd.DataFrame, str]:
        """
        Apply functional group reclassification using substring matching.

        This method implements flexible ecological classification by grouping
        similar habitat types, species, or management categories based on
        substring patterns in their names. This is essential in landscape
        ecology where:

        1. **Taxonomic grouping**: Species names may vary but belong to same family
        2. **Habitat classification**: Detailed habitat types need broader groupings
        3. **Management categories**: Administrative names need functional groupings

        The substring matching approach handles common ecological data challenges:
        - Inconsistent naming conventions across datasets
        - Hierarchical classification systems (species → genus → family)
        - Regional variations in habitat terminology

        Ecological Examples:
        - Wetland types: 'Permanent Freshwater Marsh' → 'Wetland'
        - Vegetation: 'River Red Gum Woodland' → 'Woodland'
        - Species: 'Acacia melanoxylon' → 'Acacia'

        Args:
            df: DataFrame containing ecological data to reclassify
            reclass_map: Dictionary mapping new group names to lists of
                substring patterns. Example:
                {
                    'Wetland': ['marsh', 'swamp', 'wetland'],
                    'Woodland': ['woodland', 'forest', 'trees']
                }
            keep_unmatched_types: If True, retains original names for
                unmatched types. If False, removes unmatched records.
                Critical decision for maintaining data completeness vs.
                analytical focus.
            new_class_field: Name for the new classification column

        Returns:
            tuple: (modified DataFrame, field name for subsequent grouping)

        Raises:
            None: Method handles all edge cases gracefully with warnings

        Examples:
            >>> # Wetland functional grouping
            >>> wetland_groups = {
            ...     'Permanent': ['permanent', 'deep', 'lake'],
            ...     'Seasonal': ['seasonal', 'ephemeral', 'temporary'],
            ...     'Constructed': ['dam', 'channel', 'artificial']
            ... }
            >>> df, field = scale._apply_reclassification(
            ...     df, wetland_groups, keep_unmatched_types=True
            ... )

        Notes:
            - Case-insensitive matching for robust classification
            - First match wins (order matters in reclass_map)
            - Provides detailed mapping report for transparency
            - Handles missing or null values gracefully
        """
        # =====================================================================
        # EARLY EXIT: NO RECLASSIFICATION NEEDED
        # =====================================================================
        # If no reclassification map provided or no type field exists,
        # return original data unchanged. This allows the same aggregation
        # code to work with and without functional grouping.
        if not reclass_map or not self.type_field:
            return df, self.type_field

        # =====================================================================
        # INITIALIZE RECLASSIFICATION PROCESS
        # =====================================================================
        print(
            f"Reclassifying '{self.type_field}' into new groups '{new_class_field}' using substring match from reclass map..."
        )

        # Clean up any existing classification column to avoid conflicts
        # This ensures reproducible results when re-running analysis
        if new_class_field in df.columns:
            df.drop(columns=new_class_field, inplace=True)
            print(
                f"\u26a0 Warning: removed existing '{new_class_field}' column and re-applying group rules."
            )

        # =====================================================================
        # DEFINE SUBSTRING MATCHING FUNCTION
        # =====================================================================
        # This nested function implements the core classification logic.
        # It uses case-insensitive substring matching to handle variations
        # in ecological naming conventions (e.g., 'Wetland' vs 'wetland').
        def match_reclass(value):
            """Match ecological type to functional group using substring patterns."""
            # Convert to string to handle numeric codes or mixed data types
            value_str = str(value).lower()

            # Iterate through functional groups in order (first match wins)
            for group, substrings in reclass_map.items():
                for substr in substrings:
                    if substr.lower() in value_str:
                        return group

            # Return None if no pattern matches (handled below)
            return None

        # =====================================================================
        # APPLY CLASSIFICATION TO ALL RECORDS
        # =====================================================================
        # Use pandas apply() for vectorized operation across all ecological types
        df[new_class_field] = df[self.type_field].apply(match_reclass)

        # =====================================================================
        # HANDLE UNMATCHED ECOLOGICAL TYPES
        # =====================================================================
        # Identify types that didn't match any classification pattern.
        # This is common in ecological data due to:
        # - Rare or unusual habitat types
        # - Regional naming variations
        # - Data entry inconsistencies
        # - New types not in original classification scheme

        unmatched = df[df[new_class_field].isna()][self.type_field].unique()

        if len(unmatched) > 0:
            if not keep_unmatched_types:
                # OPTION 1: Remove unmatched types (focused analysis)
                # Use when analysis requires only well-defined functional groups
                # Common in comparative studies or standardized assessments
                df = df[~df[new_class_field].isna()].copy()
                print(
                    f"Dropped {len(unmatched)} types that did not match any re-grouping rules. "
                    f"Keeping {len(df)} re-grouped rows.\n"
                    f"Set option 'keep_unmatched_types=True' to retain unmatched types."
                )
            else:
                # OPTION 2: Retain unmatched types with original names
                # Use when comprehensive coverage is more important than
                # standardized grouping. Preserves all ecological information.
                df[new_class_field] = df[new_class_field].fillna(df[self.type_field])
                print(
                    f"Retained {len(unmatched)} original types that did not match any re-grouping rules. "
                    f"Set option 'keep_unmatched_types=False' to drop these rows."
                )

        # =====================================================================
        # GENERATE CLASSIFICATION REPORT FOR TRANSPARENCY
        # =====================================================================
        # Provide detailed mapping report to ensure classification accuracy.
        # This is critical in ecological research for:
        # - Verifying classification logic
        # - Documenting methodological decisions
        # - Enabling reproducible research
        # - Identifying potential classification errors

        mapping = (
            df[[self.type_field, new_class_field]]
            .drop_duplicates()
            .sort_values(self.type_field)
        )

        print(
            f"\nClassification mapping applied: '{self.type_field}' → '{new_class_field}'"
        )
        print("=" * 60)

        # Format output for readability with aligned columns
        max_width = max(
            mapping[self.type_field].astype(str).map(len).max(),
            mapping[new_class_field].astype(str).map(len).max(),
        )

        # Display each mapping with clear formatting
        for _, row in mapping.iterrows():
            original = f"'{row[self.type_field]}'"
            classified = f"'{row[new_class_field]}'"
            print(f"  {original.ljust(max_width + 3)} → {classified}")

        print(f"\nTotal unique types: {len(mapping)}")
        print(f"Functional groups created: {mapping[new_class_field].nunique()}")

        # =====================================================================
        # RETURN CLASSIFIED DATA AND FIELD NAME
        # =====================================================================
        # Return both the modified DataFrame and the field name to use for
        # subsequent grouping operations. This tuple approach ensures the
        # calling code knows which field contains the functional groups.
        return df, new_class_field

    def _aggregate_joined(
        self,
        df: GeoDataFrame,
        target_scale: "SpatialScale",
        metric_columns: list[str],
        method: str = "weighted_mean",
        reclass_map: dict[str, str] = None,
        new_class_field: str = "regrouped",
        group_by: list[str] = None,
        keep_unmatched_types: bool = True,
        weighting_field: Optional[str] = None,
        result_name: str = None,
    ) -> pd.DataFrame:
        """
        Aggregate ecological metrics from spatially joined data using
        ecologically-appropriate weighting schemes.

        Aggregated_Value = Σ(Metric_i × Weight_i) / Σ(Weight_i)

        Where weights vary by method:
        - Area-weighted: Weight_i = Area_i (hectares)
        - Length-weighted: Weight_i = Length_i (meters)
        - Count-based: Weight_i = 1 (equal weighting)
        - Frequency-weighted: Weight_i = Frequency_i / Total_frequency

        Args:
            df: Spatially joined GeoDataFrame containing source features
                overlaid with target boundaries. Must include both source
                metrics and target identifiers.
            target_scale: SpatialScale object representing the aggregation
                target (e.g., catchments, management units)
            metric_columns: List of ecological metrics to aggregate.
                Examples: ['NDVI_mean', 'species_richness', 'pH', 'condition_score']
            method: Aggregation method selection:
                - 'area_weighted': Area-proportional weighting (default)
                - 'length_weighted': Length-proportional weighting
                - 'count': Equal weighting for all features
                - 'sum': Additive aggregation without averaging
                - 'frequency_weighted': Compositional weighting by type frequency
            reclass_map: Optional functional group classification.
                Dictionary mapping group names to substring lists.
                Example: {'Wetland': ['marsh', 'swamp'], 'Forest': ['woodland', 'trees']}
            new_class_field: Name for reclassified functional groups column
            group_by: Additional grouping variables for stratified analysis.
                Examples: ['management_zone', 'conservation_status']
            keep_unmatched_types: Whether to retain ecological types that
                don't match reclassification patterns
            weighting_field: Override default weighting field.
                Useful for custom ecological importance weights
            result_name: Storage key for aggregated results in target_scale.results

        Returns:
            pd.DataFrame: Aggregated ecological metrics without geometry.
                Columns include target identifiers, functional groups (if used),
                and aggregated metrics with method-specific suffixes.

        Raises:
            ValueError: If method unsupported or required columns missing

        Notes:
            - Preserves spatial relationships through proper weighting
            - Handles edge cases (no overlaps, missing data) gracefully
            - Provides detailed logging for methodological transparency
            - Results cached for computational efficiency
            - Maintains data provenance through clear column naming
        """
        # =====================================================================
        # INITIALIZE AGGREGATION PROCESS
        # =====================================================================
        # Create a working copy to avoid modifying the original joined data.
        # This ensures the spatial join cache remains clean for reuse and
        # prevents unintended side effects in subsequent analyses.
        df = df.copy()

        # =====================================================================
        # DATA QUALITY CONTROL: REMOVE DUPLICATE COLUMNS
        # =====================================================================
        # Spatial joins can create duplicate columns when source and target
        # scales have similarly named fields. Remove duplicates to prevent
        # aggregation errors and ensure clean results.
        if df.columns.duplicated().any():
            print(
                "\u26a0 Warning: Duplicate columns found in input DataFrame. Cleaning up."
            )
        df = df.loc[:, ~df.columns.duplicated()]

        # =====================================================================
        # CONFIGURE WEIGHTING AND METRICS
        # =====================================================================
        # Determine the weighting field to use for aggregation. This is critical
        # for ecologically-meaningful results:
        # - Area weights for habitat condition metrics
        # - Length weights for corridor connectivity measures
        # - Count weights for point-based observations
        weight_field = weighting_field or self.weighting_field
        metric_columns = [str(col) for col in (metric_columns or [])]

        # =====================================================================
        # FUNCTIONAL GROUP CLASSIFICATION (OPTIONAL)
        # =====================================================================
        # Apply ecological reclassification to group similar habitat types,
        # species, or management categories. This step is crucial for:
        # - Reducing complexity in diverse ecological datasets
        # - Creating management-relevant functional groups
        # - Enabling comparative analysis across regions
        df, reclass_field = self._apply_reclassification(
            df, reclass_map, keep_unmatched_types, new_class_field
        )

        # =====================================================================
        # CONSTRUCT GROUPING HIERARCHY
        # =====================================================================
        # Build the grouping structure for aggregation. The hierarchy typically
        # follows ecological organization:
        # 1. Target spatial units (catchments, management zones)
        # 2. Functional groups (wetland types, vegetation communities)
        # 3. Additional stratification (conservation status, management priority)

        group_keys = target_scale.unique_id_fields.copy()  # Base spatial units

        if reclass_field:  # Add functional groups if reclassification applied
            group_keys.append(reclass_field)

        if group_by:  # Add additional stratification variables
            group_keys += [g for g in group_by if g not in group_keys]

        # Report the final grouping structure for transparency
        print(f"\nAggregating by {group_keys} using '{method}' method")
        print(f"Processing {len(df)} spatially joined records")
        print(f"Target metrics: {metric_columns}")

        # =====================================================================
        # AGGREGATION METHOD 1: FREQUENCY-WEIGHTED ANALYSIS
        # =====================================================================
        # This method weights each ecological type by its relative frequency
        # within each target unit. Essential for compositional analysis where
        # the diversity and relative abundance of habitat types matters.
        #
        # Example: In a catchment with 60% wetlands and 40% forest,
        # wetland metrics get 0.6 weight, forest metrics get 0.4 weight.
        #
        # Mathematical formula:
        # Weight_i = Frequency_i / Total_frequency_in_target_unit
        if method == "frequency_weighted":
            if not reclass_field:
                raise ValueError("frequency_weighted requires a type or reclass field")

            freq = df.groupby(group_keys).size().rename("freq")
            total = freq.groupby(target_scale.unique_id_fields).sum().rename("total")
            weights = freq.to_frame().join(total, on=target_scale.unique_id_fields)
            weights["weight"] = weights["freq"] / weights["total"]
            df = df.join(weights["weight"], on=group_keys)

            weighted = df[metric_columns].multiply(df["weight"], axis=0)
            for key in group_keys:
                weighted[key] = df[key]
            result = weighted.groupby(group_keys).sum().reset_index()
            result.rename(
                columns={col: f"{col}_frqwm" for col in metric_columns}, inplace=True
            )

        # =====================================================================
        # AGGREGATION METHODS 2-4: UNIFIED WEIGHTED FRAMEWORK
        # =====================================================================
        # All remaining methods follow the same mathematical pattern:
        # Aggregated_Value = Σ(Metric × Weight) / Σ(Weight)
        #
        # This unified approach ensures consistency and allows for easy
        # extension to new weighting schemes. The key difference between
        # methods is how weights are calculated:
        #
        # - COUNT: Weight = 1 (equal contribution from each feature)
        # - SUM: Weight = 1, but no division (accumulative)
        # - WEIGHTED_MEAN: Weight = ecological importance (area/length/custom)
        elif method in {"count", "sum", "weighted_mean"}:

            # =================================================================
            # STEP 1: CONFIGURE METHOD-SPECIFIC WEIGHTS AND METRICS
            # =================================================================
            # Each aggregation method requires different weighting strategies
            # based on the ecological question being addressed.
            if method == "count":
                # COUNT METHOD: Equal weighting for feature counting
                # Use case: Number of monitoring sites, species occurrences,
                # management actions. Each feature contributes equally regardless
                # of size or other characteristics.
                df["Count"] = 1
                metrics_to_use = ["Count"]
                weight_to_use = "Count"
                suffix = "count"

            elif method == "sum":
                # SUM METHOD: Additive aggregation for extensive properties
                # Use case: Total habitat area, total species count, cumulative
                # impact scores. Values accumulate across space without averaging.
                if not metric_columns:
                    raise ValueError("sum method requires metric_columns")
                df["Sum"] = 1
                metrics_to_use = metric_columns
                weight_to_use = "Sum"
                suffix = "sum"

            elif method == "weighted_mean":
                # WEIGHTED MEAN METHOD: Ecologically-proportional weighting
                # Use case: Habitat condition scores, vegetation indices, water
                # quality measures. Larger/longer/more important features have
                # greater influence on aggregated values.
                if not metric_columns:
                    raise ValueError("weighted_mean method requires metric_columns")
                if not weight_field:
                    raise ValueError("weighted_mean requires a weighting_field")
                metrics_to_use = metric_columns
                weight_to_use = weight_field
                suffix = f"{weight_field[:3].lower()}wm"  # e.g., "arewm" or "lenwm"

            # =================================================================
            # STEP 2: APPLY UNIFIED WEIGHTED AGGREGATION FORMULA
            # =================================================================
            # Mathematical implementation of: Σ(Metric × Weight) / Σ(Weight)
            # This preserves the ecological meaning of intensive properties
            # while properly accounting for spatial extent or importance.
            weighted_metrics = df[metrics_to_use].multiply(df[weight_to_use], axis=0)
            weighted_metrics.columns = [f"{col}_weighted" for col in metrics_to_use]

            # Combine grouping keys, weights, and weighted metrics for aggregation
            agg_df = pd.concat(
                [
                    df[
                        group_keys + [weight_to_use]
                    ].copy(),  # .copy() prevents SettingWithCopyWarning
                    weighted_metrics,
                ],
                axis=1,
            )

            # =================================================================
            # STEP 3: SPATIAL AGGREGATION BY GROUPING HIERARCHY
            # =================================================================
            # Group by target units and functional groups, then sum weighted
            # values. This step combines all overlapping source features
            # within each target boundary.
            grouped = agg_df.groupby(group_keys).sum().reset_index()

            # =================================================================
            # STEP 4: COMPUTE FINAL AGGREGATED VALUES
            # =================================================================
            # Convert summed weighted values to final metrics based on the
            # chosen aggregation method. This step differs between methods.
            if method == "count":
                # COUNT RESULT: Total number of features in each target unit
                # Useful for: monitoring site density, species occurrence frequency
                result = grouped[group_keys + [weight_to_use]]
                result.rename(columns={weight_to_use: "FeatureCount"}, inplace=True)

            elif method == "sum":
                # SUM RESULT: Cumulative totals for extensive properties
                # Useful for: total habitat area, cumulative impact scores
                result = grouped[
                    group_keys + [f"{col}_weighted" for col in metrics_to_use]
                ]
                result.rename(
                    columns={
                        f"{col}_weighted": f"{col}_{suffix}" for col in metrics_to_use
                    },
                    inplace=True,
                )

            elif method == "weighted_mean":
                # WEIGHTED MEAN RESULT: Properly weighted averages
                # Useful for: habitat condition, vegetation health, water quality
                # Formula ensures larger/more important features have greater influence
                for col in metrics_to_use:
                    grouped[f"{col}_{suffix}"] = (
                        grouped[f"{col}_weighted"] / grouped[weight_to_use]
                    )
                result = grouped[
                    group_keys + [f"{col}_{suffix}" for col in metrics_to_use]
                ]

        else:
            raise ValueError(
                f"Unsupported method: {method}. Choose from: "
                f"count, sum, weighted_mean, frequency_weighted"
            )

        # =====================================================================
        # FINALIZE AND STORE AGGREGATION RESULTS
        # =====================================================================
        # Store results in the target scale's results dictionary for later
        # access, visualization, and export. Results are stored without
        # geometry to save memory and enable flexible output formats.

        result_obj = SpatialScale.Result(
            f"{self.name}_{method}_{target_scale.name}",
            target_scale,
            self.name,
            group_keys,
            result,
        )

        if result_name:
            # Initialize results dictionary if needed
            if not hasattr(target_scale, "results") or target_scale.results is None:
                target_scale.results = {}

            # Validate and potentially modify result name to avoid conflicts
            result_name = SpatialScale.validate_result_name(
                result_name, set(target_scale.results.keys())
            )

            result_obj.name = result_name
            # Store aggregated results for later use
            target_scale.results[result_name] = result_obj
            # = SpatialScale.Result(
            #     result_name,
            #     target_scale,
            #     self.name,
            #     group_keys,
            #     result,
            # )
            print(f"Stored aggregated results as '{result_name}'")

        # # create a copy of the stored result and add the geometry to return
        # geometry_df = target_scale.data[target_scale.unique_id_fields + ["geometry"]]
        # result_obj.data = geometry_df.merge(
        #     result_obj.data, on=target_scale.unique_id_fields, how="left"
        # )

        return result_obj

    def aggregate_to(
        self,
        target_scale: "SpatialScale",
        metric_columns: list[str],
        how: str = "intersects",
        method: str = "weighted_mean",
        weighting_field: Optional[str] = None,
        group_by: list[str] = None,
        reclass_map: dict[str, str] = None,
        new_class_field: str = "regrouped",
        keep_unmatched_types: bool = True,
        result_name: str = None,
    ) -> GeoDataFrame:
        """
        Spatially join and aggregate base scale metrics into the target scale.

        Aggregate ecological metrics from spatially joined data using
        ecologically-appropriate weighting schemes.

        This method implements the core scaling-up logic that transfers
        fine-resolution ecological data to broader management scales while
        preserving biological meaning and statistical validity.

        ECOLOGICAL RATIONALE:
        Different ecological metrics require different aggregation approaches:

        1. **Area-weighted means**: For intensive properties that represent
           conditions per unit area (e.g., species density, vegetation health,
           water quality). Larger habitats contribute proportionally more.

        2. **Length-weighted means**: For linear features like riparian zones
           where longer segments have greater ecological influence.

        3. **Count-based aggregation**: For point observations where each
           monitoring site contributes equally regardless of spatial extent.

        4. **Frequency-weighted means**: For compositional data where the
           relative abundance of different habitat types matters.

        5. **Simple sums**: For extensive properties that accumulate across
           space (e.g., total habitat area, total species count).

        MATHEMATICAL FRAMEWORK:
        All methods use a unified weighted aggregation formula:

        Aggregated_Value = Σ(Metric_i x Weight_i) / Σ(Weight_i)

        Where weights vary by method:
        - Area-weighted: Weight_i = Area_i (hectares)
        - Length-weighted: Weight_i = Length_i (meters)
        - Count-based: Weight_i = 1 (equal weighting)
        - Frequency-weighted: Weight_i = Frequency_i / Total_frequency

        Args:
            target_scale: SpatialScale object representing the aggregation
                target (e.g., catchments, management units)
            metric_columns: List of ecological metrics to aggregate.
                Examples: ['NDVI_mean', 'species_richness', 'pH', 'condition_score']
            method: Aggregation method selection:
                - 'area_weighted': Area-proportional weighting (default)
                - 'length_weighted': Length-proportional weighting
                - 'count': Equal weighting for all features
                - 'sum': Additive aggregation without averaging
                - 'frequency_weighted': Compositional weighting by type frequency
            how: Spatial join method (e.g., "intersects", "contains").
            weighting_field: Override default weighting field.
                Useful for custom ecological importance weights
            group_by: Additional grouping variables for stratified analysis.
                Examples: ['management_zone', 'conservation_status']
            reclass_map: Optional functional group classification.
                Dictionary mapping group names to substring lists.
                Example: {'Wetland': ['marsh', 'swamp'], 'Forest': ['woodland', 'trees']}
            new_class_field: Name for reclassified functional groups column
            keep_unmatched_types: Whether to retain ecological types that
                don't match reclassification patterns
            result_name: Storage key for aggregated results in target_scale.results



        Returns:
            Result: Result object containing aggregated ecological data with
                metadata about the scaling process. Access the GeoDataFrame
                via result.data attribute.


        Examples:
            >>> # Area-weighted wetland condition aggregation
            >>> result = wetlands.aggregate_to(
            ...     target_scale=catchments,
            ...     metric_columns=['condition_score', 'NDVI_2023'],
            ...     method='area_weighted',
            ...     reclass_map={'Permanent': ['permanent', 'deep'],
            ...                  'Seasonal': ['seasonal', 'ephemeral']},
            ...     result_name='wetland_condition'
            ... )

            >>> # Frequency-weighted habitat diversity
            >>> result = habitats.aggregate_to(
            ...     target_scale=planning_units,
            ...     metric_columns=['diversity_index'],
            ...     method='frequency_weighted',
            ...     group_by=['conservation_status']
            ... )

        Notes:
            - Preserves spatial relationships through proper weighting
            - Handles edge cases (no overlaps, missing data) gracefully
            - Provides detailed logging for methodological transparency
            - Results cached for computational efficiency
            - Maintains data provenance through clear column naming


        """
        if not isinstance(target_scale, SpatialScale):
            raise TypeError("target_scale must be a SpatialScale instance")

        if not target_scale.source.exists():
            raise FileNotFoundError(
                f"Target scale file not found: {target_scale.source}"
            )
        # Semantic method aliases
        semantic_aliases = {
            "area_weighted": ("weighted_mean", "Area_Ha"),
            "length_weighted": ("weighted_mean", "Length_m"),
        }

        # Translate semantic alias to real method and weighting field
        if method in semantic_aliases:
            actual_method, default_weight = semantic_aliases[method]
            method = actual_method
            if not weighting_field:
                weighting_field = default_weight
                print(f"Using '{weighting_field}' as weight for '{method}' aggregation")
        elif weighting_field:
            print(
                f"Using '{weighting_field}' for '{weighting_field}_weighted' aggregation"
            )
        if not metric_columns:
            raise ValueError("metric_columns cannot be empty")

        metric_columns = [str(col) for col in metric_columns]

        # Validate metric columns in this scale
        self._validate_metric_data(
            self.data,
            unique_id=self.unique_id_fields,
            metric_field_names=metric_columns,
            not_allowed=[self.type_field, self.weighting_field] + (group_by or []),
        )

        # Check result name validity before proceeding
        if result_name:
            existing_keys = getattr(target_scale, "results", {}).keys()
            result_name = SpatialScale.validate_result_name(
                result_name, existing_keys=set(existing_keys)
            )

        # Spatial join
        joined_df = self.spatial_join(target_scale, how=how)

        # Aggregate and store in target scale
        return self._aggregate_joined(
            joined_df,
            target_scale,
            metric_columns=metric_columns,
            method=method,
            reclass_map=reclass_map,
            new_class_field=new_class_field,
            group_by=group_by,
            keep_unmatched_types=keep_unmatched_types,
            weighting_field=weighting_field,
            result_name=result_name,
        )

    def save_results(
        self,
        output_folder: Union[str, Path],
        file_types: Union[str, list[str]] = ".csv",
    ) -> None:
        """
        Save the aggregated results from this SpatialScale to disk.

        Results are stored internally without geometry and rejoined with
        geometry from self.data before saving.

        Args:
            output_folder: Directory where files will be saved.
            file_types: File types to save (e.g., ".csv", ".shp", ".gpkg").

        Raises:
            ValueError: If no results exist or unsupported file types provided.
        """
        if not self.results:
            raise ValueError("No results to save. Run `aggregate_to` first.")

        if output_folder == "list":
            print("Available results:")
            for result_name, result_obj in self.results.items():
                # group_by = [f for f in group_keys if f not in self.unique_id_fields]
                print(
                    f"  {result_name.ljust(30)} - scaled from: '{result_obj.from_scale_name}' to {result_obj.parent.name} with grouping: {result_obj.group_by} - {str(len(result_obj.data)).rjust(10)} rows {result_obj.data.columns.tolist()}"
                )
            return

        for result_name, result_obj in self.results.items():
            result_obj.export(output_folder, file_types)

    def _validate_metric_data(
        self,
        df: pd.DataFrame,
        unique_id: Union[str, List[str]] = None,
        metric_field_names: Optional[Union[str, list]] = None,
        not_allowed: Optional[list] = None,
        strict: bool = False,
        set_index: bool = True,
    ) -> pd.DataFrame:
        """
        Validate that metric fields in a dataframe are numeric and do not conflict
        with SpatialScale internal fields like unique_id, weighting_field, or type_field.

        Args:
            df: Input DataFrame or GeoDataFrame.
            metric_field_names: Metric fields to validate. If None, infer from numeric fields.
            group_by: Optional grouping fields to exclude from metrics.
            strict: If True, raise if non-numeric or missing values are found.
            set_index: If True, set index to self.unique_id_fields.

        Returns:
            Cleaned DataFrame or GeoDataFrame with validated numeric metrics.
        """
        df = df.copy()
        is_geo = isinstance(df, GeoDataFrame)

        # Normalize columns
        df.columns = df.columns.map(str)
        unique_id = (
            [unique_id]
            if isinstance(unique_id, str)
            else unique_id or self.unique_id_fields
        )
        reserved = {
            # unique_id,
            self.weighting_field,
            self.type_field,
            *(not_allowed or []),
        }
        reserved.update(unique_id)

        (
            [metric_field_names]
            if isinstance(metric_field_names, str)
            else metric_field_names
        )

        if not metric_field_names:
            metric_field_names = [
                col
                for col in df.columns
                if col not in reserved and pd.api.types.is_numeric_dtype(df[col])
            ]
            print(f"No metric fields specified. Inferred: {metric_field_names}")
        else:
            if isinstance(metric_field_names, str):
                metric_field_names = [metric_field_names]
            metric_field_names = [str(c) for c in metric_field_names]

        # Check for conflicts
        conflicts = [col for col in metric_field_names if col in reserved]
        if conflicts:
            raise ValueError(
                f"Metric fields conflict with reserved fields: {conflicts}"
            )

        # Check column presence
        required_cols = unique_id + metric_field_names
        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            available_cols = [col for col in df.columns if col != "geometry"]
            raise ValueError(
                f"Missing required columns: {missing}. Available columns: {available_cols}"
            )

        # Enforce numeric
        df[metric_field_names] = df[metric_field_names].apply(
            pd.to_numeric, errors="coerce"
        )
        invalid_rows = df[required_cols].isna().any(axis=1)

        if invalid_rows.any():
            count = invalid_rows.sum()
            if strict:
                raise ValueError(f"{count} invalid rows found.")
            df = df[~invalid_rows]
            print(
                f"\u26a0 Dropping {count} invalid rows (missing values or not numeric). {len(df)} rows remain."
            )

        # Enforce uniqueness of unique_id (can be multiple columns)
        duplicated_mask = df.duplicated(subset=unique_id, keep=False)
        if duplicated_mask.any():
            dup_rows = df[duplicated_mask]
            dup_keys = dup_rows[unique_id].drop_duplicates()

            # Format each duplicate key nicely
            formatted_keys = "\n".join(
                "- " + ", ".join(str(val) for val in row) for row in dup_keys.values
            )

            msg = (
                f"Duplicate compound keys found in columns: {unique_id}.\n"
                f"Each combination of values must be unique across rows (order matters).\n"
                f"Examples of duplicates:\n{formatted_keys}"
            )
            print(msg)
            raise ValueError(msg)

        result = df.set_index(unique_id) if set_index else df
        return (
            GeoDataFrame(result, geometry="geometry", crs=df.crs) if is_geo else result
        )

    def standardise(
        self,
        fields: Optional[list[str]] = None,
        suffix: str = "z",
    ) -> pd.DataFrame:
        """
        Standardise selected numeric fields (z-score) from a stored result.
        The z-score is the number of standard deviations a data point is away from the mean

        Args:
            fields: List of columns to standardise (default: all numeric excluding protected).
            suffix: Suffix to append to new standardised columns (e.g., "z" ? "2020z").

        Returns:
            DataFrame (or GeoDataFrame) with standardised fields.
        """

        protected = {self.unique_id_fields, self.weighting_field, self.type_field}

        if fields:
            fields = [str(f) for f in fields]
            invalid = [f for f in fields if f not in self.data.columns]
            if invalid:
                raise ValueError(f"Fields not found in result: {invalid}")
        else:
            # Auto-detect numeric fields not protected
            fields = [
                col
                for col in self.data.columns
                if col not in protected
                and pd.api.types.is_numeric_dtype(self.data[col])
            ]

        print(f"Standardizing fields: {fields} (suffix: '{suffix}')")

        # Perform z-score standardization
        for col in fields:
            self.data[f"{col}{suffix}"] = (
                self.data[col] - self.data[col].mean()
            ) / self.data[col].std()

        print(f"Standardized fields are: {[f'{col}{suffix}' for col in fields]}")

        return self.data

    def normalise(
        self,
        fields: Optional[list[str]] = None,
        suffix: str = "n",
    ) -> pd.DataFrame:
        """
        Normalize selected numeric fields (min-max scaling) from a stored result.
        i.e. all fields will be rescaled to between 0 and 1.

        Args:
            fields: List of columns to normalise (default: all numeric excluding protected).
            suffix: Suffix for normalised output columns (e.g., "n" ? "2020n").

        Returns:
            DataFrame (or GeoDataFrame) with normalised fields.
        """

        protected = {self.unique_id_fields, self.weighting_field, self.type_field}

        if fields:
            fields = [str(f) for f in fields]
            invalid = [f for f in fields if f not in self.data.columns]
            if invalid:
                raise ValueError(f"Fields not found in result: {invalid}")
        else:
            fields = [
                col
                for col in self.data.columns
                if col not in protected
                and pd.api.types.is_numeric_dtype(self.data[col])
            ]

        print(f"Normalising fields: {fields} (suffix: '{suffix}')")

        for col in fields:
            min_val = self.data[col].min()
            max_val = self.data[col].max()
            if min_val == max_val:
                print(f"\u26a0 Skipping {col} no variation (min = max = {min_val})")
                continue
            self.data[f"{col}{suffix}"] = (self.data[col] - min_val) / (
                max_val - min_val
            )
        print(f"Normalised fields are: {[f'{col}{suffix}' for col in fields]}")

        return self.data

    def plot(
        self,
        result: str = None,
        field: str = None,
        reclass_map: dict = None,
        reclass_field: str = "regrouped",  # this is the default can be over-ridden when applying the reclass_map
        vmin: float = None,
        vmax: float = None,
        title: str = None,
    ) -> None:
        """
        Plot the base scale spatial data.

        Args:
            base_scale: SpatialScale object containing the base data
            title: Optional title for the plot
        """

        print(f"Plotting {self.name}...")
        if result:
            # result = result.lower()
            if result not in self.results.keys():
                print(
                    f"Result not found. Using {self.name} data. Saved results are {list(self.results.keys())}"
                )
                df = self.data
            else:
                print(f"Using save result: {result}")
                result_obj = self.results[result]
                # results are saved without geometry for efficiency.  So now re-attach the geometry
                geometry_df = self.data[self.unique_id_fields + ["geometry"]]
                df = geometry_df.merge(
                    result_obj.data, on=self.unique_id_fields, how="left"
                )
            result_grouping = result_obj.group_by
        else:
            df = self.data

        if reclass_map:
            df, reclass_field = self._apply_reclassification(
                df, reclass_map, keep_unmatched_types=False
            )

        if field is None or field not in df.columns:
            if field is not None:
                print(
                    f"Field '{field}' not found in data. "
                    f"Available fields: {df.columns.tolist()}"
                )
            ax = self.data.plot(color="lightblue", edgecolor="blue")
            ax.set_title(f"Base Scale: {self.name}")
            ax.axis("off")
        else:
            # Color scale settings
            vmin = vmin or df[field].min()
            vmax = vmax or df[field].max()
            cmap = plt.cm.viridis
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

            # Group info
            if reclass_field and reclass_field in df.columns:
                grp = df[reclass_field]
            elif result_grouping:
                # No grouping field available, create single group
                # grp = pd.Series(['All'] * len(df), index=df.index)
                type_field = [
                    f for f in result_grouping if f not in self.unique_id_fields
                ]
                grp = df[type_field[0]]
                print(f"grouped by: {type_field}")
            elif self.type_field and self.type_field in df.columns:
                grp = df[self.type_field]
            else:
                raise ValueError(
                    f"Grouping field not found. Available fields: {df.columns.tolist()}"
                )

            # groups = grp.dropna().unique() if grp is not None else ["All"]
            groups = grp.unique()
            n_groups = len(groups)
            cols = 3
            rows = math.ceil(n_groups / cols)

            # Create subplots
            fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5 * rows))
            axes = axes.flatten()

            # Store axis and group info for title-setting step
            plotted_axes = []

            # Plot each group
            for i, group in enumerate(groups):
                ax = axes[i]

                # Set axis to base extent
                bounds = df.total_bounds
                ax.set_xlim(bounds[0], bounds[2])
                ax.set_ylim(bounds[1], bounds[3])

                filtered_result = df[grp == group]
                filtered_result.plot(
                    column=field,
                    cmap=cmap,
                    norm=norm,
                    ax=ax,
                    edgecolor="blue",
                    legend=False,
                )
                ax.set_title(f"{self.name} - {group}")
                # ax.axis("off")
                # keep the axes as a bounding box but hide the tickmarks and labels
                ax.tick_params(
                    axis="both",
                    which="both",
                    length=0,
                    labelbottom=False,
                    labelleft=False,
                )
                plotted_axes.append((ax, group))

            # Turn off unused axes
            for j in range(i + 1, len(axes)):
                axes[j].axis("off")

            # Shared vertical colorbar next to first row
            first_row_axes = axes[0:cols]
            bbox = first_row_axes[-1].get_position()
            cbar_ax = fig.add_axes([bbox.x1 + 0.01, bbox.y0, 0.01, bbox.y1 - bbox.y0])
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm._A = []
            cbar = fig.colorbar(sm, cax=cbar_ax)
            cbar.set_label(field)

        plt.show()

    class Result:
        """
        Container for aggregated ecological data from spatial scaling operations.

        This class encapsulates the results of hierarchical spatial aggregation,
        providing metadata about the scaling process and source data. It serves
        as the return object for SpatialScale.aggregate_to() operations.

        The Result class maintains provenance information and resulting metrics.

        Attributes:
            name (str): Identifier for this aggregation result
            parent: the SpatialScale object holding the result (the target scale)
            from_scale_name (str): Name of the source spatial scale
            group_by (List[str]): grouping variables used to stratify aggregation
            data (pd.DataFrame): Aggregated ecological metrics without geometries

        """

        def __init__(
            self,
            name: str,
            parent: "SpatialScale",
            from_scale_name: str,
            group_by: List[str],
            data: pd.DataFrame,
        ):

            self.name = name
            self.parent = parent
            self.from_scale_name = from_scale_name
            self.group_by = group_by
            self.data = data

        def __str__(self) -> str:
            """
            Generate human-readable summary of aggregation results.

            Returns:
                str: Formatted summary including metadata and data preview
            """
            return (
                f"{'-'*5} Result: {self.name} {'-'*20}\n"
                f"Scaled from: {self.from_scale_name}\n"
                f"Scaled to: {self.parent.name}\n"
                f"Grouped by: {self.group_by}\n"
                f"Records: {len(self.data)}\n"
                f"Columns: {list(self.data.columns)}\n"
                f"Data preview:\n{self.data.head()}\n"
            )

        def summary_stats(self) -> pd.DataFrame:
            """
            Generate summary statistics for aggregated ecological metrics.

            Returns:
                pd.DataFrame: Statistical summary of numeric columns
            """
            numeric_cols = self.data.select_dtypes(include=[np.number]).columns
            return self.data[numeric_cols].describe()

        def export(
            self, output_folder: Union[str, Path], file_format: str = "csv"
        ) -> None:
            """
            Export aggregated results to file.

            Args:
                output_folder: Directory where files will be saved.
                file_types: File types to save (e.g., ".csv", ".shp", ".gpkg").
            """
            if output_folder is None:
                print("No output path provided. Results not saved.")
                return

            # Normalize and prepare output folder
            output_folder = Path.cwd() / Path(output_folder)
            output_folder.mkdir(parents=True, exist_ok=True)

            # Normalize file types
            file_format = [file_format] if isinstance(file_format, str) else file_format
            file_format = [
                f.lower() if f.startswith(".") else f".{f.lower()}" for f in file_format
            ]

            valid_types = {".csv", ".shp", ".gpkg"}
            for ext in file_format:
                if ext not in valid_types:
                    raise ValueError(
                        f"Invalid file format: {ext}. Valid types: {valid_types}"
                    )
                out_file = output_folder / f"{self.name}{ext}"
                print(f"exporting {out_file}")
                if ext == ".csv":
                    # Export without geometry for CSV
                    # data_no_geom = self.data.drop(columns="geometry", errors="ignore")
                    self.data.to_csv(out_file, index=False)
                if ext in [".gpkg", ".shp"]:

                    geometry_df = self.parent.data[
                        self.parent.unique_id_fields + ["geometry"]
                    ]
                    result_df = geometry_df.merge(
                        self.data, on=self.parent.unique_id_fields, how="left"
                    )
                    if ext == f".gpkg":
                        result_df.to_file(out_file, driver="GPKG")
                    if ext == f".shp":
                        result_df.to_file(out_file, driver="ESRI Shapefile")


def plot_spatial_hierarchy(*args: SpatialScale) -> None:
    """
    Plot all spatial scales together with thumbnails and overlay plot.

    Args:
        *args: SpatialScale objects where first is base scale, rest are aggregation scales

    Example:
        plot_spatial_hierarchy(base_scale, region_scale, basin_scale)
    """
    if not args:
        raise ValueError("At least one SpatialScale object required")

    base_scale = args[0]
    agg_scale_objects = list(args[1:])

    print(
        f"Plotting {base_scale.name} with {len(agg_scale_objects)} aggregation scales..."
    )

    def draw_scale(
        ax, gdf, title=None, facecolor="none", edgecolor="black", lw=1, alpha=0.7
    ):
        gdf.plot(ax=ax, color=facecolor, edgecolor=edgecolor, lw=lw, alpha=alpha)
        if title:
            ax.set_title(title, fontsize=10)
        ax.axis("off")

    # Prepare color map for aggregation scales
    n_scales = len(agg_scale_objects)
    colors = cm.Set1(np.linspace(0, 1, n_scales)) if n_scales > 0 else []

    # --------- Thumbnails ----------
    all_scales = [(base_scale.name, base_scale, "lightblue", "blue", 0.7)]
    for i, scale_obj in enumerate(agg_scale_objects):
        all_scales.append((scale_obj.name, scale_obj, "none", colors[i], 1))

    n_thumbs = len(all_scales)
    fig, axes = plt.subplots(1, n_thumbs, figsize=(2.5 * n_thumbs, 1.25))
    if n_thumbs == 1:
        axes = [axes]

    for ax, (name, scale_obj, facecolor, edgecolor, lw) in zip(axes, all_scales):
        draw_scale(
            ax,
            scale_obj.data,
            title=name,
            facecolor=facecolor,
            edgecolor=edgecolor,
            lw=lw,
        )

    plt.tight_layout()
    plt.show()

    # --------- Main Overlay Plot ----------
    fig, ax = plt.subplots(figsize=(10, 8))

    # Base layer
    draw_scale(
        ax,
        base_scale.data,
        facecolor="lightblue",
        edgecolor="blue",
        lw=0.5,
        alpha=0.6,
    )

    # Aggregation layers
    for i, scale_obj in enumerate(agg_scale_objects):
        draw_scale(
            ax,
            scale_obj.data,
            facecolor="none",
            edgecolor=colors[i],
            lw=2,
            alpha=1.0,
        )

    # Set axis to base extent
    bounds = base_scale.data.total_bounds
    ax.set_xlim(bounds[0], bounds[2])
    ax.set_ylim(bounds[1], bounds[3])

    # Legend
    legend_elements = [Line2D([0], [0], color="lightblue", lw=4, label=base_scale.name)]
    for i, scale_obj in enumerate(agg_scale_objects):
        legend_elements.append(
            Line2D([0], [0], color=colors[i], lw=2, label=scale_obj.name)
        )

    ax.legend(loc="upper left", handles=legend_elements)
    ax.set_title("Spatial Scale Hierarchy")
    plt.show()
