"""
MDB Scaling-up Framework - Core Functions

This file contains all the main functions that do the actual work of:
1. Loading and validating your spatial data
2. Calculating appropriate weights (area, length, or count)
3. Aggregating small features into larger management units
4. Handling different geometry types (points, lines, polygons)
5. Managing functional groups and classifications

For novice users:
- You typically don't need to modify this file
- All settings are controlled through config.py
- The functions here are called automatically by the notebook
- Each function has detailed comments explaining what it does
"""

import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import re
from matplotlib.lines import Line2D
from typing import List, Optional, Union
from pathlib import Path
from geopandas import GeoDataFrame
from shapely.validation import make_valid, explain_validity
from pyproj import CRS


class SpatialScale:
    """
    A class that represents one spatial scale in your analysis.

    Think of this as a "container" that holds:
    - Your spatial data (points, lines, or polygons)
    - Information about how to measure each feature (area, length, count)
    - Rules for grouping similar features together

    This class handles all the technical details of:
    - Loading shapefiles
    - Detecting geometry types (point/line/polygon)
    - Calculating appropriate measures (area in hectares, length in meters, etc.)
    - Validating that your data is properly formatted

    For novice users:
    - One SpatialScale = one shapefile in your analysis
    - The "base scale" contains your smallest features
    - "Aggregation scales" contain larger boundaries
    - You create these through the config.py settings
    """

    def __init__(
        self,
        name: str,
        source: Union[str, Path],
        unique_id_field: Union[str, List[str]],
        weighting_field: Optional[str] = None,
        metric_fields: Optional[List[str]] = None,
        measure_multiplier: Optional[float] = None,
        type_field: Optional[str] = None,
        default_crs: str = None,
        is_base_scale: bool = False,
    ):
        """
        Create a new spatial scale for analysis.

        This function loads your shapefile and prepares it for analysis by:
        1. Reading the geographic data
        2. Detecting what type of geometry it contains (points/lines/polygons)
        3. Calculating appropriate measures (area/length/count)
        4. Validating that everything looks correct

        Parameters explained in plain English:
            name: A short name for this scale (used in output files)
            source: Path to your shapefile
            unique_id_field: Column that uniquely identifies each feature
            weighting_field: Column containing size/importance measure (or None to auto-calculate)
            metric_fields: List of columns containing measures to be used for analysis
            measure_multiplier: Number to multiply measures by (usually leave as None)
            type_field: Column containing feature types for grouping (or None for no grouping)
            is_base_scale: True if this is your smallest scale, False for larger scales
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

        # Store unique_id as string (convert from list if needed)
        self.unique_id_field = (
            unique_id_field if isinstance(unique_id_field, str) else unique_id_field[0]
        )

        # -----------------------------------------------------
        # Determine which columns are required from the shapefile
        # Only load what we need to make processing faster
        # -----------------------------------------------------
        required_cols = {self.unique_id_field}
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

        # Handle CRS - default to EPSG:4326 if not specified or invalid
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

        # Read and validate shapefile
        self.data = self._validate_geometries(
            gpd.read_file(self.source).to_crs(self.default_crs), fix_invalid=True
        )

        # Convert columns to strings (if not already)
        self.data.columns = self.data.columns.map(str)

        # Convert required_cols to strings to match DataFrame columns
        required_cols = [str(c) for c in required_cols]

        # Select columns based on whether metric_fields is specified
        if metric_fields:
            # Only keep required columns when metric_fields is specified
            cols_to_keep = list(set(required_cols) | {"geometry"})
            self.data = self.data.loc[:, cols_to_keep]
        print(f"   No specific metric fields chosen so loading all columns.")
        cols = [col for col in self.data.columns if col != "geometry"]
        print(f"   Loaded: {cols}")
        # If metric_fields is None, retain all columns (no filtering)

        # Handle complex geometries by breaking them into simple parts
        # E.g., a MultiPolygon becomes multiple separate Polygons
        # This ensures each feature is treated as a separate unit for analysis
        multi_geom_count = self.data.geometry.geom_type.isin(
            ["MultiPolygon", "MultiLineString", "MultiPoint"]
        ).sum()
        if multi_geom_count > 0:
            print(
                f"   Breaking {multi_geom_count} complex geometries into simple parts..."
            )
            self.data = self.data.explode(index_parts=False).reset_index(drop=True)

        # Check that all required columns exist in the shapefile
        missing_cols = [c for c in required_cols if c not in self.data.columns]
        if missing_cols:
            raise ValueError(
                f"Required columns {missing_cols} not found in {source}\n"
                f"Available columns: {list(self.data.columns)}\n"
                f"Please check your column names in config.py"
            )

        # Automatically detect what type of geometry this shapefile contains
        # This determines how we calculate measures and perform aggregation
        first_geom = self.data.geometry.iloc[0]
        if first_geom.geom_type == "Polygon":
            geometry_type = "polygon"  # Areas (wetlands, habitats, management units)
        elif first_geom.geom_type == "LineString":
            geometry_type = "line"  # Linear features (rivers, transects, corridors)
        elif first_geom.geom_type == "Point":
            geometry_type = "point"  # Point locations (monitoring sites, observations)
        else:
            raise ValueError(f"Unsupported geometry type: {first_geom.geom_type}")

        print(f"   Detected geometry type: {geometry_type} ({len(self.data)} features)")

        # Calculate default weighting field based on geometry type
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
        tabular_data_unique_id_field: Optional[str] = None,
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
            unique_id_field: Column in the CSV matching source unique ID. If None, defaults to self.unique_id.
            pivot_index: If pivoting, list of columns to use as index.
            pivot_columns: If pivoting, column name whose values become new columns.
            pivot_values: If pivoting, column name containing values to fill.
            how: Type of join to perform (default "left").

        Raises:
            ValueError if required columns not found.
        """

        # Load CSV data
        tabular_data = pd.read_csv(data_path)

        # Determine the ID field to use
        if tabular_data_unique_id_field is None:
            tabular_data_unique_id_field = self.unique_id_field

        pivot_row_id = [pivot_row_id] if isinstance(pivot_row_id, str) else pivot_row_id

        # Check unique_id field exists
        if tabular_data_unique_id_field not in tabular_data.columns:
            raise ValueError(
                f"unique_id_field '{tabular_data_unique_id_field}' not found in data columns: {list(tabular_data.columns)}"
            )
        if (
            not pivot_row_id
            and tabular_data[tabular_data_unique_id_field].duplicated().any()
        ):
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
            self.unique_id_field,
            not_allowed=[self.weighting_field, self.type_field],
        )

        # Single merge operation handles both same and different column names
        self.data = self.data.merge(
            tabular_data_clean,
            how=how,
            left_on=self.unique_id_field,
            right_on=tabular_data_unique_id_field,
            suffixes=("", "_external"),
            validate="one_to_one",
        )

        print(
            f"Joined external data from {data_path} on '{tabular_data_unique_id_field}'"
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

            # Update target's unique_id field if it was renamed
            updated_uid = rename_map.get(target_scale.unique_id_field)
            if updated_uid:
                target_scale.unique_id_field = updated_uid
                print(
                    f"\u26a0 Warning: Updated target scale unique_id_field to '{updated_uid}'"
                )

            # Run the spatial join (inner join only, to maintain overlaps)
            joined = gpd.sjoin(self.data, target_scale.data, how="inner", predicate=how)
            self._join_cache[cache_key] = joined

        else:
            print(f"Using cached spatial join: {cache_key}")

        return self._join_cache[cache_key]

    def aggregate_joined(
        self,
        df: GeoDataFrame,
        target_scale: "SpatialScale",
        metric_columns: list[str],
        method: str = "weighted_mean",
        reclass_map: dict[str, str] = None,
        group_by: list[str] = None,
        result_name: str = None,
        keep_unmatched: bool = True,
        weighting_field: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Aggregate metrics from a spatial join result into the target spatial scale.

        All aggregation types are handled using a unified weighted framework:
        - "count": metric=1, weight=1
        - "sum": metric=values, weight=1
        - "weighted_mean": metric=values, weight=from field
        - "frequency_weighted": metric=values, weight=frequency / total per group

        Args:
            df: Joined GeoDataFrame of source + target features.
            target_scale: The SpatialScale to aggregate into.
            metric_columns: List of fields to aggregate.
            method: Aggregation method.
            reclass_map: Optional mapping from types to groups.
            group_by: Optional extra grouping fields.
            result_name: Key to store result in target_scale.results.
            keep_unmatched: Whether to retain unmapped types if reclass_map is provided.
            weighting_field: Optional override of default weight field.

        Returns:
            DataFrame: Aggregated result (without geometry).
        """
        df = df.copy()
        weight_field = weighting_field or self.weighting_field
        metric_columns = [str(col) for col in (metric_columns or [])]

        # --- Handle type reclassification ---
        reclass_field = None
        if reclass_map and self.type_field:
            print(
                f"Reclassifying '{self.type_field}' into new groups '_grp_' using supplied reclass map..."
            )
            flat_map = {
                t: g
                for g, types in reclass_map.items()
                for t in (types if isinstance(types, (list, set, tuple)) else [types])
            }
            df["_grp_"] = df[self.type_field].map(flat_map)
            unmatched = df[df["_grp_"].isna()][self.type_field].unique()
            if len(unmatched) > 0:
                if not keep_unmatched:
                    df = df[~df["_grp_"].isna()].copy()
                    print(f"Dropped {len(unmatched)} unmatched types")
                else:
                    print(f"Retained {len(unmatched)} unmatched types")
            df["_grp_"] = df["_grp_"].fillna(df[self.type_field])
            reclass_field = "_grp_"
        elif self.type_field:
            reclass_field = self.type_field

        # --- Grouping keys ---
        group_keys = [target_scale.unique_id_field]
        if reclass_field:
            group_keys.append(reclass_field)
        if group_by:
            group_keys += [g for g in group_by if g not in group_keys]

        print(f"Aggregating by {group_keys} using '{method}'")

        # --- Frequency-weighted method ---
        if method == "frequency_weighted":
            if not reclass_field:
                raise ValueError("frequency_weighted requires a type or reclass field")

            freq = df.groupby(group_keys).size().rename("freq")
            total = freq.groupby(target_scale.unique_id_field).sum().rename("total")
            weights = freq.to_frame().join(total, on=target_scale.unique_id_field)
            weights["weight"] = weights["freq"] / weights["total"]
            df = df.join(weights["weight"], on=group_keys)

            weighted = df[metric_columns].multiply(df["weight"], axis=0)
            for key in group_keys:
                weighted[key] = df[key]
            result = weighted.groupby(group_keys).sum().reset_index()
            result.rename(
                columns={col: f"{col}_frqwm" for col in metric_columns}, inplace=True
            )

        # --- Unified aggregation framework ---
        # All methods follow the same pattern: metric * weight, then aggregate
        # This leverages the fact that __init__ creates appropriate weight fields:
        # - Count method: uses Count=1 (from aggregation_method="count")
        # - Sum method: uses Sum=1 (from aggregation_method="sum")
        # - Geometry method: uses Area_Ha, Length_m, or Count based on geometry type
        elif method in {"count", "sum", "weighted_mean"}:

            # Step 1: Determine metrics and weights based on method
            if method == "count":
                # Count aggregation: each feature = 1, weight = 1 (simple counting)
                # Create Count field on-demand for this aggregation
                df["Count"] = 1
                metrics_to_use = ["Count"]
                weight_to_use = "Count"
                suffix = "count"

            elif method == "sum":
                # Sum aggregation: metric values, weight = 1 (unweighted sum)
                # Create Sum field on-demand for this aggregation
                if not metric_columns:
                    raise ValueError("sum method requires metric_columns")
                df["Sum"] = 1
                metrics_to_use = metric_columns
                weight_to_use = "Sum"
                suffix = "sum"

            elif method == "weighted_mean":
                # Weighted mean: metric values, weight = area/length/count from __init__
                # Uses Area_Ha, Length_m, or Count created in __init__ based on geometry type
                if not metric_columns:
                    raise ValueError("weighted_mean method requires metric_columns")
                if not weight_field:
                    raise ValueError("weighted_mean requires a weighting_field")
                metrics_to_use = metric_columns
                weight_to_use = weight_field
                suffix = f"{weight_field[:4]}wm"  # e.g., "Area_wm" or "Leng_wm"

            # Step 2: Apply the unified weighted aggregation formula
            # Formula: (metric1*weight + metric2*weight + ...) / (weight1 + weight2 + ...)
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

            # Step 3: Group and sum all weighted values
            grouped = agg_df.groupby(group_keys).sum().reset_index()

            # Step 4: Calculate final results based on method
            if method == "count":
                # For count: just return the sum of weights (each weight=1, so this gives count)
                result = grouped[group_keys + [weight_to_use]]
                result.rename(columns={weight_to_use: "FeatureCount"}, inplace=True)

            elif method == "sum":
                # For sum: return the weighted sums (weight=1, so this gives unweighted sums)
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
                # For weighted mean: divide weighted sums by weight sums to get averages
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

        # --- Store results (without geometry) ---
        if not hasattr(target_scale, "results") or target_scale.results is None:
            target_scale.results = {}

        result_name = SpatialScale.validate_result_name(
            result_name, set(target_scale.results.keys())
        )
        target_scale.results[result_name] = result
        print(f"Stored result in {target_scale.name}.results['{result_name}']")
        return result

    def aggregate_to(
        self,
        target_scale: "SpatialScale",
        metric_columns: list[str],
        method: str = "weighted_mean",
        reclass_map: dict[str, str] = None,
        group_by: list[str] = None,
        result_name: str = None,
        keep_unmatched: bool = True,
        how: str = "intersects",
        weighting_field: Optional[str] = None,
    ) -> GeoDataFrame:
        """
        Spatially join and aggregate base scale metrics into the target scale.

        This wraps the spatial join and aggregation into one call.

        Args:
            target_scale: The SpatialScale to aggregate into.
            metric_columns: List of numeric fields to aggregate.
            method: Aggregation method ("sum", "mean", "weighted_mean", etc.).
            reclass_map: Optional mapping for grouping types (used if self.type_field exists).
            group_by: Optional extra fields to group by (e.g., basin, catchment).
            result_name: Key to store result in target_scale.results.
            keep_unmatched: If False, drops types not in reclass_map.
            how: Spatial join method (e.g., "intersects", "contains").
            weighting_field: Optional override of this scale's weighting field.

        Returns:
            GeoDataFrame: Aggregated result with geometries from the target scale.
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
            unique_id=self.unique_id_field,
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
        return self.aggregate_joined(
            joined_df,
            target_scale,
            metric_columns=metric_columns,
            method=method,
            reclass_map=reclass_map,
            group_by=group_by,
            result_name=result_name,
            keep_unmatched=keep_unmatched,
            weighting_field=weighting_field,
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
            for result_name, result_df in self.results.items():
                print(
                    f"  {result_name.ljust(30)} - {str(len(result_df)).rjust(10)} rows {result_df.columns.tolist()}"
                )
            return

        if output_folder is None:
            print("No output path provided. Results not saved.")
            return

        # Normalize and prepare output folder
        output_folder = Path.cwd() / Path(output_folder).name
        output_folder.mkdir(parents=True, exist_ok=True)

        # Normalize file types
        file_types = [file_types] if isinstance(file_types, str) else file_types
        file_types = [
            f.lower() if f.startswith(".") else f".{f.lower()}" for f in file_types
        ]
        valid_types = {".csv", ".shp", ".gpkg"}

        for ext in file_types:
            if ext not in valid_types:
                raise ValueError(
                    f"Invalid file type: {ext}. Valid types: {', '.join(valid_types)}"
                )

        # Get base geometry for join
        geometry_df = self.data[[self.unique_id_field, "geometry"]]

        for result_name, result_df in self.results.items():
            # Join geometry on demand
            full_df = geometry_df.merge(result_df, on=self.unique_id_field, how="left")

            # Clean column names (remove _target if possible)
            rename_map = {
                col: col.replace("_target", "")
                for col in full_df.columns
                if col.endswith("_target")
                and col.replace("_target", "") not in full_df.columns
            }
            if rename_map:
                full_df.rename(columns=rename_map, inplace=True)

            for ext in file_types:
                save_path = output_folder / f"{self.name}_{result_name}{ext}"

                if ext == ".csv":
                    full_df.drop(columns="geometry", errors="ignore").to_csv(
                        save_path, index=False
                    )
                elif ext == ".gpkg":
                    full_df.to_file(save_path, layer=self.name, driver="GPKG")
                elif ext == ".shp":
                    full_df.to_file(save_path, driver="ESRI Shapefile")

                print(f"Saved {ext.upper()}: {save_path}")

    def _validate_metric_data(
        self,
        df: pd.DataFrame,
        unique_id: str = None,
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
            set_index: If True, set index to self.unique_id_field.

        Returns:
            Cleaned DataFrame or GeoDataFrame with validated numeric metrics.
        """
        df = df.copy()
        is_geo = isinstance(df, GeoDataFrame)

        # Normalize columns
        df.columns = df.columns.map(str)
        unique_id = unique_id or self.unique_id_field
        reserved = {
            unique_id,
            self.weighting_field,
            self.type_field,
            *(not_allowed or []),
        }

        if metric_field_names is None:
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
        required_cols = [unique_id] + metric_field_names
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
            print(f"\u26a0 Dropping {count} invalid rows.")
            df = df[~invalid_rows]

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
            suffix: Suffix to append to new standardised columns (e.g., "z" → "2020z").

        Returns:
            DataFrame (or GeoDataFrame) with standardised fields.
        """

        protected = {self.unique_id_field, self.weighting_field, self.type_field}

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
            suffix: Suffix for normalised output columns (e.g., "n" → "2020n").

        Returns:
            DataFrame (or GeoDataFrame) with normalised fields.
        """

        protected = {self.unique_id_field, self.weighting_field, self.type_field}

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
                print(f"\u26a0 Skipping {col} — no variation (min = max = {min_val})")
                continue
            self.data[f"{col}{suffix}"] = (self.data[col] - min_val) / (
                max_val - min_val
            )
        print(f"Normalised fields are: {[f'{col}{suffix}' for col in fields]}")

        return self.data

    def plot(self) -> None:
        """
        Plot the base scale spatial data.

        Args:
            base_scale: SpatialScale object containing the base data
            title: Optional title for the plot
        """
        ax = self.data.plot(color="lightblue", edgecolor="blue")
        ax.set_title(f"Base Scale: {self.name}")
        ax.axis("off")
        plt.show()


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
