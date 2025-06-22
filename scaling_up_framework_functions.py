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

import os
import pandas as pd
import numpy as np
import geopandas as gpd
from typing import List, Optional, Union
from pathlib import Path

# =============================================================================
# IMPORT SETTINGS FROM CONFIG FILE
# =============================================================================
# Try to load settings from config.py, use defaults if file not found
try:
    from config import MILLENNIUM_DROUGHT_YEARS, DEFAULT_CRS
except ImportError:
    # Fallback values if config not available (shouldn't normally happen)
    print("Warning: Could not load config.py, using default values")
    MILLENNIUM_DROUGHT_YEARS = [2002, 2003, 2006, 2007, 2008, 2009]
    DEFAULT_CRS = "EPSG:3577"


def validate_csv_data(df: pd.DataFrame, metric_name: str, unique_id: str) -> pd.DataFrame:
    """
    Validate that CSV data contains required columns and has proper format.
    
    Args:
        df: Input dataframe to validate
        metric_name: Name of the metric column that should exist
        unique_id: Name of the unique ID column that should exist
        
    Returns:
        Validated dataframe
        
    Raises:
        ValueError: If required columns are missing or data format is invalid
    """
    # Check for required columns
    required_cols = [unique_id, metric_name]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Required columns {missing_cols} not found in CSV.\n"
                        f"Available columns: {list(df.columns)}\n"
                        f"Please ensure your CSV contains '{unique_id}' and '{metric_name}' columns.")

    # Check for year column (optional)
    has_year = "year" in df.columns
    if has_year:
        print(f"Year column detected - will perform temporal analysis")
        # Validate year column
        if df["year"].isna().any():
            print("Warning: Some year values are missing - these rows will be excluded")
            df = df.dropna(subset=["year"])
        # Convert year to integer
        df["year"] = df["year"].astype(int)
    else:
        print(f"No year column found - will perform single-time aggregation")

    # Check for missing values in required columns
    for col in required_cols:
        missing_count = df[col].isna().sum()
        if missing_count > 0:
            print(f"Warning: {missing_count} missing values in '{col}' column - these rows will be excluded")
            df = df.dropna(subset=[col])

    # Validate metric column is numeric
    if not pd.api.types.is_numeric_dtype(df[metric_name]):
        try:
            df[metric_name] = pd.to_numeric(df[metric_name], errors='coerce')
            nan_count = df[metric_name].isna().sum()
            if nan_count > 0:
                print(f"Warning: {nan_count} non-numeric values in '{metric_name}' converted to NaN")
                df = df.dropna(subset=[metric_name])
        except:
            raise ValueError(f"Column '{metric_name}' contains non-numeric data that cannot be converted")

    print(f"Validation complete: {len(df)} valid records")
    return df


def pivot_year(df: pd.DataFrame, metric_name: str, unique_id: str = "UID") -> pd.DataFrame:
    """
    Pivot ecological time series data to years as columns and calculate baseline statistics.
    Only used when year column is present in the data.
    
    Args:
        df: Input dataframe with year, metric columns
        metric_name: Name of metric column to pivot
        unique_id: Primary key column name
        
    Returns:
        DataFrame with pivoted years and calculated baseline statistics
    """
    if "year" not in df.columns:
        # No year column - return data as-is
        return df.set_index(unique_id)  #[[metric_name]]
    
    # Get all available years
    years = df["year"].unique().tolist()
   
    pivoted_data = df.pivot(index=unique_id, columns="year", values=metric_name)
    
    # Calculate statistics for all years
    years_data = pivoted_data[years]
    
    baseline_stats = {
        "count": years_data.count(axis=1, numeric_only=True),
        "baseline": years_data.mean(axis=1, numeric_only=True),
        "stddev": years_data.std(axis=1, numeric_only=True),
        "max": years_data.max(axis=1, numeric_only=True),
        "median": years_data.median(axis=1, numeric_only=True)
    }
    
    # Calculate MAD (Median Absolute Deviation) - robust measure of variability
    # MAD is preferred over standard deviation for ecological data with outliers
    # Common in species counts, vegetation indices with extreme values
    mad_values = []
    for year in years:
        mad_values.append(abs(baseline_stats["median"] - pivoted_data[year]))
        # set column
    baseline_stats["mad"] = pd.concat(mad_values, axis=1).median(axis=1, numeric_only=True)
    
    print (baseline_stats.keys())
    
    #set series headings for each baseline_stats
    for key, value in baseline_stats.items():
        baseline_stats[key] = value.rename(f"{metric_name}_{key}")
    
    # Combine all results
    result_dfs = [pivoted_data] + list(baseline_stats.values())
    return pd.concat(result_dfs, axis=1)


def aggregate_measure_weighted(
    df: pd.DataFrame,
    aggregator: pd.DataFrame,
    boundary_fields: List[str],
    measure_field: str,
    functional_groups: Optional[List[str]] = None,
    metric_columns: List[str] = None,
    aggregation_method: str = "geometry",
    is_base_scale: bool = False,
) -> pd.DataFrame:
    """
    Aggregate ecological data across spatial scales using appropriate weighting methods.

    This is the core aggregation function that scales up ecological metrics from
    fine-scale units (sites, patches, segments) to larger management units
    (wetland complexes, valleys, basins) using ecologically appropriate weighting.

    Ecological Examples:
        - Vegetation NDVI: Area-weighted mean across habitat patches → valley condition
        - Stream health: Length-weighted mean across river segments → catchment health
        - Species counts: Sum across monitoring sites → regional abundance
        - Habitat patches: Count of patches → landscape fragmentation metric

    Weighting Logic:
        - Points: Equal weighting (count=1) or sum of attribute values
        - all:  frequency of occurrence within the aggregation boundary
        - Lines: Length-weighted for metrics affected by stream/corridor extent
        - Polygons: Area-weighted for metrics affected by habitat/patch size

    Args:
        df: Parameter values per base scale unit (sites, patches, segments)
        aggregator: Aggregator spatial units (complexes, valleys) or None for base scale
        boundary_fields: Fields to group by for aggregation (e.g., ["ValleyName"])
        measure_field: Field containing area/length/count for weighting
        functional_groups: Additional grouping fields (e.g., ["HabitatType"])
        is_base_scale: True if this is the base scale (no spatial aggregation)

    Returns:
        Measure-weighted aggregated values with ecological meaning preserved
    """
    if functional_groups is None:
        functional_groups = []
    print(f"Functional groups: {functional_groups}")

    if metric_columns is None:
        metric_columns = (
            df.select_dtypes(include=[np.number])
            .columns.difference([measure_field])
            .tolist()
        )
    print(f"metric_columns : {metric_columns}")

    grouping_columns = boundary_fields + functional_groups
    print(f"grouping_columns : {grouping_columns}")

    # Ensure aggregator is valid before accessing columns
    if aggregation_method not in ["geometry", "count", "sum", "frequency"]:
        raise ValueError(
            f"Invalid aggregation method: {aggregation_method}\n"
            f"Supported methods: 'geometry', 'count', 'sum', 'frequency'"
        )

    if aggregation_method == "count":
        # Count aggregation - a simple count of features within the aggregation boundary
        # This counts how many base units contribute to each aggregation boundary
        # Group by aggregation boundaries and count base units
        return aggregator.groupby(grouping_columns).size().to_frame(name="count")

    joined_data = df.join(aggregator, how="inner").replace([np.inf, -np.inf], np.nan)
    print(f"Aggregating {len(joined_data)} records across {len(aggregator)} boundaries")

    # -----------------------------------------------------------------
    # Sum aggregation
    # -----------------------------------------------------------------

    # Sum aggregation - sums attribute values across aggregation boundaries
    # This is used when the attribute values themselves are meaningful
    # E.g., total species abundance, carbon storage, management costs
    # It does not consider the size of the features, just their values

    if aggregation_method == "sum":
        # Group by aggregation boundaries and sum base unit values
        # This sums the base unit values for each aggregation boundary
        aggregated_results = (
            joined_data[metric_columns + grouping_columns]
            .groupby(grouping_columns)[metric_columns]
            .sum(numeric_only=True)
        )

    elif aggregation_method == "frequency":

        # -----------------------------------------------------------------
        # Frequency Aggregation
        # -----------------------------------------------------------------

        if not functional_groups:
            raise ValueError(
                "Frequency weighting requires functional_groups to define subgroups within boundaries."
            )

        # Step 1: Frequency counts per boundary–group
        frequency_counts = joined_data.groupby(grouping_columns).size().rename("freq")

        # Step 2: Total counts per boundary
        total_counts = frequency_counts.groupby(boundary_fields).sum().rename("total")

        # Step 3: Join to get frequency weights
        frequency_weights = frequency_counts.to_frame().join(
            total_counts, on=boundary_fields
        )
        frequency_weights["weight"] = (
            frequency_weights["freq"] / frequency_weights["total"]
        )

        # Step 4: Join weights back to the data
        joined_data = joined_data.join(frequency_weights["weight"], on=grouping_columns)

        # Step 5: Apply weights to metrics
        weighted_metrics = joined_data[metric_columns].multiply(
            joined_data["weight"], axis=0
        )

        # Step 6: Aggregate weighted metrics to boundary
        aggregated_results = weighted_metrics.groupby(boundary_fields).sum()

    else:
        # -----------------------------------------------------------------
        # Geometry Aggregation
        # -----------------------------------------------------------------

        # Geometry aggregation - gives more weight to larger features and requires a measure field
        # to quantity the area/length of each base unit contributing to the aggregation

        # Weighting is achieved by multiplying metrics by area/length/count
        # This ensures larger features contribute more to landscape metrics
        # than smaller patches, preserving ecological meaning
        if measure_field not in joined_data.columns:
            raise ValueError(
                f"Measure field '{measure_field}' not found in joined data.\n"
                f"Available columns: {list(joined_data.columns)}\n"
                f"Please ensure the measure field is correctly specified."
            )
        joined_data[metric_columns] = joined_data[metric_columns].multiply(
            joined_data[measure_field], axis=0
        )
        # print ( joined_data[metric_columns])
        # Group by aggregation boundaries and sum weighted values
        # Sums both the weighted metrics and the measures (total area/length/count)
        aggregated_results = (
            joined_data[metric_columns + [measure_field] + grouping_columns]
            .groupby(grouping_columns)
            .sum(numeric_only=True)
        )
        # Normalize by total measure to get proper weighted averages
        # Converts weighted sums back to meaningful ecological units
        # E.g., (sum of area*NDVI) / (total area) = area-weighted mean NDVI
        aggregated_results[metric_columns] = aggregated_results[metric_columns].div(
            aggregated_results[measure_field], axis=0
        )
    if aggregation_method in ["sum", "geometry", "frequency"]:
        return aggregated_results[
            metric_columns
            + ([measure_field] if measure_field in aggregated_results.columns else [])
        ]
    elif aggregation_method == "count":
        return aggregated_results


def aggregate_count(
    df: pd.DataFrame,
    aggregator: pd.DataFrame,
    boundary_fields: List[str],
    functional_groups: Optional[List[str]] = None,
    is_base_scale: bool = False,
) -> pd.DataFrame:
    """
    Aggregate ecological data across spatial scales using a simple count.

    This is the core aggregation function that scales up ecological metrics from
    fine-scale units (sites, patches, segments) to larger management units
    (wetland complexes, valleys, basins) using a simple count of base units.

    Ecological Examples:
        - Habitat patches: Count of patches → landscape fragmentation metric
        - Monitoring sites: Count of sites → regional abundance

    Args:
        df: Parameter values per base scale unit (sites, patches, segments)
        aggregator: Aggregator spatial units (complexes, valleys) or None for base scale
        boundary_fields: Fields to group by for aggregation (e.g., ["ValleyName"])
        functional_groups: Additional grouping fields (e.g., ["HabitatType"])
        is_base_scale: True if this is the base scale (no spatial aggregation)

    Returns:
        Count of base units aggregated across aggregation boundaries
    """
    if functional_groups is None:
        functional_groups = []

    # Ensure df is valid before accessing columns
    if df is None or df.empty:
        raise ValueError("Input DataFrame is None or empty")

    metric_columns = df.select_dtypes(include=[np.number]).columns.values.tolist()

    grouping_columns = boundary_fields + functional_groups

    # # Base scale - no spatial aggregation, just group by functional types
    # # Used for base-level analysis (e.g., NDVI by habitat type within sites)
    if is_base_scale:
        # grp column already exists from base scale loading (GROUP_RULES applied)
        # numeric_cols = joined_data.select_dtypes(include=[np.number]).columns.tolist()
        # return joined_data[metric_columns + [measure_field] + grouping_columns].set_index(grouping_columns, append=True)
        return joined_data

    print(f"Aggregating {len(df)} records across {len(aggregator)} boundaries")

    joined_data = df.join(aggregator, how="inner").replace([np.inf, -np.inf], np.nan)
    # Group by aggregation boundaries and sum weighted values
    # Sums both the weighted metrics and the measures (total area/length/count)
    aggregated_results = (
        joined_data[metric_columns + grouping_columns]
        .groupby(grouping_columns)
        .count(numeric_only=True)
    )
    # print(aggregated_results)

    return aggregated_results[metric_columns]


def standardise_z(df: pd.DataFrame) -> pd.DataFrame:
    """
    Z-score standardization of input dataframe.
    
    Args:
        df: Input dataframe to standardize
        
    Returns:
        Z-score standardized dataframe
    """
    return (df - df.mean()) / df.std()


def normalise(df: pd.DataFrame) -> pd.DataFrame:
    """
    Min-max normalization to scale values between 0 and 1.
    
    Args:
        df: Input dataframe to normalize
        
    Returns:
        Normalized dataframe with values between 0-1
    """
    df_min = df.min().min() if hasattr(df.min(), 'min') else df.min()
    df_max = df.max().max() if hasattr(df.max(), 'max') else df.max()
    
    return (df - df_min) / (df_max - df_min)


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
        measure_field: Optional[str] = None,
        measure_multiplier: Optional[float] = None,
        type_field: Optional[str] = None,
        aggregation_method: str = "geometry",
        extra_columns: Optional[List[str]] = None,
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
            measure_field: Column containing size/importance measure (or None to auto-calculate)
            measure_multiplier: Number to multiply measures by (usually leave as None)
            aggregation_method: How to weight features ("geometry", "count", or "sum")
            type_field: Column containing feature types for grouping (or None for no grouping)
            extra_columns: Additional columns to load from shapefile
            is_base_scale: True if this is your smallest scale, False for larger scales
        """
        self.name = name
        self.source = Path(source)
        self.is_base_scale = is_base_scale

        # Check that the shapefile actually exists
        if not self.source.exists():
            raise FileNotFoundError(f"Cannot find shapefile: {source}\n"
                                  f"Please check that the file exists and the path is correct.")

        # Convert unique_id_field to a list format (technical requirement)
        # This allows for compound keys if needed (multiple columns)
        self.unique_id = [unique_id_field] if isinstance(unique_id_field, str) else unique_id_field

        # Build list of columns to load from the shapefile
        # Only load what we need to make processing faster
        cols = self.unique_id.copy()  # Always need the ID field
        if measure_field:
            cols.append(measure_field)    # Size/importance measure
        if type_field:
            cols.append(type_field)       # Feature types for grouping
        if extra_columns:
            cols.extend([c for c in extra_columns if c not in cols])  # Any additional fields

        # Load the shapefile and convert to standard coordinate system
        # This ensures all spatial calculations are accurate
        print(f"Loading {name} data from {self.source.name}...")
        self.data = gpd.read_file(self.source, columns=cols).to_crs(DEFAULT_CRS)

        # Handle complex geometries by breaking them into simple parts
        # E.g., a MultiPolygon becomes multiple separate Polygons
        # This ensures each feature is treated as a separate unit for analysis
        multi_geom_count = self.data.geometry.geom_type.isin(['MultiPolygon', 'MultiLineString', 'MultiPoint']).sum()
        if multi_geom_count > 0:
            print(f"  Breaking {multi_geom_count} complex geometries into simple parts...")
            self.data = self.data.explode(index_parts=False).reset_index(drop=True)

        self.data = validate_geometries(self.data)

        # Check that all required columns exist in the shapefile
        missing_cols = [c for c in cols if c not in self.data.columns]
        if missing_cols:
            raise ValueError(f"Required columns {missing_cols} not found in {source}\n"
                           f"Available columns: {list(self.data.columns)}\n"
                           f"Please check your column names in config.py")

        # Automatically detect what type of geometry this shapefile contains
        # This determines how we calculate measures and perform aggregation
        first_geom = self.data.geometry.iloc[0]
        if first_geom.geom_type == 'Polygon':
            geometry_type = "polygon"    # Areas (wetlands, habitats, management units)
        elif first_geom.geom_type == 'LineString':
            geometry_type = "line"       # Linear features (rivers, transects, corridors)
        elif first_geom.geom_type == 'Point':
            geometry_type = "point"      # Point locations (monitoring sites, observations)
            if aggregation_method not in ["count", "sum"]:
                print(
                    f"  WARNING: Point data detected. Setting aggregation method to 'count'."
                )
                aggregation_method = "count"
        else:
            raise ValueError(f"Unsupported geometry type: {first_geom.geom_type}")

        print(f"  Detected geometry type: {geometry_type} ({len(self.data)} features)")

        # Calculate appropriate measure based on aggregation method and geometry
        # This determines how ecological features are weighted in aggregation
        if not measure_field:
            if aggregation_method == "count":
                # Count aggregation: Each feature contributes equally (presence/absence)
                # Use for: Number of habitat patches, monitoring sites, species occurrences
                self.data["Count"] = 1
                measure_field = "Count"
            elif aggregation_method == "sum":
                # Sum aggregation: Add attribute values, ignoring feature size
                # Use for: Total species abundance, carbon storage, management costs
                # Default to 1, but user should provide actual values in measure_field
                self.data["Sum"] = 1  
                measure_field = "Sum"
            elif geometry_type == "polygon":
                # Area-weighted aggregation for habitat/vegetation patches
                # Larger patches contribute more to landscape-level metrics
                # Use for: Vegetation condition, habitat quality, land cover metrics
                self.data["Area_Ha"] = self.data.geometry.area / 10000  # Convert m² to hectares
                measure_field = "Area_Ha"
            elif geometry_type == "line":
                # Length-weighted aggregation for linear features
                # Longer segments contribute more to network-level metrics
                # Use for: Stream health, corridor connectivity, infrastructure
                self.data["Length_m"] = self.data.geometry.length
                measure_field = "Length_m"
            elif geometry_type == "point":
                # Point features default to count (each location = 1 unit)
                # Use for: Monitoring site density, species occurrence frequency
                self.data["Count"] = 1
                measure_field = "Count"

        # Apply measure multiplier if specified (e.g., convert units)
        if measure_multiplier and measure_field:
            print(f"  Applying multiplier {measure_multiplier} to {measure_field}")
            self.data[measure_field] *= measure_multiplier

        self.measure_field = measure_field
        self.type_field = type_field
        self.geometry_type = geometry_type

    def __str__(self) -> str:
        """Create a human-readable description of this spatial scale."""
        scale_type = "BASE SCALE" if self.is_base_scale else "AGGREGATION SCALE"
        return f"{self.name} ({scale_type}) - {len(self.data)} features from {self.source.name}"


def validate_geometries(gdf: gpd.GeoDataFrame, fix_invalid: bool = True) -> gpd.GeoDataFrame:
    """
    Validate and optionally fix invalid geometries.
    
    Args:
        gdf: GeoDataFrame to validate
        fix_invalid: Whether to attempt fixing invalid geometries
        
    Returns:
        GeoDataFrame with validated geometries
    """
    from shapely.validation import make_valid, explain_validity
    
    invalid_mask = ~gdf.geometry.is_valid
    invalid_count = invalid_mask.sum()
    
    if invalid_count == 0:
        print("All geometries are valid.")
        return gdf
    
    print(f"Found {invalid_count} invalid geometries")
    
    if fix_invalid:
        print("Attempting to fix invalid geometries...")
        gdf.loc[invalid_mask, 'geometry'] = gdf.loc[invalid_mask, 'geometry'].apply(make_valid)
        
        # Check if fixes worked
        still_invalid = (~gdf.geometry.is_valid).sum()
        if still_invalid == 0:
            print("All geometries successfully fixed.")
        else:
            print(f"Warning: {still_invalid} geometries remain invalid after repair.")
    
    return gdf


def plot_base_scale(base_scale, base_scale_config: dict) -> None:
    """
    Plot the base scale spatial data.
    
    Args:
        base_scale: SpatialScale object containing the base data
        base_scale_config: BASE_SCALE configuration dictionary
    """
    import matplotlib.pyplot as plt
    
    base_scale.data.plot(figsize=(10, 8), color='lightblue', edgecolor='blue')
    plt.title(f'Base Scale: {base_scale_config["name"]}')
    plt.show()


def plot_spatial_hierarchy(base_scale, base_scale_config: dict, agg_scale_objects: dict) -> None:
    """
    Plot all spatial scales together with proper legend and extent.
    
    Args:
        base_scale: SpatialScale object containing the base data
        base_scale_config: BASE_SCALE configuration dictionary
        agg_scale_objects: Dictionary of aggregation scale objects
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import matplotlib.cm as cm
    import numpy as np
    
    fig, ax = plt.subplots(figsize=(12, 10))
    base_scale.data.plot(ax=ax, color='lightblue', alpha=0.6, edgecolor='blue', label=base_scale_config["name"])
    
    # Generate colors based on number of aggregation scales
    n_scales = len(agg_scale_objects)
    colors = cm.Set1(np.linspace(0, 1, n_scales)) if n_scales > 0 else []
    
    for i, (name, scale_obj) in enumerate(agg_scale_objects.items()):
        scale_obj.data.plot(ax=ax, color='none', edgecolor=colors[i], linewidth=2, label=name)
    
    # Set axis limits to base scale extent
    ax.set_xlim(base_scale.data.total_bounds[0], base_scale.data.total_bounds[2])
    ax.set_ylim(base_scale.data.total_bounds[1], base_scale.data.total_bounds[3])
    
    # Create custom legend
    legend_elements = [Line2D([0], [0], color='lightblue', lw=4, label=base_scale_config["name"])]
    for i, name in enumerate(agg_scale_objects.keys()):
        legend_elements.append(Line2D([0], [0], color=colors[i], lw=2, label=name))
    ax.legend(handles=legend_elements)
    
    ax.set_title('Spatial Scale Hierarchy')
    plt.show()


def load_metric_from_spatial(base_scale, metric_name: str, unique_id: str) -> pd.DataFrame:
    """
    Load metric data directly from spatial data attributes.
    
    Args:
        base_scale: SpatialScale object containing the spatial data
        metric_name: Name of the attribute column to use as metric
        unique_id: Name of the unique ID column
        
    Returns:
        DataFrame with unique_id and metric columns
    """
    print(f"Loading '{metric_name}' from spatial data...")
    
    # Check if metric column exists
    if metric_name not in base_scale.data.columns:
        raise ValueError(f"Column '{metric_name}' not found in spatial data.\n"
                        f"Available columns: {list(base_scale.data.columns)}")
    
    # Create DataFrame from spatial data
    data = base_scale.data[[unique_id, metric_name]].copy()
    data = data.reset_index(drop=True)
    
    # Validate metric column is numeric
    if not pd.api.types.is_numeric_dtype(data[metric_name]):
        data[metric_name] = pd.to_numeric(data[metric_name], errors='coerce')
        nan_count = data[metric_name].isna().sum()
        if nan_count > 0:
            print(f"Warning: {nan_count} non-numeric values converted to NaN")
            data = data.dropna(subset=[metric_name])
    
    print(f"Loaded {len(data)} records from spatial data")
    return data


def assign_functional_group(type_value: str, base_scale_name: str, group_rules: dict, include_unmatched: bool = True) -> Optional[str]:
    """
    Assign functional group based on base scale type field using configurable rules.
    
    Args:
        type_value: Type classification string from base scale
        base_scale_name: Name of the base scale to get appropriate rules
        group_rules: Dictionary of grouping rules by base scale
        include_unmatched: If True, return original type_value for unmatched types; if False, return None
        
    Returns:
        Assigned functional group name, original type_value, or None based on configuration
    """
    if base_scale_name not in group_rules:
        return None

    type_lower = type_value.lower()
    scale_rules = group_rules[base_scale_name]

    # Check each group - convert keywords to lowercase for case-insensitive matching
    for group_name, keywords in scale_rules.items():
        keywords_lower = [keyword.lower() for keyword in keywords]
        if any(keyword in type_lower for keyword in keywords_lower):
            return group_name

    # Handle unmatched types based on configuration
    return type_value if include_unmatched else None
