"""
Configuration Settings for MDB Scaling-up Framework

This file contains all the settings that control how spatial data is processed
and aggregated across different scales.

What this framework does:
- Takes small spatial features (like individual wetland patches or monitoring sites)
- Groups them into larger management units (wetland complexes,  valleys or river basins)
- Calculates weighted averages that account for the contribution of each feature either 
  by size (area or length) or by magnitude of an attribute (e.g. an abundance or diversity measure)

Key concepts:
- BASE_SCALE: Your smallest spatial units (the building blocks)
- AGGREGATION_SCALES: Larger boundaries that contain multiple base units
- WEIGHTING: How much each small unit contributes to the larger average
"""

from pathlib import Path

# =============================================================================
# FILE PATHS
# =============================================================================
# Tell the program where to find the data files
# Use relative paths (starting with ./) so the code works on any computer
SPATIAL_PATH = Path("./input/spatial/")  # Folder containing shapefiles
DATA_PATH = Path("./input/csv/")         # Folder containing your metric data (CSV files)
OUTPUT_PATH = Path("./output/")          # Folder where results will be saved

# =============================================================================
# METRIC DATA SOURCES
# =============================================================================
# The framework can load metric data from two sources:
#
# SOURCE 1: CSV FILES (for temporal or external data)
# Your CSV files must contain these columns:
# - unique_id column: Must match BASE_SCALE["unique_id"] field name
# - metric column: Contains the values to aggregate (e.g., "NDVI", "SpeciesCount")
# - year column: Optional, named exactly "year" (enables temporal analysis)
#
# SOURCE 2: BASE_SCALE SPATIAL DATA ATTRIBUTES (for single time point data)
# If no CSV file is found, the framework will look for the metric as a column
# in the BASE_SCALE shapefile itself.
#
# EXAMPLE: Using ANAE as BASE_SCALE
# - load_metric_data("NDVI") looks for:
#   1. CSV file containing "NDVI" in filename
#   2. If no CSV: column named "NDVI" in ANAE shapefile
#
# EXAMPLE CSV STRUCTURES:
# With temporal data:
# UID,NDVI,year
# 001,0.65,2020
# 001,0.72,2021
# 
# Without temporal data:
# UID,SpeciesCount
# 001,15
# 002,23
#
# EXAMPLE SPATIAL ATTRIBUTE:
# ANAE shapefile with columns: UID, ANAE_TYPE, Area_Ha, BiomassKg
# - load_metric_data("BiomassKg") uses BiomassKg column directly

# =============================================================================
# SPATIAL DATA FILES
# =============================================================================
# Shapefiles containing your geographic boundaries for different spatial scales
# Each file defines boundaries that will contain multiple base scale units
SPATIAL_FILES = {
    "ANAE": SPATIAL_PATH / "ANAEv3_BWS.shp.zip",        # Base scale wetland polygons
    "Ramsar": SPATIAL_PATH / "ramsar_wetlands.shp.zip",   # International wetland sites
    "DIWA": SPATIAL_PATH / "DIWA_complex.shp.zip",       # Directory of Important Wetlands
    "Valley": SPATIAL_PATH / "BWSRegions.shp.zip",        # Basin valley boundaries
    "NorthSouthBasin": SPATIAL_PATH / "Northern_Southern_Basin.shp.zip",          # Northern and Southern Basin boundaries
}

#-----------------------------------------------------------------
# ANAE wetland polygons as the Base or atomic unit
#-----------------------------------------------------------------

# Base scale configuration - the atomic unit for all aggregations
# IMPORTANT: The unique_id field name here must match the column name in your CSV files
BASE_SCALE = {
    "name": "ANAE",
    "file": SPATIAL_PATH / "ANAEv3_BWS.shp.zip",
    "unique_id": "UID",
    "measure_field": "Area_Ha",  # Area_Ha for polygons, Length_m for lines
    "type_field": "ANAE_TYPE",
    "aggregation_method": "geometry"  # "geometry" (area/length/count), "count" (force count), or "sum" (sum values)
    # geometry_type auto-detected from data
}

# Grouping configuration
# INCLUDE_UNMATCHED_TYPES (True/False)
#  True: The value of type_field is used as the GROUP for any featured not found in the GROUP_RULES definition ensuring all inputs 
#  False: Features where the value of type_field are NOT FOUND in any groups are ignored (excluded from analysis)
INCLUDE_UNMATCHED_TYPES = True  # True: use original type_field, False: ignore (group=None)

# Functional group rules for different base scales
# The rules dictionary maps the base scale name to a dictionary of functional group rules.
# Each functional group rule is a dictionary where the keys are the functional group names
# keyword matching of types to groups is always case insensitive
# the structure is group_name: keywords in the values of the configured type_field

GROUP_RULES = {
    "ANAE": {
        "river red gum woodland": ["river red gum", "woodland"],
        "river red gum swamps and forests": ["river red gum"],
        "black box": ["black box"],
        "coolibah": ["coolibah"],
        "lignum": ["lignum"],
        "cooba": ["cooba"],
        "shrubland": ["f2.4: shrubland riparian zone or floodplain"],
        "submerged lake": ["permanent lake", "permanent wetland", "aquatic bed"],
        "tall reed beds": ["tall emergent marsh"],
        "grassy meadows": ["grass", "meadow"],
        "herbfield": ["forb marsh", "temporary wetland", "temporary lake"],
        "clay pan": ["clay"]
    },
    "Rivers": {
        "Permanent Lowland": ["permanent lowland"],
        "Temporary Lowland": ["temporary lowland"],
        "Other temporary": ["temporary transitional, "],
    }
}


#-----------------------------------------------------------------
# Geofabric v3  river lines (Classified by the ANAE data set)
#-----------------------------------------------------------------

# Base scale configuration - the atomic unit for all aggregations
# BASE_SCALE = {
#     "name": "Rivers",
#     "file": SPATIAL_PATH / "ANAEv3_Rivers_24mar2021_CEWfrequency_added_14feb2025.zip",
#     "unique_id": "UID",
#     "measure_field": "Length_km",  # Area_Ha for polygons, Length_m for lines
#     "type_field": "ANAE_TYPE",
#     "aggregation_method": "geometry"  # "geometry" (area/length/count), "count" (force count), or "sum" (sum values)
#     # geometry_type auto-detected from data
# }

# Functional group rules for different base scales





# Analysis parameters
MILLENNIUM_DROUGHT_YEARS = [2002, 2003, 2006, 2007, 2008, 2009]
DEFAULT_CRS = "EPSG:3577"  # Australian Albers
BASELINE_YEARS = list(range(1989, 2023))  # Excluding drought years

# =============================================================================
# UNUSED PARAMETERS (Legacy - can be removed)
# =============================================================================
# These parameters are not currently used by the framework
# They were included for potential future features
# YEAR_WINDOW_WIDTH = 5        # For moving window analysis (not implemented)
# TREND_WINDOW_WIDTH = 2       # For trend analysis (not implemented) 
# MIN_POLYGON_SIZE_HA = 1.0    # For filtering small polygons (not implemented)

# Output file naming
OUTPUT_PATTERN = "{metric}_{aggregator}_{window}yr{tag}.csv"