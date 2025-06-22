"""
Test Configuration for MDB Scaling-up Framework

Simple test setup with rectangular geometries for validation
"""

from pathlib import Path

# =============================================================================
# TEST FILE PATHS
# =============================================================================
SPATIAL_PATH = Path("./test/")
DATA_PATH = Path("./test/")
OUTPUT_PATH = Path("./test/output/")

# =============================================================================
# TEST SPATIAL DATA FILES
# =============================================================================
SPATIAL_FILES = {
    "BaseUnits": SPATIAL_PATH / "base_units.shp",
    "Regions": SPATIAL_PATH / "regions.shp", 
    "Basin": SPATIAL_PATH / "basin.shp"
}

# =============================================================================
# TEST BASE SCALE CONFIGURATION
# =============================================================================
BASE_SCALE = {
    "name": "BaseUnits",
    "file": SPATIAL_PATH / "base_units.shp",
    "unique_id": "UID",
    "measure_field": "Area_Ha",
    "type_field": "HabitatTyp",  # Truncated by shapefile format
    "aggregation_method": "geometry"
}

# =============================================================================
# TEST GROUPING CONFIGURATION
# =============================================================================
INCLUDE_UNMATCHED_TYPES = True

# Test functional group rules
GROUP_RULES = {
    "BaseUnits": {
        "vegetation": ["Forest", "Grassland"],
        "water": ["Wetland"]
    }
}

# =============================================================================
# TEST ANALYSIS PARAMETERS
# =============================================================================
MILLENNIUM_DROUGHT_YEARS = [2002, 2003, 2006, 2007, 2008, 2009]
DEFAULT_CRS = "EPSG:4326"
BASELINE_YEARS = list(range(2020, 2023))

# =============================================================================
# UNUSED PARAMETERS (Legacy)
# =============================================================================
# YEAR_WINDOW_WIDTH = 5
# TREND_WINDOW_WIDTH = 2
# MIN_POLYGON_SIZE_HA = 1.0

# Output file naming
OUTPUT_PATTERN = "{metric}_{aggregator}_{window}yr{tag}.csv"