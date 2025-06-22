# MDB Scaling-up Framework

A generic spatial aggregation framework for ecological data analysis across multiple spatial scales.

## Overview

This framework aggregates ecological data from fine-scale base units (sites) to progressively larger spatial scales (wetland complexes, valleys, river basins) using appropriate weighting methods based on geometry type and research objectives.

## Aggregation Methods in Wetland Ecology

### Area-Weighted Aggregation
**Why it's useful:** Larger wetland patches contribute more to landscape-level metrics because they support more species, store more carbon, and provide greater ecosystem services.

**Wetland Example:** Calculating average vegetation health (NDVI) across a valley
- Small marsh (2 ha, NDVI=0.8) + Large swamp (50 ha, NDVI=0.6)
- Area-weighted result: (2×0.8 + 50×0.6) ÷ 52 = 0.61
- Simple average would be: (0.8 + 0.6) ÷ 2 = 0.70 (misleading)

**Configuration:**
```python
BASE_SCALE = {
    "aggregation_method": "geometry",  # Uses area weighting for polygons
    "measure_field": "Area_Ha"  # Optional - auto-calculated if not provided
}
```

### Length-Weighted Aggregation
**Why it's useful:** Longer stream segments have greater influence on watershed health, connectivity, and water quality than shorter segments.

**Wetland Example:** Assessing stream health across a catchment
- Short pristine creek (2 km, quality=0.9) + Long degraded river (20 km, quality=0.4)
- Length-weighted result: (2×0.9 + 20×0.4) ÷ 22 = 0.45
- Simple average would be: (0.9 + 0.4) ÷ 2 = 0.65 (overestimates health)

**Configuration:**
```python
BASE_SCALE = {
    "aggregation_method": "geometry",  # Uses length weighting for lines
    "measure_field": "Length_km"  # Optional - auto-calculated if not provided
}
```

### Count-Based Aggregation
**Why it's useful:** For presence/absence data, species richness, or when each feature represents equal ecological importance regardless of size.

**Wetland Example:** Counting rare species occurrences across wetland complexes
- Each monitoring site contributes equally regardless of wetland size
- Large wetland (5 sites, 3 species each) = 15 total species
- Small wetland (2 sites, 4 species each) = 8 total species

**Configuration:**
```python
BASE_SCALE = {
    "aggregation_method": "count",  # Each feature counts as 1
    # OR
    "aggregation_method": "sum",   # Sum attribute values
    "measure_field": "SpeciesCount"
}
```

## Key Features

- **Auto-detection** of geometry types (Point, Line, Polygon)
- **Multi-geometry support** with automatic exploding
- **Ecologically appropriate weighting** geometry (area or length)], count and sum
- **Optional grouping** with configurable rules or original type field fallback
- **Temporal analysis** with baseline statistics and drought exclusion
- **Quality validation** with geometry repair capabilities
- **Multi-scale aggregation** across spatial hierarchies

## Ecological Applications

This framework supports diverse ecological research questions:

- **Biodiversity assessment** across spatial scales
- **Habitat connectivity** analysis  
- **Water quality** monitoring aggregation
- **Species distribution** modeling
- **Conservation planning** at multiple scales
- **Ecosystem service** quantification
- **Climate change** impact assessment

## Spatial Scale Hierarchy

The code is currently set up for the hierarchy of scales pictured below but can be changed to use additional or fewer intermediate scales easily.  A generic SpatialUnit Class is used to define the scales using any compatible set of polygons as a shape file or ESRI ArcGIS polygon Feature Class.  Here we use Basin Valleys, and two different Wetland polygon data sets (DIWA and Ramsar).  THe Base Scale are the building blocks that will be aggregated across the different spatial scale.  The code here uses the ANAE wetland and floodplain polygon data set, and also the Geofabric river lines.  The Code framework accepts and autodetects the Base Scale data as points, lines or polygons.

```[]
River Basin (Largest)
├── Valleys
│   ├── Wetland Complexes  
│   │   ├── Individual Sites (Base Scale)
│   │   │   ├── Points (monitoring stations, species observations)
│   │   │   ├── Lines (river segments, transects)
│   │   │   └── Polygons (vegetation patches, habitat areas)
```

## Geometry Types & Aggregation Methods

### Point Data Examples

**Ecological Use Cases:**

- Bird observation sites
- Water quality monitoring stations  
- Species occurrence records
- Infrastructure locations

**Aggregation Methods:**

- **Count**: Number of monitoring sites per valley
- **Sum**: Total bird abundance across all sites in a wetland complex
- **Geometry (Count)**: Same as count for points

### Line Data Examples  

**Ecological Use Cases:**

- River segments
- Transect lines
- Migration corridors
- Fence lines/barriers

**Aggregation Methods:**

- **Count**: Number of river segments per valley (ignoring length)
- **Sum**: Total fish species richness across segments (ignoring length)
- **Geometry (Length)**: Total river length per wetland complex

### Polygon Data Examples

**Ecological Use Cases:**

- Vegetation patches
- Habitat areas
- Management units
- Disturbance zones

**Aggregation Methods:**

- **Count**: Number of vegetation patches per valley (ignoring size)
- **Sum**: Total species count across patches (ignoring patch size)
- **Geometry (Area)**: Total habitat area per wetland complex

## Practical Examples

### Example 1: Bird Monitoring Program

```[]
Base Scale: Point observations at 50 sites
Aggregation: Wetland complexes (5 complexes)
Method: SUM - Total bird abundance per complex
Result: Complex A = 1,250 birds, Complex B = 890 birds
```

### Example 2: River Health Assessment  

```[]
Base Scale: 200 river segments (lines)
Aggregation: Valleys (8 valleys)
Method: LENGTH - Water quality weighted by segment length
Result: Valley 1 = 85km of good quality river, Valley 2 = 45km
```

### Example 3: Vegetation Condition

```[]
Base Scale: 1,000 habitat polygons
Aggregation: River basin (1 basin)
Method: AREA - NDVI weighted by polygon area
Result: Basin-wide vegetation condition = 0.67 (normalized NDVI)
```

## Configuration Examples

### Point Data (Species Observations)

```python
BASE_SCALE = {
    "name": "ObservationSites",
    "file": "species_points.shp",
    "unique_id": "SiteID", 
    "measure_field": "SpeciesCount",
    "type_field": "HabitatType",  # Optional - for grouping
    "aggregation_method": "sum"  # Sum species across sites
}
```

### Line Data (River Segments)

```python
BASE_SCALE = {
    "name": "RiverSegments", 
    "file": "river_lines.shp",
    "unique_id": "SegmentID",
    "measure_field": "Length_km",
    "type_field": "StreamOrder",  # Optional - for grouping
    "aggregation_method": "geometry"  # Length-weighted aggregation
}
```

### Polygon Data (Vegetation Patches)

```python
BASE_SCALE = {
    "name": "VegetationPatches",
    "file": "habitat_polygons.shp", 
    "unique_id": "PatchID",
    "measure_field": "Area_Ha",
    "type_field": "VegetationType",  # Optional - for grouping
    "aggregation_method": "geometry"  # Area-weighted aggregation
}
```

## When to Use Each Aggregation Method

### Count Aggregation

**Use when:** Feature presence/density matters more than size

- **Points**: Number of monitoring sites per region
- **Lines**: Number of stream segments per catchment  
- **Polygons**: Number of habitat patches per landscape

### Sum Aggregation  

**Use when:** Attribute totals matter, ignoring feature size

- **Points**: Total species abundance across sites
- **Lines**: Total species richness across segments
- **Polygons**: Total carbon storage across patches

### Geometry Aggregation

**Use when:** Feature size/extent affects the metric

- **Points**: Automatically becomes count
- **Lines**: Stream length affects water volume calculations
- **Polygons**: Habitat area affects population capacity

## Output Files

Results are saved as CSV files following the pattern:

```[]
{metric}_{spatial_scale}.csv
```

Examples:

- `NDVI_Sites.csv` - Base scale vegetation condition
- `NDVI_WetlandComplexes.csv` - Complex-level aggregation  
- `NDVI_Valleys.csv` - Valley-level aggregation
- `NDVI_Basin.csv` - Basin-wide aggregation

## CSV Data Requirements

### Required Columns
Your CSV files must contain:
- **unique_id column**: Must match `BASE_SCALE["unique_id"]` in config.py
- **metric column**: Numeric values to aggregate (e.g., NDVI, SpeciesCount)

### Optional Columns
- **year column**: Named exactly "year" (case-sensitive)

### Behavior Differences

**With Year Column (Temporal Analysis):**
```csv
UID,NDVI,year
001,0.65,2020
001,0.72,2021
002,0.58,2020
```
- Calculates baseline statistics excluding drought years
- Provides temporal trends and variability measures
- Outputs include: baseline mean, std dev, max, median, MAD

**Without Year Column (Single Time Point):**
```csv
UID,SpeciesCount
001,15
002,23
003,8
```
- Simple spatial aggregation across scales
- No temporal statistics calculated
- Faster processing for snapshot analyses

## Installation & Usage

1. Install requirements:

```bash
pip install -r requirements.txt
```

2. Configure your data in `config.py`

3. Prepare your CSV data (see requirements above)

4. Run the analysis:

```python
jupyter notebook scaling_up_framework.ipynb
```

## Grouping Behavior

### type_field Configuration

- **type_field = None**: All features aggregated together (single output per spatial scale)
- **type_field defined + GROUP_RULES**: Uses functional groups for aggregation
- **type_field defined, no GROUP_RULES**: Uses original type field values for grouping

### GROUP_RULES (Optional)

Define functional groupings to combine similar feature types:

```python
GROUP_RULES = {
    "VegetationPatches": {
        "wetland_vegetation": ["marsh", "swamp", "bog"],
        "riparian_forest": ["river red gum", "black box"],
        "grassland": ["meadow", "pasture"]
    }
}
```
