# Murray-Darling Basin Scaling-up Framework

A framework for hierarchical data aggregation across spatial scales with a specific focus on scaling up the evaluation of rivers, wetlands and floodplains within large river systems like the Murray-Darling Basin.

## Table of Contents

- [Overview](#overview)
- [Scaling Concepts in Landscape Ecology](#scaling-concepts-in-landscape-ecology)
- [Framework Capabilities](#framework-capabilities)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage Examples](#usage-examples)
- [API Reference](#api-reference)
- [Contributing](#contributing)

## Overview

The Scaling-up Framework addresses a fundamental challenge in landscape ecology: **how to transfer ecological information across spatial scales while preserving biological relevance and statistical validity**.

This framework:

- Aggregate fine-resolution ecological data to management-relevant scales
- Maintains ecological meaning through appropriate weighting schemes
- Handles complex spatial hierarchies common in landscape ecology
- Ensure statistical rigor in multi-scale analyses

## The Scale Challenge

In landscape ecology, **scale** refers to both the spatial extent (how large an area) and resolution (how fine the detail) of ecological observations. The challenge arises because:

1. **Ecological processes operate at multiple scales simultaneously**
   - Individual organisms respond to local habitat conditions
   - Populations are influenced by landscape connectivity
   - Communities are shaped by regional climate patterns

2. **Management decisions occur at different scales than data collection**
   - Data: Individual wetlands, vegetation plots, monitoring sites
   - Management: River reaches, catchments, bioregions

3. **Scale-dependent patterns emerge**
   - Rare species may appear absent at fine scales but present at broad scales
   - Habitat fragmentation effects vary with observation scale
   - Water quality patterns differ between local and regional assessments

### Hierarchical Organization

Ecological systems commonly exhibit **hierarchical organization** where smaller units nest within larger ones:

```
Individual Wetlands → River Reaches → Sub-catchments → Catchments → Bioregions
    (Fine scale)                                                    (Broad scale)
```

### Aggregation Challenges

Simply averaging fine-scale data to broader scales can be **ecologically misleading** because:

1. **Unequal representation**: Large habitats should influence regional patterns more than small ones
2. **Functional differences**: Different habitat types contribute differently to ecosystem services
3. **Spatial autocorrelation**: Nearby observations are not independent
4. **Edge effects**: Boundaries between management units may not reflect ecological boundaries

## Framework Capabilities


### 1. Multi-Scale Spatial Data Management

**Automated Data Loading and Validation**

- Supports shapefiles, GeoPackage, and other spatial formats
- Automatic coordinate reference system standardization
- Geometry validation and repair for robust analysis
- Efficient memory management for large ecological datasets

**Hierarchical Scale Relationships**

```python
# Define spatial hierarchy
wetlands = SpatialScale("wetlands", "wetlands.shp", "WETLAND_ID", is_base_scale=True)
reaches = SpatialScale("reaches", "river_reaches.shp", "REACH_ID") 
catchments = SpatialScale("catchments", "catchments.shp", "CATCH_ID")
```

### 2. Ecologically-Informed Aggregation Methods

**Area-Weighted Aggregation**

- Appropriate for habitat quality metrics, vegetation indices
- Larger habitats contribute proportionally more to regional patterns
- Accounts for habitat extent in ecosystem service calculations

**Length-Weighted Aggregation**

- Designed for riparian corridor analysis, stream assessments
- Longer river segments have greater influence on reach-level patterns
- Maintains connectivity information in aggregated data

**Frequency-Weighted Analysis**

- Ideal for species occurrence data, management action frequency
- Accounts for relative abundance of different habitat types
- Preserves information about habitat diversity

**Custom Weighting Schemes**

- User-defined weights based on ecological importance
- Conservation priority weighting
- Threat level or management urgency weighting

### 3. Functional Group Classification

**Intelligent Reclassification**

```python
# Group wetland types for analysis using keyword search
wetland_groups = {
    'Permanent': ['Permanent Wetland', 'Deep Freshwater Marsh', 'Permanent Lake'],
    'Seasonal': ['Seasonal Wetland', 'Ephemeral Lake', 'Floodplain Wetland'],
    'Constructed': ['Farm Dam', 'Irrigation Channel', 'Constructed Wetland']
}
```

**Substring-Based Matching**

- Flexible classification using partial text matching
- Handles variations in naming conventions across datasets
- Maintains traceability of classification decisions

### 4. Statistical Validation and Quality Control

**Data Quality Assurance**

- Duplicate record detection and handling
- Missing data identification and reporting
- Outlier detection with ecological context
- Spatial topology validation


### 5. Visualization and Reporting

**Multi-Scale Visualization**

- Hierarchical map displays showing scale relationships
- Choropleth maps to visualise metrics across scales


## Requirements

python
pandas
geopandas
matplotlib
jupyter notebooks

```

# Install dependencies
pip install -r requirements.txt

```

## Quick Start

### Basic Workflow

```python
from scaling_framework_functions import SpatialScale

# 1. Load base scale data (individual wetlands)
wetlands = SpatialScale(
    name="wetlands",
    source="data/MDB_wetlands.shp",
    unique_id_fields="WETLAND_ID",
    type_field="WETLAND_TYPE",
    is_base_scale=True
)

# 2. Load aggregation scale (catchments)  
catchments = SpatialScale(
    name="catchments",
    source="data/MDB_catchments.shp", 
    unique_id_fields="CATCH_ID"
)

# 3. Aggregate wetland condition to catchment scale
wetlands.aggregate_to(
    target_scale=catchments,
    metric_columns=['condition_score', 'NDVI_2023'],
    method='area_weighted',
    result_name='wetland_condition'
)

# 4. Save results
catchments.save_results('output/', file_types=['.csv', '.gpkg'])
```

## API Reference

### SpatialScale Class

#### Constructor

```python
SpatialScale(name, source, unique_id_fields, weighting_field=None, 
            metric_fields=None, type_field=None, is_base_scale=False)
```

#### Key Methods

**aggregate_to(target_scale, metric_columns, method='area_weighted', ...)**

- Spatially join and aggregate data to target scale
- Methods: 'area_weighted', 'length_weighted', 'count', 'sum', 'frequency_weighted'

**join_data(data_path, pivot_columns=None, ...)**  

- Join external tabular data with spatial features
- Supports pivoting from long to wide format

**save_results(output_folder, file_types=['.csv'])**

- Export aggregated results to multiple formats
- Supports CSV, Shapefile, GeoPackage

**plot(result=None, field=None, reclass_map=None)**

- Generate maps and visualizations
- Multi-panel plots by functional groups

#### Aggregation Methods

| Method | Use Case | Weighting | Example |
|--------|----------|-----------|---------|
| `area_weighted` | Habitat metrics | Polygon area | Wetland condition scores |
| `length_weighted` | Corridor analysis | Line length | Riparian vegetation health |
| `count` | Site-based data | Equal weighting | Species occurrence records |
| `sum` | Additive quantities | None | Total habitat area |
| `frequency_weighted` | Compositional data | Type frequency | Habitat type diversity |

## Usage Examples

refer included jupyter notebooks.

`test_framework.ipynb` will install simple test data and run through some example calculations and plots.

`scaling_framework.ipynb` is an application of the framework in the Murray-Darling Basin to scale up from ANAE ecosystem polygons to larger scales represented by Ramsar sites, Directory of Important Wetlands (DIWA), Valleys, and Northern and Southern Basin regions


## Author adn Licence

Open Source  - Creative Commons v4 with attribution CCby4

Citation:


`Brooks (2025). MDB Scaling-up Framework: Flow-MER Program. Commonwealth Environmental Water Holder, Australian Government Department of Climate Change, Energy, the Environment and Water. Sourced from <github URL>`


Shane Brooks

![Brooks Logo](brooks-logo.png)

