# Generic Scaling-up Framework

A framework for hierarchical data aggregation across spatial scales with a specific focus on scaling up the evaluation of rivers, wetlands and floodplains within large river systems like the Murray-Darling Basin.

## Overview

The Scaling-up Framework addresses a fundamental challenge in landscape ecology: **how to transfer ecological information across spatial scales while preserving biological relevance and statistical validity**.

This framework:

- Aggregate fine-resolution ecological data to management-relevant scales
- Maintains ecological meaning through appropriate weighting schemes
- Handles complex spatial hierarchies common in landscape ecology
- Ensure statistical rigor in multi-scale analyses

### The Scale Challenge

In landscape ecology, **scale** refers to both the spatial extent (how large an area) and resolution (how fine the detail) of ecological observations.

Ecological systems commonly exhibit **hierarchical organization** where smaller units nest within larger ones:

```[]
Individual Wetlands → River Reaches → Sub-catchments → Catchments → Bioregions
    (Fine scale)                                                    (Broad scale)
```

Issues of scale complicate the study and management of natural systems because:

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

### Aggregation Challenges

Simply averaging fine-scale data to broader scales can be **ecologically misleading** because:

1. **Unequal representation**: Large habitats should influence regional patterns more than small ones
2. **Functional differences**: Different habitat types contribute differently to ecosystem services
3. **Spatial autocorrelation**: Nearby observations are not independent
4. **Edge effects**: Boundaries between management units may not reflect ecological boundaries

## This Framework

The code uses spatial data sets (shapefiles or geopackages)  to construct a spatial hierarchies.
Data attached to the smallest scale features (the atomic building blocks) are then applied to larger scales (using spatial joins) and weighted averages.

### Metric aggregation

Uses common weighting schemes that leverage the weighted mean formula

```[]
        Aggregated_Value = Σ(Metric_i × Weight_i) / Σ(Weight_i)

        Where weights vary by method:
        - Area-weighted: Weight_i = Area_i (hectares)
        - Length-weighted: Weight_i = Length_i (meters)
        - Count-based: Weight_i = 1 (equal weighting)
        - Sum (of attributes): Weight_i = metric_value_i
        - Frequency-weighted: Weight_i = Frequency_i / Total_frequency
```

| Method | Use Case | Weighting | Example |
|--------|----------|-----------|---------|
| `area_weighted` | Habitat metrics | Polygon area | Wetland or vegetation condition scores. e.g. larger habitats contribute proportionally more to regional patterns|
| `length_weighted` | Corridor analysis | Line length | Riparian vegetation health e.g. to resolve local river segments up to reach-level patterns |
| `count` | Site-based data | Equal weighting | Species occurrence records, management action frequency |
| `sum` | Additive quantities | None | User-defined weights e.g. conservation priority, habitat diversity or abundance measures, risk, priority |
| `frequency_weighted` | Compositional data | Type frequency | Habitat type diversity, management action frequency|

### Data Loading and Validation

- Supports shapefiles, GeoPackage, and other spatial formats to define the spatial scales
- Geometry validation and repair
- Selectable intersection rules (intersects, within, contains etc)
- Load additional data from CSV to attach to the spatial entities in wide cross-tabular or  long format.

### 3. Functional Group Classification

- A reclassification dictionary can be applied to the loaded data to generate new types in which to group the data. It can be useful to simplify complex data sets.  We use this to aggregate ecosystem type classes into broader subsets e.g. different types of forest becomes 'trees', or lakes, swamps and marshes could be grouped to a single 'wetland' class
- Handles variations in naming conventions across datasets
- Maintains traceability of classification decisions

### 4. Visualization and Reporting

- some simple helper routines are provided to visualise scale boundary relationships and metric values across scales.

## Quick Start

refer included jupyter notebooks.

`test_framework.ipynb` will install simple test data and run through some example calculations and plots.

`scaling_framework.ipynb` is an application of the framework in the Murray-Darling Basin to scale up from ANAE ecosystem polygons to larger scales represented by Ramsar sites, Directory of Important Wetlands (DIWA), Valleys, and Northern and Southern Basin regions

## Requirements

python
pandas
geopandas
matplotlib
jupyter notebooks

```[]

# Install dependencies
pip install -r requirements.txt

```

## Author adn Licence

Open Source  - Creative Commons v4 with attribution CCby4

Citation:

`Brooks (2025). MDB Scaling-up Framework: Flow-MER Program. Commonwealth Environmental Water Holder, Australian Government Department of Climate Change, Energy, the Environment and Water. Sourced from <github URL>`

Shane Brooks

![Brooks Logo](brooks-logo.png)
