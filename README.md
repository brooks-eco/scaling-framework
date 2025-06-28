# Generic Scaling-up Framework

A framework for hierarchical data aggregation across spatial scales with a specific focus on scaling up the evaluation of rivers, wetlands and floodplains within large river systems like the Murray-Darling Basin.

## Overview

A fundamental challenge in landscape ecology is **how to transfer ecological information across spatial scales while preserving biological relevance and statistical validity**.

This framework is a code library developed as a small research project within the [Flow-MER project](https://www.flow-mer.org.au/)

The framework is used to:

- Aggregate fine-resolution ecological data to larger management-relevant scales
- Maintains ecological meaning through appropriate weighting schemes
- Handles complex spatial hierarchies common in landscape ecology
- Ensure statistical rigor in calculation of weighted metrics to summaries larger areas from smaller scale building blocks.

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

Averaging fine-scale data to broader scales can be **ecologically misleading** because of:

1. **Unequal representation**: Large habitats should influence regional patterns more than small ones
2. **Functional differences**: Different habitat types contribute differently to ecosystem services
3. **Spatial autocorrelation**: Nearby observations are not independent
4. **Edge effects**: Boundaries between management units may not reflect ecological boundaries

## This Framework

The code uses spatial data sets (shapefiles or geopackages) and frames the spatial relationships between (e.g. smaller patches within larger regions).
Data attached to the smallest scale features (the building blocks) are then applied to larger scales using spatial joins and weighted averages.
Larger scales can then be aggregated up the chain to even larger scales - e.g. wetlands -> catchments -> whole river basin

### Metric aggregation

Metric aggregation from small scales to larger scales leverages the `weighted mean formula` in different ways.

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

`test_framework.ipynb` will install simple generated test data and run through some example calculations and plots.  This should run out of the box.  The helper script `create_test_data.py` generates some simple rectangles (small rectangles within larger rectangles) and includes examples of:

- using a reclass_map to reclassify the types
- plotting the spatial hierarchy
- loading some extra data from CSV
- calculating area and frequency weighted averages from the base scale up to region scale
- visualising the results

WORK IN PROGRESS - `scaling_framework.ipynb` is an application of the framework in the Murray-Darling Basin to scale up from ANAE ecosystem polygons to larger scales represented by Ramsar sites, Directory of Important Wetlands (DIWA), Valleys, and Northern and Southern Basin regions.
The source data sets are large and not included in this repository

- Brooks S (2021) Australian National Aquatic Ecosystem (ANAE) Classification of the Murray-Darling Basin v3.0. Wetland Polygons. Accessed 16 May 2023, [https://fed.dcceew.gov.au/datasets/1e57385ab8374f51b4b518a8cf571dbc/about](https://fed.dcceew.gov.au/datasets/1e57385ab8374f51b4b518a8cf571dbc/about)
- DAWE (2020) Ramsar Wetlands of Australia. Accessed 26 February 2021, [http://www.environment.gov.au/fed/catalog/search/resource/details.page?uuid=%7BF49BFC55-4306-4185-85A9-A5F8CD2380CF%7D](http://www.environment.gov.au/fed/catalog/search/resource/details.page?uuid=%7BF49BFC55-4306-4185-85A9-A5F8CD2380CF%7D)
- MDBA (2018) Basin-Wide Environmental Watering Strategy Regions for Vegetation Outcomes. Accessed 29 June 2025, [https://data.gov.au/data/dataset/basin-wide-environmental-watering-strategy-regions-for-vegetation-outcomes](https://data.gov.au/data/dataset/basin-wide-environmental-watering-strategy-regions-for-vegetation-outcomes).
- DIWA Directory of Important Wetlands, polygon shapefile (doesn't seem to be accessible any more), [https://www.dcceew.gov.au/water/wetlands/australian-wetlands-database](https://www.dcceew.gov.au/water/wetlands/australian-wetlands-database)

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

Brooks (2025). MDB Scaling-up Framework: Flow-MER Program. Commonwealth Environmental Water Holder, Australian Government Department of Climate Change, Energy, the Environment and Water. Sourced from [https://github.com/brooks-eco/scaling-framework](https://github.com/brooks-eco/scaling-framework)

Shane Brooks

![Brooks Logo](brooks-logo.png)
