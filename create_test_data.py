"""
Create test shapefiles for the scaling-up framework
"""
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
from pathlib import Path

# Create test directory
test_dir = Path("./test")
test_dir.mkdir(exist_ok=True)

# Define coordinate system
crs = "EPSG:4326"

# Create base scale rectangles (small units)
base_polygons = [
    Polygon([(0, 0), (2, 0), (2, 1), (0, 1)]),    # Bottom left
    Polygon([(2, 0), (4, 0), (4, 1), (2, 1)]),    # Bottom right
    Polygon([(0, 1), (2, 1), (2, 2), (0, 2)]),    # Top left
    Polygon([(2, 1), (4, 1), (4, 2), (2, 2)]),    # Top right
]

base_data = {
    'UID': ['B001', 'B002', 'B003', 'B004'],
    'HabitatType': ['Forest', 'Wetland', 'Forest', 'Grassland'],
    'Area_Ha': [200, 200, 200, 200],
    'NDVI': [0.8, 0.6, 0.75, 0.5],
    'geometry': base_polygons
}

# Ensure all data is properly formatted
base_gdf = gpd.GeoDataFrame(base_data, crs=crs)
# Force save all columns by explicitly setting dtypes
base_gdf['UID'] = base_gdf['UID'].astype(str)
base_gdf['HabitatType'] = base_gdf['HabitatType'].astype(str)
base_gdf['Area_Ha'] = base_gdf['Area_Ha'].astype(float)
base_gdf['NDVI'] = base_gdf['NDVI'].astype(float)

# Save with explicit schema to ensure all columns are preserved
base_gdf.to_file(test_dir / "base_units.shp", driver='ESRI Shapefile')

# Create aggregation scale 1 - regions (2 rectangles)
region_polygons = [
    Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]),    # Left region
    Polygon([(2, 0), (4, 0), (4, 2), (2, 2)]),    # Right region
]

region_data = {
    'RegionID': ['R001', 'R002'],
    'RegionName': ['West', 'East'],
    'geometry': region_polygons
}

region_gdf = gpd.GeoDataFrame(region_data, crs=crs)
region_gdf.to_file(test_dir / "regions.shp")

# Create aggregation scale 2 - basin (1 large rectangle)
basin_polygon = [Polygon([(0, 0), (4, 0), (4, 2), (0, 2)])]

basin_data = {
    'BasinID': ['BASIN001'],
    'BasinName': ['Test Basin'],
    'geometry': basin_polygon
}

basin_gdf = gpd.GeoDataFrame(basin_data, crs=crs)
basin_gdf.to_file(test_dir / "basin.shp")

# Create test CSV data
csv_data = {
    'UID': ['B001', 'B002', 'B003', 'B004'] * 3,
    'NDVI': [0.8, 0.6, 0.75, 0.5, 0.82, 0.58, 0.77, 0.52, 0.78, 0.62, 0.73, 0.48],
    'year': [2020, 2020, 2020, 2020, 2021, 2021, 2021, 2021, 2022, 2022, 2022, 2022]
}

csv_df = pd.DataFrame(csv_data)
csv_df.to_csv(test_dir / "test_NDVI_data.csv", index=False)

print("Test data created successfully!")
print(f"Files created in {test_dir}:")
print("- base_units.shp (4 rectangles)")
print("- regions.shp (2 regions)")
print("- basin.shp (1 basin)")
print("- test_NDVI_data.csv (temporal data)")