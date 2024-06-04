from io import StringIO
import ee
import folium

start_date = '2020-07-24'
end_date = '2020-07-25'
# Initialize the Earth Engine module.
ee.Initialize()

def add_ee_layer(self, ee_object, vis_params, name):
    """Function to add Earth Engine layers to Folium map."""
    try:
        if isinstance(ee_object, ee.image.Image):
            map_id_dict = ee.Image(ee_object).getMapId(vis_params)
        elif isinstance(ee_object, ee.imagecollection.ImageCollection):
            ee_object = ee_object.mosaic()
            map_id_dict = ee.Image(ee_object).getMapId(vis_params)
        elif isinstance(ee_object, ee.geometry.Geometry):
            ee_object = ee.FeatureCollection([ee.Feature(ee_object)])
            map_id_dict = ee_object.getMapId(vis_params)
        elif isinstance(ee_object, ee.featurecollection.FeatureCollection):
            map_id_dict = ee_object.getMapId(vis_params)
        
        folium.raster_layers.TileLayer(
            tiles=map_id_dict['tile_fetcher'].url_format,
            attr='Map data &copy; <a href="https://earthengine.google.com/">Google Earth Engine</a>',
            name=name,
            overlay=True,
            control=True
        ).add_to(self)
    except:
        print("Could not add layer to map.")

folium.Map.add_ee_layer = add_ee_layer

# Define your area of interest (AOI) with coordinates.
aoi = ee.Geometry.Rectangle([-120, 49, -95.00, 60])

# Load the Landsat 8 ImageCollection.
l8_collection = ee.ImageCollection('LANDSAT/LC08/C01/T1') \
    .filterDate(start_date, end_date) \
    .filterBounds(aoi) \
    .sort('CLOUD_COVER') \
    .first()

# Calculate NDVI.
ndvi = l8_collection.normalizedDifference(['B5', 'B4']).rename('NDVI')

# Visualization parameters.
vis_params = {
    'min': 0,
    'max': 1,
    'palette': ['blue', 'white', 'green']
}

# Set the location for the center of the map.
map_center = aoi.centroid().coordinates().getInfo()[::-1]  # Reverse to match [lat, lon] format.

# Create a Folium map centered at the AOI.
my_map = folium.Map(location=map_center, zoom_start=5)

# Add the NDVI layer to the map.
my_map.add_ee_layer(ndvi, vis_params, 'NDVI')

# Add a layer control panel to the map.
my_map.add_child(folium.LayerControl())


# Showing map
my_map




########################## CHRS Precipitation ########################


import netCDF4 as nc

# Replace with the path to your NetCDF file
file_path = '/users/your/localfile/PDIR_2024-06-03054028pm/PDIR_2024-06-03054028pm.nc'

# Open the NetCDF file
dataset = nc.Dataset(file_path, 'r')

# Print the dataset details
print(dataset)

# List all variables
print("\nVariables in the dataset:")
for var in dataset.variables:
    print(var)

# Access specific variables (replace 'variable_name' with the actual variable name)
# Example: variable_name = 'temperature'
variable_name = 'your_variable_name'
if variable_name in dataset.variables:
    data = dataset.variables[variable_name][:]
    print(f"\nData for variable '{variable_name}':\n", data)
else:
    print(f"\nVariable '{variable_name}' not found in the dataset.")

# Close the dataset
dataset.close()


import xarray as xr

net_cdf=xr.open_dataset(file_path)


df=net_cdf.to_dataframe()

df=df.reset_index()

df=df.drop(columns='datetime')

import geopandas as gpd

gpd.GeoDataFrame(df, geometry = gpd.points_from_xy(x = df['lon'], y=df['lat'], crs='EPSG:4326'))\
    .explore(column='precip', legend =True)



######
