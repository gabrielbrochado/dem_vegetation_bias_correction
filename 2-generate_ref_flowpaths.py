# IMPORT LIBS

import os
import random
import rasterio
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, box, LineString, MultiLineString, Polygon
import pyproj
from pyproj import Transformer
from rasterio.features import shapes
from shapely.geometry import shape
from rasterio.features import shapes
from shapely.ops import split
import math
from shapely.validation import make_valid
from rasterio.windows import Window
import tempfile

# FUNCTIONS


def vectorize_forest_raster(raster_path):
    """
    Vectorize the forested areas from a raster file and return a GeoDataFrame of forest polygons.

    Parameters:
    - raster_path (str): Path to the raster file containing the forest mask.

    Returns:
    - forest_gdf (GeoDataFrame): GeoDataFrame containing the forest polygons.
    """
    # Load the raster data
    with rasterio.open(raster_path) as src:
        # Read the raster data
        raster = src.read(1)
        raster_affine = src.transform

        # Binarize the raster (forest height > 0)
        binary_raster = (raster > 0).astype(np.uint8)

        # Extract forest polygons from the binarized raster
        forest_polygons = []
        for shape_pair in shapes(binary_raster, mask=(binary_raster == 1), transform=raster_affine):
            geom, value = shape_pair
            if value == 1:
                forest_polygons.append(shape(geom))

        # Create a GeoDataFrame for the forest polygons
        forest_gdf = gpd.GeoDataFrame(geometry=forest_polygons, crs=src.crs)

    return forest_gdf

def calculate_forest_class(gdf, forest_gdf):
    """
    Classify each line in the GeoDataFrame based on whether the majority of its length intersects with forested areas
    and calculate the forest percentage.

    Parameters:
    - gdf (GeoDataFrame): GeoDataFrame containing line geometries.
    - forest_gdf (GeoDataFrame): GeoDataFrame containing forest polygons.

    Returns:
    - output_gdf (GeoDataFrame): GeoDataFrame with added 'class' and 'forest_percentage' fields.
    """
    # Initialize the 'class' column to 0 and 'forest_percentage' to 0 for all flowpaths
    gdf['class'] = 0
    gdf['forest_percentage'] = 0.0

    # Check if each flowpath intersects with any forest polygon
    for idx, flowpath in gdf.iterrows():
        # Calculate the intersection of the flowpath with forest polygons
        forest_intersections = forest_gdf.intersection(flowpath.geometry)

        # Convert flowpath to GeoSeries for CRS transformation and length calculation
        flowpath_series = gpd.GeoSeries([flowpath.geometry], crs=gdf.crs)
        total_length = flowpath_series.to_crs('EPSG:4087').length[0]

        # Calculate the total length of the intersected segments
        forest_length = forest_intersections.to_crs('EPSG:4087').length.sum()

        # Calculate the forest percentage
        forest_percentage = forest_length / total_length if total_length > 0 else 0

        # Store the forest percentage in the GeoDataFrame
        gdf.at[idx, 'forest_percentage'] = forest_percentage

        # Determine class based on the length percentage
        if forest_percentage > 0.5:
            gdf.at[idx, 'class'] = 1

    return gdf

def find_vertices_population(lines_gdf, output_path = None):

    ### CREATE VERTEX LIST WITHOUT END-POINT DUPLICATES ####

    # Load drainage lines
    gdf = lines_gdf
  
    # Ensure all geometries are appropriate line types
    if not all(gdf.geometry.type.isin(['LineString', 'MultiLineString'])):
        raise ValueError("All geometries must be LineStrings or MultiLineStrings")
    
    # Lists to store vertices and their line IDs
    start_points = []
    other_vertices = []
    start_point_ids = []
    other_vertex_ids = []
    
    # Extract vertices
    for index, row in gdf.iterrows():
        geom = row.geometry
    
        # If MultiLineString
        if isinstance(geom, MultiLineString):
            for line in geom.geoms:
                coords = list(line.coords)
                start_points.append(coords[0])
                start_point_ids.append(index)
                for coord in coords[1:]:
                    other_vertices.append(coord)
                    other_vertex_ids.append(index)
    
        # If LineString
        elif isinstance(geom, LineString):
            coords = list(geom.coords)
            start_points.append(coords[0])
            start_point_ids.append(index)
            for coord in coords[1:]:
                other_vertices.append(coord)
                other_vertex_ids.append(index)
    
    # Remove points in other_vertices that are in the start points list
    filtered_other = [(pt, id) for pt, id in zip(other_vertices, other_vertex_ids) if pt not in start_points]
    
    # Generate gdf with points
    all_vertices = start_points + [pt for pt, id in filtered_other]
    all_ids = start_point_ids + [id for pt, id in filtered_other]
    
    all_points_gdf = gpd.GeoDataFrame({
        'geometry': gpd.points_from_xy(*zip(*all_vertices)),
        'line_id': all_ids 
    })
    
    all_points_gdf = all_points_gdf.set_crs('EPSG:4326')
    
    #### fnid end points with dangles ####
    
    # Function to extract the end point of a directed line
    def extract_end_point(line):
        return line.boundary.geoms[-1]
    
    # Extract end points from all lines
    endpoints = [extract_end_point(line) for line in gdf.geometry]
    
    # Create a GeoDataFrame from endpoints
    endpoints_gdf = gpd.GeoDataFrame(geometry=endpoints, crs=gdf.crs)
    
    # Function to check if a point is a dangle
    def is_dangle(point, lines):
        intersecting_lines = lines[lines.intersects(point)]
        return len(intersecting_lines) == 1
    
    # Identify dangle points
    dangle_points = [point for point in endpoints_gdf.geometry if is_dangle(point, gdf)]
    
    # Create a GeoDataFrame from dangle points
    dangle_points_gdf = gpd.GeoDataFrame(geometry=dangle_points, crs=gdf.crs)
    
    ### REMOVE END-POINT DANGLES #####
    
    common_points = gpd.sjoin(all_points_gdf, dangle_points_gdf, how="inner", predicate="intersects")
    filtered_points_gdf = all_points_gdf.drop(common_points.index)
    
    return filtered_points_gdf

def create_box_around_point(point, radius_in_meters):
    # Create a transformer to convert from EPSG:4326 to a local projection (UTM)
    transformer_to_utm = Transformer.from_crs("EPSG:4326", "EPSG:4087", always_xy=True)
    transformer_to_wgs84 = Transformer.from_crs("EPSG:4087", "EPSG:4326", always_xy=True)

    # Convert the point to UTM
    point_utm = transformer_to_utm.transform(point.x, point.y)
    
    # Create a Point object in UTM coordinates
    point_utm_obj = Point(point_utm)

    # Create a circle around the point with the specified radius in UTM coordinates
    circle_utm = point_utm_obj.buffer(radius_in_meters)

    # Convert the circle back to EPSG:4326
    circle_wgs84 = gpd.GeoSeries([circle_utm]).set_crs("EPSG:4087").to_crs("EPSG:4326").geometry.iloc[0]

    return circle_wgs84

def drop_z(geometry):
    
    if geometry.has_z:
        # Create a new LineString from the XY parts of each coordinate
        return LineString([(x, y) for x, y, z in geometry.coords])
    return geometry

def construct_downslope_path(line, vertex_point, gdf):
    # Split the line using the vertex
    split_line = split(line, vertex_point)
    
    # The downslope part is the part after the vertex
    downslope_part = split_line.geoms[-1]
    
    # Function to find the next downslope line
    def find_next_line(current_line, gdf):
        end_point = Point(current_line.coords[-1])
        touching_lines = gdf[gdf.touches(end_point)]
        for index, row in touching_lines.iterrows():
            start_point = Point(row.geometry.coords[0])
            if end_point.equals(start_point):
                return row.geometry
        return None
    
    # Construct the full downslope path
    full_downslope_path = [downslope_part]
    current_line = downslope_part
    while True:
        next_line = find_next_line(current_line, gdf)
        if next_line is None or next_line.equals(current_line):
            break
        full_downslope_path.append(next_line)
        current_line = next_line
    
    # Combine all parts into a single LineString
    combined_path = LineString([point for line in full_downslope_path for point in line.coords])
    
    return combined_path

# PROCESS

base_path = 'C:/Users/gabro/Downloads' # <- the one used to unzip the suplementary materials zip

for area in [1,2,3,4]:

    input_shp = base_path + f"/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/HID_Trecho_Drenagem_L_merge_clean_strahler.shp"
    prefix = input_shp.split("/")[:-1]
    drainage_gdf = gpd.read_file(input_shp).to_crs('EPSG:4326')
    drainage_gdf = drainage_gdf.explode(ignore_index=True)
    drainage_gdf['geometry'] = drainage_gdf['geometry'].apply(drop_z)
    drainage_gdf['line_id'] =  drainage_gdf.index
    vertices_pop = find_vertices_population(drainage_gdf)
    
    results = []

    if area == 1:
        raster_path = base_path + f"/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/aux/adjusted_fh2020_align_12.tif"
    elif area == 2:
        raster_path = base_path + f"/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/aux/adjusted_fh2020_align_13.tif"
    elif area == 3:
        raster_path = base_path + f"/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/aux/adjusted_fh2020_align_10.tif"
    elif area == 4:
        raster_path = base_path + f"/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/aux/adjusted_fh2020_align_13.tif"

    forest_gdf = vectorize_forest_raster(raster_path)
    
    
    for radius in [1000, 2000, 3000]:
    
        # Set the random seed for reproducibility
        seed_value = 42
        random.seed(seed_value)
        np.random.seed(seed_value)
        
        clipped_paths_gdf =  clipped_paths_gdf =  gpd.GeoDataFrame(columns=['geometry'], geometry='geometry')
        random_vertices_gdf = gpd.GeoDataFrame()
        max_fails = 500
        fail_count = 0
    
        while fail_count <= max_fails:
            # Sample a point
            sample_point = vertices_pop.sample(1)
            point_geom = sample_point.geometry.values[0]
            line_id = sample_point.line_id.values[0]
    
            # Create box around point
            window = create_box_around_point(point_geom, radius)
            window_gdf = gpd.GeoDataFrame({'geometry': [window]}, crs='EPSG:4326')
            
            # FClip
            clipped_drainage_gdf = gpd.clip(drainage_gdf, window_gdf).explode(ignore_index= True)
    
            # Find line on the clipped drainage that intersects the point and have the same line id
            line_gdf = clipped_drainage_gdf[(clipped_drainage_gdf.geometry.intersects(point_geom)) & (clipped_drainage_gdf.line_id == line_id)]
            line = line_gdf.geometry.values[0]
    
            # Find downslope path
            path = construct_downslope_path(line, point_geom, clipped_drainage_gdf)
    
            # Create gdf with the path
            path_gdf = gpd.GeoDataFrame({'geometry' : [path]}, geometry = 'geometry').set_crs('EPSG:4326')
                
            # Calculate the number of intersections with existing paths in 'clipped_paths_gdf'
            n_intersections = len(clipped_paths_gdf[clipped_paths_gdf.geometry.intersects(path)])
    
            # Check if n_intersection == 0 and path intersects the window border
            if (n_intersections == 0)  and (path.intersects(window.boundary)):
                clipped_paths_gdf = pd.concat([clipped_paths_gdf, path_gdf])
                random_vertices_gdf = pd.concat([random_vertices_gdf, sample_point])
                fail_count = 0
    
            else:
                fail_count = fail_count + 1
                continue
    
        random_vertices_gdf['point_id'] = range(1, len(clipped_paths_gdf)+1)
        clipped_paths_gdf['point_id'] = range(1, len(clipped_paths_gdf)+1)
    
        clipped_paths_gdf.to_file(os.path.dirname(input_shp) + f"/b{radius}_ref_flowpaths_r_class_p.gpkg", driver="GPKG")

        gdf_path = os.path.dirname(input_shp) + f"/b{radius}_ref_flowpaths_r_class_p.gpkg"
        gdf = gpd.read_file(gdf_path)
        
        # Run the function with precomputed forest polygons
        output_gdf = calculate_forest_class(gdf, forest_gdf)
        
        # Save the output GeoDataFrame
        output_gdf_path = os.path.dirname(input_shp) + f"/b{radius}_ref_flowpaths_r_class_p.gpkg"
        output_gdf.to_file(output_gdf_path, driver='GPKG', mode='w')