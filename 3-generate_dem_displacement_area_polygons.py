# IMPORT LIBS

import os
import pandas as pd
import geopandas as gpd
import rasterio
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.validation import make_valid
from shapely.geometry import Point
from pyproj import Transformer
import math
from shapely.geometry import LineString

# INPUTS
base_path = 'C:/Users/gabro/Downloads'  # <- the one used to unzip the supplementary materials zip
output_filename = base_path + '/Suplementary_Materials_Brochado_Renno_2024/stacked_dataframe_r_class_p.csv'

# FUNCTIONS

def extract_starting_points(lines_gdf):
    """
    Create a new GeoDataFrame with the starting points of each line in the input GeoDataFrame.

    Parameters:
    -----------
    lines_gdf : GeoDataFrame
        GeoDataFrame containing line geometries.

    Returns:
    --------
    start_points_gdf : GeoDataFrame
        New GeoDataFrame with starting points (Point geometries) of each line in the input GeoDataFrame.
    """
    # Extract the starting point of each line geometry
    start_points = lines_gdf.geometry.apply(lambda line: Point(line.coords[0]))

    # Create a new GeoDataFrame with the same attributes but replace the geometry with starting points
    start_points_gdf = lines_gdf.copy()
    start_points_gdf['geometry'] = start_points

    return start_points_gdf



def get_next_pixel(row, col, direction):
    directions = {
        1: (0, 1),    # E
        2: (1, 1),    # SE
        4: (1, 0),    # S
        8: (1, -1),   # SW
        16: (0, -1),  # W
        32: (-1, -1), # NW
        64: (-1, 0),  # N
        128: (-1, 1)  # NE
    }
    move = directions.get(direction, (0, 0))
    return row + move[0], col + move[1]



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


def point_ldd_flowpath(ldd_array, start_point, transform):
    # Convert start point to row, col
    row, col = rasterio.transform.rowcol(transform, start_point.x, start_point.y)
    
    # Ensure starting point is added as the first point in the flowpath
    flowpath_coords = [(start_point.x, start_point.y)]

    num_rows, num_cols = ldd_array.shape

    while True:
        if not (0 <= row < num_rows and 0 <= col < num_cols):
            break  # Ensure the indices are within the raster bounds

        direction = ldd_array[row, col]
        if direction == 0:
            break

        next_row, next_col = get_next_pixel(row, col, direction)
        
        x_next, y_next = rasterio.transform.xy(transform, next_row, next_col)
        next_point = (x_next, y_next)
        
        flowpath_coords.append(next_point)

        row, col = next_row, next_col

    # Forcefully make the starting point of the LineString the input Point
    flowpath_coords[0] = (start_point.x, start_point.y)

    # Handle the case where flowpath_coords has only one element creating a dummy short line
    if len(flowpath_coords) == 1:
        # Create a second point 1 meter away in the x-axis
        x_new = flowpath_coords[0][0] + 1 / 111320  # 1 meter in degrees, approx conversion
        y_new = flowpath_coords[0][1]
        flowpath_coords.append((x_new, y_new))

    # Create LineString from flowpath coordinates
    flowpath_line = LineString(flowpath_coords)
    
    return flowpath_line


def points_gdf_clipped_ldd_flowpaths(gdf, ldd_array, transform, box_size):
    flowpaths = []
    threshold = 5 * (1 / 111320)
    for idx, row in gdf.iterrows():
        start_point = row.geometry
        point_id =  row.point_id
        flowpath = point_ldd_flowpath(ldd_array, start_point, transform)

        # if the flowpath is a short dummy line skip this iteration
        if flowpath.length <= threshold:
            continue
        
        # Create GeoDataFrame for the single flowpath
        single_flowpath_gdf = gpd.GeoDataFrame(geometry=[flowpath], crs=gdf.crs)
        
        # Create box around point
        window = create_box_around_point(start_point, box_size)
        window_gdf = gpd.GeoDataFrame({'geometry': [window]}, crs='EPSG:4326')
        
        # Clip the flowpath using the window
        clipped_flowpath_gdf = gpd.clip(single_flowpath_gdf, window_gdf).explode(ignore_index=True)
        
        # Choose the line that intersects the starting point (multlines can be created while clipping)
        clipped_flowpath_gdf = clipped_flowpath_gdf[clipped_flowpath_gdf.geometry.intersects(start_point)]

        # Copy point_id
        clipped_flowpath_gdf['point_id'] = point_id
        
        
        # Append the valid clipped flowpath to the list
        if not clipped_flowpath_gdf.empty:
            flowpaths.append(clipped_flowpath_gdf)
    
    # Create a new GeoDataFrame with the flowpaths
    flowpaths_gdf = gpd.GeoDataFrame(pd.concat(flowpaths, ignore_index=True), crs='EPSG:4326')
    
    return flowpaths_gdf


def generate_polygons_between_lines(gdf1, gdf2):
    # Ensure both GeoDataFrames have 'point_id' field
    if 'point_id' not in gdf1.columns or 'point_id' not in gdf2.columns:
        raise ValueError("Both GeoDataFrames must have a 'point_id' field.")
    
    polygons = []

    for idx, row in gdf1.iterrows():
        point_id = row['point_id']
        line1 = row.geometry
        classification =  row['class']
        # mean_strahler =  row['mean_strahler']
        
        # Find the matching line in gdf2
        matching_row = gdf2[gdf2['point_id'] == point_id]
        if not matching_row.empty:
            line2 = matching_row.iloc[0].geometry

            coords1 = list(line1.coords)
            coords2 = list(line2.coords)
            polygon_coords = coords1 + coords2[::-1]
            polygon = Polygon(polygon_coords)

            # Make the polygon valid (to deal with self-intersecting polygon)
            if not polygon.is_valid:
                polygon = make_valid(polygon)
            
            # Store polygon as gdf
            single_polygon_gdf = gpd.GeoDataFrame([{'geometry': polygon}], crs=gdf1.crs)

            # Calculate area between lines
            single_polygon_gdf['area'] = single_polygon_gdf.explode(index_parts=True).to_crs('EPSG:6933').area.sum()
            single_polygon_gdf['point_id'] = point_id
            single_polygon_gdf['class'] = classification
            
            polygons.append(single_polygon_gdf)
        else:
            print(f"Warning: No matching line found in gdf2 for point_id {point_id}")

    # Combine all single polygon GeoDataFrames into one
    polygons_gdf = gpd.GeoDataFrame(pd.concat(polygons, ignore_index=True), crs=gdf1.crs)
    
    return polygons_gdf



def find_ldd_polygons(ldd_path, 
                      box_size, 
                      ref_paths_gdf,
                      flowpaths_file = None, 
                      polygons_file = None):
    
    ldd_raster = rasterio.open(ldd_path)
    ldd_array = ldd_raster.read(1)
    transform = ldd_raster.transform

    vertices_gdf = extract_starting_points(ref_paths_gdf)
    
    flowpaths_gdf = points_gdf_clipped_ldd_flowpaths(vertices_gdf, ldd_array, transform, box_size)
    polygons_gdf = generate_polygons_between_lines(ref_paths_gdf, flowpaths_gdf)
    
    if flowpaths_file:
        flowpaths_gdf.to_file(flowpaths_file, driver="GPKG")
    if polygons_file:     
        polygons_gdf.to_file(polygons_file, driver="GPKG")
    
    return polygons_gdf

def add_suffix_to_path(file_path, suffix):
    directory, filename = os.path.split(file_path)
    base, ext = os.path.splitext(filename)
    new_filename = f"{base}_{suffix}{ext}"
    return os.path.join(directory, new_filename)

def th_hydro_dem(burned_dem, hydro_dem):
    # Define Terrahidro path
    th_path = 'C:/Users/gabro/TerraHidro-5.2.0'

    # Execute remove pits
    cmd = rf'cd\ && cd {th_path} && th removepits "{burned_dem}" "{hydro_dem}"'
    print(cmd)
    flag = os.system(cmd)

    if flag != 0:
        print('Error executing th removepits command.')
    else:
        print('Hydro DEM done.')

def th_d8(hydro_dem, hydro_dem_d8):
    # Define Terrahidro path
    th_path = 'C:/Users/gabro/TerraHidro-5.2.0'

    # Execute d8
    cmd = rf'cd\ && cd {th_path} && th d8 "{hydro_dem}" "{hydro_dem_d8}"'
    flag = os.system(cmd)

    if flag != 0:
        print('Error executing th d8 command.')
    else:
        print('Calculating d8 done.')

def process_dem(burned_dem):
    # Generate the output filenames with the appropriate suffixes
    hydro_dem = add_suffix_to_path(burned_dem, "hydro")
    hydro_dem_d8 = add_suffix_to_path(hydro_dem, "d8")

    # Run the th_hydro_dem function
    th_hydro_dem(burned_dem, hydro_dem)

    # Run the th_d8 function
    th_d8(hydro_dem, hydro_dem_d8)
    print('hydro_dem_d8')
    return hydro_dem_d8

# PROCESS

# Initialize an empty list to store DataFrames for stacking
gdf_list = []

for area in [1,2,3,4]:  # Iterate through the desired areas
    for radius in [1000, 2000, 3000]:  # Iterate through the desired radii
        # Define file paths based on the area and radius
        glo = f'{base_path}/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/glo-30.tif'
        fabdem = f'{base_path}/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/fabdem.tif'
        new_corrected_dem = f'{base_path}/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/new_corrected_dem.tif'
        ref_flowpaths_path = f'{base_path}/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/b{radius}_ref_flowpaths_r_class_p.gpkg'

        # PROCESS
        
        ## Find DEM LDDs
        glo_ldd = process_dem(glo)
        fabdem_ldd = process_dem(fabdem)
        glo_corrected_ldd = process_dem(new_corrected_dem)
        
        ## Construct flowpaths from the reference flowpaths
        ref_flowpaths_gdf = gpd.read_file(ref_flowpaths_path)
        
        glo_polygons = find_ldd_polygons(
            glo_ldd,
            radius,
            ref_flowpaths_gdf,
            polygons_file=ref_flowpaths_path.split("_ref")[0] + f"_glo_polygons_r_class_p.gpkg",
        )

        fabdem_polygons = find_ldd_polygons(
            fabdem_ldd,
            radius,
            ref_flowpaths_gdf,
            polygons_file=ref_flowpaths_path.split("_ref")[0] + f"_fabdem_polygons_r_class_p.gpkg",
        )

        glo_corrected_polygons = find_ldd_polygons(
            glo_corrected_ldd,
            radius,
            ref_flowpaths_gdf,
            polygons_file=ref_flowpaths_path.split("_ref")[0] + f"_corrected_polygons_r_class_p.gpkg",
        )

        # Finding common valid IDs across all DataFrames
        valid_ids = list(
            set(glo_polygons['point_id'].values)
            .intersection(fabdem_polygons['point_id'].values)
            .intersection(glo_corrected_polygons['point_id'].values)
        )
        
        # Filter and reset index for each DataFrame using valid IDs
        glo_filtered = glo_polygons[glo_polygons['point_id'].isin(valid_ids)].reset_index(drop=True)
        glo_corrected_filtered = glo_corrected_polygons[glo_corrected_polygons['point_id'].isin(valid_ids)].reset_index(drop=True)
        fabdem_filtered = fabdem_polygons[fabdem_polygons['point_id'].isin(valid_ids)].reset_index(drop=True)

        # Combine into a single DataFrame for the current area and radius
        areas_gdf = pd.DataFrame({
            'point_id': glo_filtered['point_id'],
            'glo_corrected_area': glo_corrected_filtered['area'],
            'glo_area': glo_filtered['area'],
            'fabdem_area': fabdem_filtered['area'],
            'class': glo_filtered['class'],
            'area': area,  # Add area and radius columns for identification
            'radius': radius
        })

        # Append the resulting DataFrame to the list
        gdf_list.append(areas_gdf)

# Stack all DataFrames vertically
stacked_areas_gdf = pd.concat(gdf_list, ignore_index=True)

# Save the stacked DataFrame to CSV
stacked_areas_gdf.to_csv(output_filename, index=False)

print(f"Saved stacked DataFrame to: {output_filename}")
