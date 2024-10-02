# INPUTS

# Input paths           
base_path = "C:/Users/gabro/Downloads/teste8"
raw_dem_path =  f'{base_path}/glo-30.tif'
water_mask_path = f'{base_path}/glo-30_wbm.tif'
shapefile_path = f'{base_path}/extent.gpkg'

# Output paths
output_path_final = f'{base_path}/new_corrected_dem.tif'
aux_path = f'{base_path}/aux'


# INSTALL LIBS

!pip install geopandas --user
!pip install rasterio
!pip install geedim
!pip install scikit-image

# IMPORT LIBS

import ee
import os
import math
import shutil
import geemap
import geedim
import folium
import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np
import tempfile
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from scipy.ndimage import label as label_function
from scipy.ndimage import grey_dilation, convolve
from rasterio.merge import merge
import tempfile

# FUNCTIONS

## General Functions


def read_geotiff(file_path):
    with rasterio.open(file_path) as src:
        data = src.read(1)
        profile = src.profile
    return data, profile

# Function to write data to a GeoTIFF file
def write_geotiff(file_path, data, profile):
    with rasterio.open(file_path, 'w', **profile) as dst:
        dst.write(data, 1)
        
import os
from osgeo import gdal
import rasterio

def get_raster_info(raster_path):
    with rasterio.open(raster_path) as src:
        bounds = src.bounds
        xres = src.res[0]
        yres = src.res[1]
        srs = src.crs.to_string()
    return bounds, xres, yres, srs

def align_raster(target_path, ref_path, resampling, dtype, output_path):
    # Get the information from the reference raster using rasterio
    bounds, xres, yres, srs = get_raster_info(ref_path)

    # If output file exists, delete it (to handle overwriting)
    if os.path.exists(output_path):
        os.remove(output_path)

    # Open the target raster
    src_ds = gdal.Open(target_path)
    if src_ds is None:
        raise ValueError("Could not open source raster file")

    # Create the options string for gdalwarp
    warp_options = gdal.WarpOptions(format='GTiff',
                                    outputBounds=[bounds.left, bounds.bottom, bounds.right, bounds.top],
                                    xRes=xres, yRes=yres,
                                    dstSRS=srs,
                                    resampleAlg=resampling,
                                    outputType=gdal.GetDataTypeByName(dtype))

    # Perform the warp operation
    result = gdal.Warp(output_path, src_ds, options=warp_options)
    if result is None:
        raise ValueError("Error executing gdal_warp.")

    # Clean up
    result = None
    src_ds = None

import numpy as np
import rasterio
from scipy.ndimage import gaussian_gradient_magnitude

def custom_mean_filter(data, kernel_size=3):
    """
    Apply a mean filter to a 2D array using an efficient approach,
    handling NaN values by ignoring them in mean calculations. This version
    preserves the original edge values by not applying the mean filter to them.

    Parameters:
    - data: 2D numpy array to filter.
    - kernel_size: Size of the square kernel. Must be an odd integer.

    Returns:
    - filtered_data: 2D numpy array after applying the mean filter,
                     with original edge values preserved.
    """
    from numpy.lib.stride_tricks import sliding_window_view
    
    if kernel_size % 2 == 0:
        raise ValueError("Kernel size must be an odd integer.")
    
    # Initialize filtered_data array with original data to preserve edges
    filtered_data = data.copy().astype(float)
    
    # Create sliding windows
    window_shape = (kernel_size, kernel_size)
    windows = sliding_window_view(data, window_shape)
    
    # Calculate the mean, ignoring NaNs, across the windows
    mean_values = np.nanmean(windows, axis=(-2, -1))
    
    # Calculate start index for placing mean_values back into filtered_data
    start_idx = kernel_size // 2
    
    # Place mean_values into filtered_data, preserving the original edge values
    filtered_data[start_idx:-start_idx, start_idx:-start_idx] = mean_values
    
    return filtered_data

## Data Download Functions


def ee_export_image(image, area, output_path, nrows, ncols):
    # Initialize a list to hold the downloaded tiles
    tiles = []

    # Get bounds of the area of interest
    bounds = area.bounds()
    coords = bounds.coordinates().getInfo()[0]  # Extracting coordinates from nested list
    xmin, ymin = coords[0]
    xmax, ymax = coords[2]
    xdist = (xmax - xmin) / ncols
    ydist = (ymax - ymin) / nrows

    # Download DEM for each tile
    for i in range(nrows):
        for j in range(ncols):
            # Define tile region
            xmin_tile = xmin + j * xdist
            xmax_tile = xmin + (j + 1) * xdist
            ymin_tile = ymax - (i + 1) * ydist
            ymax_tile = ymax - i * ydist
            region = ee.Geometry.Rectangle([xmin_tile, ymin_tile, xmax_tile, ymax_tile])

            # Create temporary file
            temp_file = tempfile.mktemp(suffix='.tif')

            # Download DEM for the current tile
            geemap.ee_export_image(image, filename=temp_file, region=region)

            # Add the downloaded tile to the list
            tiles.append(temp_file)

    # Merge the downloaded tiles
    merge_tiles(tiles, output_path)

    # Remove temporary files
    for tile in tiles:
        os.remove(tile)

def merge_tiles(tiles, output_path):
    src_files_to_mosaic = [rasterio.open(tile) for tile in tiles]

    # Merge tiles
    mosaic, out_trans = merge(src_files_to_mosaic)

    # Update metadata
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({"driver": "GTiff",
                     "height": mosaic.shape[1],
                     "width": mosaic.shape[2],
                     "transform": out_trans})

    # Write the mosaic to disk
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(mosaic)


def download_image_geemap(image, features, out_path, crs=None, crs_transform=None, scale=None, resampling='near', dtype=None, overwrite=True, num_threads=None, max_tile_size=None, max_tile_dim=None, shape=None, scale_offset=False, unmask_value=None, column=None, **kwargs):
    """
    Wrapper function that takes a full file path and uses the download_ee_image_tiles function.

    Parameters:
        image (ee.Image): The image to be downloaded.
        features (ee.FeatureCollection): The features to loop through to download image.
        full_file_path (str): The full file path for the output file.
        crs (str): Reproject image(s) to this EPSG or WKT CRS.
        crs_transform (list): Affine transform in the specified CRS.
        scale (float): Resample image(s) to this pixel scale (size) (m).
        resampling (str): Resampling method, can be 'near', 'bilinear', 'bicubic', or 'average'.
        dtype (str): Convert to this data type.
        overwrite (bool): Overwrite the destination file if it exists.
        num_threads (int): Number of tiles to download concurrently.
        max_tile_size (int): Maximum tile size (MB).
        max_tile_dim (int): Maximum tile width/height (pixels).
        shape (tuple): (height, width) dimensions to export (pixels).
        scale_offset (bool): Whether to apply any EE band scales and offsets to the image.
        unmask_value (float): The value to use for pixels that are masked in the input image.
        column (str): The column name to use for the filename.

    Returns:
        None
    """
    # Extract the directory and file name from the full file path
    out_dir, file_name = os.path.split(out_path)
    prefix, _ = os.path.splitext(file_name)

    # Ensure the output directory exists
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Call the original function with the modified parameters
    geemap.download_ee_image_tiles(
        image, features, out_dir=out_dir, prefix=prefix, crs=crs, crs_transform=crs_transform, scale=scale,
        resampling=resampling, dtype=dtype, overwrite=overwrite, num_threads=num_threads,
        max_tile_size=max_tile_size, max_tile_dim=max_tile_dim, shape=shape,
        scale_offset=scale_offset, unmask_value=unmask_value, column=column, **kwargs
    )

    # Rename the output file to the specified full file path
    output_filepath_with_suffix = os.path.join(out_dir, f"{file_name.split('.')[0]}1.tif")
    if os.path.exists(output_filepath_with_suffix):
        os.rename(output_filepath_with_suffix, out_path)
    else:
        raise FileNotFoundError(f"The expected output file {output_filepath_with_suffix} does not exist.")

def download_glo_dem(roi, out_path):

    glo = ee.ImageCollection("COPERNICUS/DEM/GLO30").mosaic().setDefaultProjection('EPSG:4326', None, 30)
    image= glo.select(0)

    feature = ee.Feature(roi)
    features= ee.FeatureCollection([feature])

    download_image_geemap(image,
                                features,
                                out_path= out_path,
                                crs= 'EPSG:4326',
                                overwrite=True,
                                scale= 30)

def download_glo_rivers_and_lakes(roi, out_path):

    glo = ee.ImageCollection("COPERNICUS/DEM/GLO30").mosaic().setDefaultProjection('EPSG:4326', None, 30)
    image= glo.select(4).gt(0)
    bands = image.bandNames().getInfo()
    byteImage = image.cast({band: 'byte' for band in bands})

    feature = ee.Feature(roi)
    features= ee.FeatureCollection([feature])

    download_image_geemap(byteImage,
                                    features,
                                    out_path= out_path,
                                    crs= 'EPSG:4326',
                                    overwrite=True,
                                    scale= 30)

## Method Specific Functions

import numpy as np
import rasterio
import subprocess
import os
import tempfile
from scipy.ndimage import binary_dilation, binary_erosion

def create_slope_stack(dem_path, forest_height_path, k_values, dem_smooth_kernel=None, dem_folder=None, slope_folder=None):
    # Load the DEM and forest height data
    with rasterio.open(dem_path) as src_dem:
        dem = src_dem.read(1)
        dem_profile = src_dem.profile
        
    with rasterio.open(forest_height_path) as src_fh:
        forest_height = src_fh.read(1)
        
    forest_height_smooth = custom_mean_filter(forest_height, kernel_size=5)
    
    # Initialize an empty list to hold the slope data arrays
    slope_stack = []
    
    # Ensure the save_folder exists if it is provided
    if dem_folder and not os.path.exists(dem_folder):
        os.makedirs(dem_folder)
        
    # Ensure the save_folder exists if it is provided
    if slope_folder and not os.path.exists(slope_folder):
        os.makedirs(slope_folder)
    
    for k in k_values:
        # Correct the DEM and smooth borders based on the forest height and k value
        corrected_dem = dem - (k * forest_height_smooth)
        
        # Further smooth the corrected DEM if smooth_kernel is specified
        if dem_smooth_kernel is not None and dem_smooth_kernel >= 3:
            corrected_dem = custom_mean_filter(corrected_dem, kernel_size=dem_smooth_kernel)
        
        # Determine the path for the smoothed corrected DEM
        if dem_folder:
            corrected_dem_path = os.path.join(dem_folder, f'mean_corrected_dem_0{int(k*100)}.tif')
        else:
            corrected_dem_path = tempfile.mktemp(suffix='.tif')
        
        # Write the smoothed corrected DEM
        with rasterio.open(corrected_dem_path, 'w', **dem_profile) as dst:
            dst.write(corrected_dem, 1)
        
        # Determine the path for the slope file
        if slope_folder:
            slope_file_path = os.path.join(slope_folder, f'mean_slope_0{int(k*100)}.tif')
        else:
            slope_file_path = tempfile.mktemp(suffix='.tif')
        
        # Calculate slope using gdaldem slope
        subprocess.run(['gdaldem', 'slope', corrected_dem_path, slope_file_path, '-of', 'GTiff', '-b', '1', '-s', '1.0', '-p'], check=True)
        
        # Load and append the slope data to the slope stack
        with rasterio.open(slope_file_path) as src_slope:
            slope_data = src_slope.read(1)
            slope_stack.append(slope_data)
        
        # Clean up temporary files if no save_folder is specified
        if not dem_folder:
            os.remove(corrected_dem_path)
            
        if not slope_folder:
            os.remove(slope_file_path)
    
    # Return the stack of slope data for all k_values
    return np.stack(slope_stack)

import numpy as np
from scipy.ndimage import binary_erosion, generate_binary_structure

def find_border(mask, K):
    """
    Identifies the border pixels of a binary mask and sets the edge pixels to zero based on K.

    Parameters:
    - mask: 2D numpy array representing the binary mask.
    - K: The size of the neighborhood for finding local maxima (must be an odd integer).

    Returns:
    - A modified mask with border pixels identified and edges set to zero according to K.
    """
    if K % 2 == 0:
        raise ValueError("K must be an odd integer to have a central pixel.")
    
    border = mask.astype(int) - binary_erosion(mask).astype(int)
    
    # Calculate the padding size needed based on K
    pad_size = K // 2
    
    # Set the border of the array to zero according to K
    if pad_size > 0:
        border[:pad_size, :] = 0
        border[-pad_size:, :] = 0
        border[:, :pad_size] = 0
        border[:, -pad_size:] = 0
    
    return border


def find_pixel_neighborhood_maxima(data, row, col, K):
    
    from scipy.ndimage import maximum_position
    
    if K % 2 == 0:
        raise ValueError("K must be an odd integer to have a central pixel.")
    
    half_k = K // 2
    # Calculate the bounds of the neighborhood, ensuring they are within the data's dimensions
    row_start = max(0, row - half_k)
    row_end = min(data.shape[0], row + half_k + 1)
    col_start = max(0, col - half_k)
    col_end = min(data.shape[1], col + half_k + 1)
    
    neighborhood = data[row_start:row_end, col_start:col_end]
    # Use maximum_position to find the position of the maximum within the neighborhood
    max_pos_local = maximum_position(neighborhood)
    
    # Convert local maximum position back to global position within the original array
    # No subtraction by 1 needed; adjust according to the start of the neighborhood
    max_row_abs = row_start + max_pos_local[0]
    max_col_abs = col_start + max_pos_local[1]

    return (max_row_abs, max_col_abs)

def find_borders_neighborhood_maxima(data, forest_mask, K):
    border = find_border(forest_mask, K)
    row_indices, col_indices = np.where(border == 1)
    maxima_positions = [find_pixel_neighborhood_maxima(data, row, col, K) for row, col in zip(row_indices, col_indices)]
    row_maxima, col_maxima = zip(*maxima_positions)
    maxima_map = np.zeros(dem.shape, dtype=int)
    maxima_map[row_maxima, col_maxima] = 1
    return maxima_map


def find_best_k_pixel(slope_stack, k_values, mask):
    """
    For each pixel indicated by the mask, find the k value that results in the minimum slope.

    Parameters:
    - slope_stack: 3D numpy array, where each layer corresponds to slope values for a specific k value.
    - k_values: List of k values corresponding to each layer in the slope_stack.
    - mask: 2D binary numpy array indicating the pixels to be evaluated.

    Returns:
    - A 2D numpy array with the same shape as mask, where each pixel contains the k value 
      resulting in the minimum slope for that pixel. Pixels outside the mask are set to np.nan.
    """
    # Find the index of the minimum slope for each pixel
    min_slope_indices = np.argmin(slope_stack, axis=0)
    
    # Map indices to k_values
    best_k_map = np.array(k_values)[min_slope_indices]
    
    # Apply the mask, setting pixels outside the mask to np.nan
    best_k_map = np.where(mask ==1, best_k_map, np.nan)

    return best_k_map

from scipy.ndimage import uniform_filter

def find_best_k_pixel_avg(slope_stack, k_values, mask, neighborhood_size):
    """
    For each pixel indicated by the mask, find the k value that results in the minimum average slope
    in the neighborhood around the pixel.

    Parameters:
    - slope_stack: 3D numpy array, where each layer corresponds to slope values for a specific k value.
    - k_values: List of k values corresponding to each layer in the slope_stack.
    - mask: 2D binary numpy array indicating the pixels to be evaluated.
    - neighborhood_size: The size of the neighborhood around each pixel (must be an odd integer).

    Returns:
    - A 2D numpy array with the same shape as mask, where each pixel contains the k value 
      resulting in the minimum average slope for that pixel. Pixels outside the mask are set to np.nan.
    """
    if neighborhood_size % 2 == 0:
        raise ValueError("neighborhood_size must be an odd integer to have a central pixel.")

    # Apply a mean filter to each layer in the slope stack
    smoothed_slope_stack = np.array([uniform_filter(slope_layer, size=neighborhood_size) for slope_layer in slope_stack])

    # Find the index of the minimum slope for each pixel
    min_slope_indices = np.argmin(smoothed_slope_stack, axis=0)
    
    # Map indices to k_values
    best_k_map = np.array(k_values)[min_slope_indices]
    
    # Apply the mask, setting pixels outside the mask to np.nan
    best_k_map = np.where(mask == 1, best_k_map, np.nan)

    return best_k_map

import rasterio
import numpy as np
from scipy.ndimage import binary_dilation

def filter_samples(samples, mask_file1):
    # Read the water bodies masks
    with rasterio.open(mask_file1) as src1:
        mask1 = src1.read(1)


    # Dilate the merged mask
    dilated_mask = binary_dilation(mask1, iterations=1)

    # Set the sample pixels that are within the dilated mask to no data (0)
    samples[np.where(dilated_mask)] = np.nan

    return samples

import rasterio
import numpy as np
from scipy.spatial import cKDTree, KDTree
import matplotlib.pyplot as plt
from scipy.ndimage import label as label_function
from scipy.ndimage import grey_dilation, convolve

import numpy as np
import rasterio
from scipy.ndimage import label, grey_dilation, convolve

import numpy as np
from scipy.spatial import KDTree

def match_labels(labeled_mask1, mask2):
    # Find the coordinates of pixels in mask 1 that are part of a component
    mask1_coords = np.column_stack(np.nonzero(labeled_mask1))
    
    # Create a KDTree for efficient nearest neighbor search
    tree = KDTree(mask1_coords)
    
    # Initialize the result mask with zeros, ensuring it can hold the label values (integers)
    result_labeled_mask1 = labeled_mask1.copy()
    result_mask2 = np.zeros_like(mask2, dtype=labeled_mask1.dtype)
    
    # Iterate over each pixel in mask 2 where mask2 is nonzero (indicating areas of interest)
    nonzero_mask2_coords = np.column_stack(np.nonzero(mask2))
    
    for coord in nonzero_mask2_coords:
        # Find the nearest pixel in mask 1 that is part of a component
        _, idx = tree.query(coord)
        
        # Copy the label of the nearest component from mask 1 to mask 2
        result_mask2[coord[0], coord[1]] = labeled_mask1[mask1_coords[idx][0], mask1_coords[idx][1]]
    
    # Identify labels in mask 1 that have no samples in mask 2
    unique_labels_in_mask1 = set(np.unique(labeled_mask1))
    unique_labels_in_mask2 = set(np.unique(result_mask2))
    labels_with_no_samples = unique_labels_in_mask1 - unique_labels_in_mask2
    
    # Create a KDTree for only valid labels in result_mask2
    valid_coords = np.column_stack(np.nonzero(result_mask2))
    valid_tree = KDTree(valid_coords)
    
    for label in labels_with_no_samples:
        # Find the coordinates of pixels with this label in mask 1
        coords_with_no_samples = np.column_stack(np.nonzero(labeled_mask1 == label))
        
        for coord in coords_with_no_samples:
            # Find the nearest label that has samples in result_mask2
            _, idx = valid_tree.query(coord)
            nearest_label = result_mask2[valid_coords[idx][0], valid_coords[idx][1]]
            
            # Update the label in labeled_mask1
            result_labeled_mask1[coord[0], coord[1]] = nearest_label
    
    return result_labeled_mask1, result_mask2


def average_interpolation_by_label(labeled_mask, samples_data):
    binary_samples_data = np.where(~np.isnan(samples_data), 1, 0)
    labeled_mask, matched_samples_labels = match_labels(labeled_mask, binary_samples_data)
    valid_samples = ~np.isnan(samples_data)
    sample_values = samples_data[valid_samples]
    sample_labels = matched_samples_labels[valid_samples]
    unique_labels = np.unique(labeled_mask)
    final_output_data = np.zeros(labeled_mask.shape)

    for label in unique_labels:
        if label == 0:
            continue

        current_label_mask = (labeled_mask == label).astype(int)
        current_sample_values = sample_values[sample_labels == label]
        average_value = np.mean(current_sample_values)
        final_output_data[current_label_mask == 1] = average_value

    return final_output_data

import rasterio
import numpy as np
from scipy.spatial import cKDTree

def create_custom_year_forest_mask(lossyear_path, fh2019_path, input_year, output_path, k):
    # Open the lossyear dataset
    with rasterio.open(lossyear_path) as src_lossyear:
        lossyear = src_lossyear.read(1)
        profile = src_lossyear.profile
    
    # Open the forest height dataset
    with rasterio.open(fh2019_path) as src_fh2019:
        fh2019 = src_fh2019.read(1)
    
    # Mask out the invalid lossyear values (lossyear > 19)
    valid_lossyear_mask = lossyear <= 20
    lossyear = np.where(valid_lossyear_mask, lossyear, 0)
    
    # Find pixels where the forest height is zero and lossyear is greater than or equal to the input year
    valid_lossyear_gte_input_year = lossyear >= input_year
    
    # Coordinates of pixels to be interpolated
    interp_coords = np.column_stack(np.where(valid_lossyear_gte_input_year))
    
    # Coordinates and values of non-zero height pixels where lossyear is 0
    valid_non_zero_mask = (fh2019 > 0) & (lossyear == 0)
    non_zero_height_coords = np.column_stack(np.where(valid_non_zero_mask))
    non_zero_height_values = fh2019[valid_non_zero_mask]
    
    # Build a KDTree for the non-zero height pixels with lossyear == 0
    tree = cKDTree(non_zero_height_coords)
    
    # Interpolate using the average of the k nearest neighbors
    interpolated_values = []
    for coord in interp_coords:
        dists, idxs = tree.query(coord, k=k)
        neighbor_values = non_zero_height_values[idxs]
        interpolated_value = np.mean(neighbor_values)
        interpolated_values.append(interpolated_value)
    
    # Replace the zero forest height pixels with the interpolated values
    interpolated_forest_height = fh2019.copy()
    for (x, y), value in zip(interp_coords, interpolated_values):
        interpolated_forest_height[x, y] = value
    
    # Save the interpolated forest height to a new GeoTIFF
    profile.update(dtype=rasterio.float32, count=1)
    
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(interpolated_forest_height.astype(rasterio.float32), 1)

def calculate_slope(dem_path, slope_path):
    gdal.DEMProcessing(slope_path, dem_path, 'slope', format='GTiff', band=1, scale=1.0, slopeFormat='percent')


def apply_bilateral_filter_on_changed_pixels(original_dem_path, treated_dem_path, output_path, sigma_dist, sigma_int):
    # Read the original DEM
    with rasterio.open(original_dem_path) as src:
        original_dem = src.read(1)
        profile = src.profile

    # Read the treated DEM
    with rasterio.open(treated_dem_path) as src:
        treated_dem = src.read(1)

    # Identify the changed pixels
    changed_pixels_mask = (original_dem != treated_dem)

    # Create a temporary file for the filtered DEM
    temp_filtered_dem_path = tempfile.mktemp(suffix='.tif')
    wbt.bilateral_filter(
        i=treated_dem_path,
        output=temp_filtered_dem_path,
        sigma_dist=sigma_dist,
        sigma_int=sigma_int
    )

    # Read the filtered DEM
    with rasterio.open(temp_filtered_dem_path) as src:
        bilateral_filtered_dem = src.read(1)

    # Apply the bilateral filter only on the changed pixels
    result_dem = np.where(changed_pixels_mask, bilateral_filtered_dem, original_dem)

    # Write the result to a new GeoTIFF
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(result_dem, 1)

    # Clean up temporary file
    os.remove(temp_filtered_dem_path)

def copy_water_body_elevations(original_dem_path, treated_dem_path, mask_path, output_path):
    # Read the original DEM
    with rasterio.open(original_dem_path) as src:
        original_dem = src.read(1)
        profile = src.profile

    # Read the treated DEM
    with rasterio.open(treated_dem_path) as src:
        treated_dem = src.read(1)

    # Read the lakes mask
    with rasterio.open(mask_path) as src:
        mask = src.read(1)

    # Copy elevations of water bodies from the original DEM to the treated DEM
    result_dem = np.where(mask, original_dem, treated_dem)

    # Write the result to a new GeoTIFF
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(result_dem, 1)

# DOWNLOAD DATA FROM GEE

ee.Authenticate() 
ee.Initialize()

### Create output dirs 
if not os.path.exists(base_path):
    os.mkdir(base_path)

if not os.path.exists(aux_path):
    os.mkdir(aux_path)

### Method parameters
s = 3 # neighborhood used for slope analysis
k_values = np.linspace(0, 1, 21) # pre determined set of k values

### Get ROI from extent
gdf = gpd.read_file(shapefile_path).to_crs(epsg=4326)
roi = geemap.gdf_to_ee(gdf).first().geometry()
Map = geemap.Map()

### Download forest height and forest loss
t = 4 # number of tiles to divide the area
fh2020_path =  f'{aux_path}/fh2020.tif'
aligned_fh2020_path =  f'{aux_path}/fh2020_align.tif'
fh2020 = ee.Image("projects/glad/GLCLU2020/Forest_height_2020").unmask(0)
ee_export_image(fh2020, roi, fh2020_path, nrows=t, ncols=t)
align_raster(fh2020_path, raw_dem_path, resampling='near', dtype='Byte', output_path=aligned_fh2020_path)

lossyear_path = f'{aux_path}/lossyear.tif'
lossyear_align_path = f'{aux_path}/lossyear_align.tif'
lossyear = ee.Image('UMD/hansen/global_forest_change_2023_v1_11').select(3)
ee_export_image(lossyear, roi, lossyear_path, nrows=t, ncols=t)
align_raster(lossyear_path, raw_dem_path, resampling='near', dtype='Float32', output_path=lossyear_align_path)

### optional water mask download
# water_mask_path = f'{aux_path}/glo_wbm_raw.tif'
# download_glo_rivers_and_lakes(roi, water_mask_path)
# aligned_water_mask_path = f"{aux_path}/glo-30_wbm.tif"
# align_raster(water_mask_path, raw_dem_path, resampling='near', dtype='Byte', output_path=aligned_water_mask_path)

# DEM CORRECTION

corrected_dems = []
dem_slopes = []
average_slopes = []

for input_year in range(10, 16):
    fh2019_path = f'{aux_path}/fh2020_align.tif'
    lossyear_path = f'{aux_path}/lossyear_align.tif'
    adj_forest_height_path = f'{aux_path}/adjusted_fh2020_align_{input_year}.tif'
    create_custom_year_forest_mask(lossyear_path, fh2019_path, input_year, adj_forest_height_path, k=128)

    #raw_dem_path =  f'{aux_path}/glo-30.tif'
    # aligned_ch2024_path =  f'{aux_path}/fh2020_align_clip.tif'
    aligned_ch2024_path =  adj_forest_height_path
    
    with rasterio.open(raw_dem_path) as src:
        dem = src.read(1)
    
    with rasterio.open(aligned_ch2024_path) as src:
        forest_height = src.read(1)
        profile = src.profile
    
    
    ### Smoothing
    forest_height_smooth = custom_mean_filter(forest_height, kernel_size=5)
    forest_mask = (forest_height > 0).astype(int)
    forest_mask_smooth = (forest_height_smooth > 0).astype(int)
    
    ### Create slope stack
    k_values = np.linspace(0, 1, 21)
    slope_stack = create_slope_stack(raw_dem_path, 
                                     aligned_ch2024_path, 
                                     k_values)
    
    ### Calculate corrected DEM
    maxima = find_borders_neighborhood_maxima(slope_stack[0,:,:], forest_mask, 3)
    best_k_maxima = find_best_k_pixel_avg(slope_stack, k_values, maxima, 3)
    # best_k_maxima = find_best_k_pixel(slope_stack, k_values, maxima)
    best_k_maxima = np.where(best_k_maxima < 0.05, np.nan, best_k_maxima)
    best_k_maxima_filtered = filter_samples(best_k_maxima, water_mask_path)

    
    samples_data = best_k_maxima_filtered
    data = forest_height

    ### Get labeled mask (choose method)
    labeled_mask = grey_dilation(label_function(data)[0], size=5)

    ### Binarize the samples data
    binary_samples_data = np.where(~np.isnan(samples_data), 1, 0)
    
    ### Apply the matching function to get the labeled masks
    labeled_mask, matched_samples_labels = match_labels(labeled_mask, binary_samples_data)

    ### AVERAGE
    interpolated_best_k = average_interpolation_by_label(labeled_mask, samples_data)
    
    ### Correct DEM
    corrected_dem = dem - (interpolated_best_k * forest_height_smooth)
    
    corrected_dem_path = f'{aux_path}/glo_avg_fh2020_year_20{input_year}.tif'
    profile.update(dtype=rasterio.float32)
    profile.update(nodata=-9999)
    with rasterio.open(corrected_dem_path , 'w', **profile) as dst:
        dst.write(corrected_dem, 1)

    corrected_dems.append(corrected_dem_path)

    slope_path = f'{aux_path}/corrected_dem_20{input_year}_slope.tif'
    calculate_slope(corrected_dem_path, slope_path)
    dem_slopes.append(slope_path)

    with rasterio.open(slope_path) as src:
        slope_data = src.read(1)
        average_slope = np.nanmean(slope_data)
        average_slopes.append(average_slope)

### Find the corrected DEM with the lowest average slope
min_slope_index = np.argmin(average_slopes)
best_corrected_dem_path = corrected_dems[min_slope_index]

### Copy the best corrected DEM to the final path
final_corrected_dem_path = f'{aux_path}/best_corrected_dem_unfinished.tif'
with rasterio.open(best_corrected_dem_path) as src:
    dem_data = src.read(1)
    profile = src.profile
    profile.update(dtype=rasterio.float32, nodata=-9999)
    with rasterio.open(final_corrected_dem_path, 'w', **profile) as dst:
        dst.write(dem_data, 1)

print(f"The best corrected DEM for this area is from the year: 20{10 + min_slope_index}")

import whitebox
### Bilateral smoothing with whitebox tools
wbt = whitebox.WhiteboxTools()
wbt.set_verbose_mode(False)


smooth_dem_path = f'{aux_path}/best_corrected_dem_smooth.tif'
apply_bilateral_filter_on_changed_pixels(
    original_dem_path= raw_dem_path,
    treated_dem_path= final_corrected_dem_path,
    output_path=smooth_dem_path,
    sigma_dist=3,
    sigma_int=5.0)

### Recover water bodies elevation

copy_water_body_elevations(raw_dem_path,
                           smooth_dem_path,
                           water_mask_path,
                           output_path_final)