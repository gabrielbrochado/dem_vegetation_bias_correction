import pandas as pd
import geopandas as gpd

#Paths to the files 
base_path = 'C:/Users/gabro/Downloads' # <- where suplementary materials was unziped
file_path = base_path + '/Suplementary_Materials_Brochado_Renno_2024/selected_stacked_dataframe_r_class_p.csv'

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(file_path)

for area in [1,2,3,4]:
    for radius in [1000,2000,3000]:

        # Filter by area and radius
        filtered_df = df[(df['area'] == area) & (df['radius'] == radius)].copy()
        
        # Get selected point ids
        point_ids = filtered_df['point_id'].values

        # Load polygons 
        glo_gpkg_path = base_path + f'/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/b{radius}_glo_polygons_r_class_p.gpkg'
        fabdem_gpkg_path = base_path + f'/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/b{radius}_fabdem_polygons_r_class_p.gpkg'
        corrected_gpkg_path = base_path + f'/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/b{radius}_corrected_polygons_r_class_p.gpkg'

        glo_gpkg_gdf = gpd.read_file(glo_gpkg_path)
        fabdem_gpkg_gdf = gpd.read_file(fabdem_gpkg_path)
        corrected_gpkg_gdf = gpd.read_file(corrected_gpkg_path)

        # Filter gdfs by selected point ids
        selected_glo_gpkg_gdf = glo_gpkg_gdf[glo_gpkg_gdf['point_id'].isin(point_ids)].copy()
        selected_fabdem_gpkg_gdf = fabdem_gpkg_gdf[fabdem_gpkg_gdf['point_id'].isin(point_ids)].copy()
        selected_corrected_gpkg_gdf = corrected_gpkg_gdf[corrected_gpkg_gdf['point_id'].isin(point_ids)].copy()

        # Write filtered gdfs
        selected_glo_gpkg_gdf.to_file(base_path + f'/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/selected_b{radius}_glo_polygons_r_class_p.gpkg', driver="GPKG")
        selected_fabdem_gpkg_gdf.to_file(base_path + f'/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/selected_b{radius}_fabdem_polygons_r_class_p.gpkg', driver="GPKG")
        selected_corrected_gpkg_gdf.to_file(base_path + f'/Suplementary_Materials_Brochado_Renno_2024/Area_{area}/selected_b{radius}_corrected_polygons_r_class_p.gpkg', driver="GPKG")