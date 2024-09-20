import geopandas as gpd
import pandas as pd
from pathlib import Path


URLS = {
    'building_footprints':"https://gisapp.wycokck.org/gisdata/shp/bldg_ftpt_py.zip",
    'zoning':"https://gisapp.wycokck.org/gisdata/shp/zoning_py.zip",
    'vacant_parcels':"https://gisapp.wycokck.org/gisdata/shp/vacant_py.zip",
    'parcels':"https://gisapp.wycokck.org/gisdata/shp/parcel_py.zip",
    'impervious_land_cover':"https://gisapp.wycokck.org/gisdata/shp/impervious_landcover_py.zip",
    'land_bank_parcels':'https://gisapp.wycokck.org/gisdata/shp/landbank_py.zip',
    
}

if __name__ == "__main__":
    
    community_cutout = gpd.read_file(snakemake.input.community)
    
    for name, url in URLS.items():
        gdf = gpd.read_file(url)
        gdf = gdf.to_crs(snakemake.config['geographic_crs'])
        
        gdf_cutout = gdf.sjoin(community_cutout, predicate='within')
        print(f"data/spatial_data/armourdale/{name}.gpkg")
        
        gdf.to_file(f"data/spatial_data/armourdale/{name}.gpkg", driver="GPKG")
        