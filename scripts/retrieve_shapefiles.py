import geopandas as gpd
from pyogrio.errors import DataSourceError
from tqdm import tqdm


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
    
    pbar = tqdm(URLS.items())
    for name, url in pbar:
        pbar.set_description(f"{name}")
        try:
            gdf = gpd.read_file(url)
            gdf = gdf.to_crs(snakemake.config['geographic_crs'])
            
            gdf_cutout = gdf.sjoin(community_cutout, predicate='within')
            
            gdf.to_file(f"data/spatial_data/armourdale/{name}.gpkg", driver="GPKG")
        except DataSourceError:
            print(f"Failed to download {name} from {url}")