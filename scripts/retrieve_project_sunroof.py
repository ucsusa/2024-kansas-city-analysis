import pandas as pd
import geopandas as gpd

URL = "https://storage.googleapis.com/project-sunroof/csv/latest/project-sunroof-census_tract.csv"

if __name__ == "__main__":
    # download project sunroof data
    df = pd.read_csv(URL)
    df.to_csv(snakemake.output.project_sunroof)

    # load the blockgroup shapefile
    print("processing local data")
    census_tract = gpd.read_file(snakemake.input.blockgroups)
    census_tract['region_name'] = (
        census_tract['STATEFP'] +
        census_tract['COUNTYFP'] +
        census_tract['TRACTCE']).astype("int64")
    census_tract = census_tract.dissolve('region_name').reset_index()[
        ['region_name', 'geometry']]

    # merge dataframes
    solar_gdf = census_tract.merge(df, on='region_name')
    solar_gdf.to_file(snakemake.output.local_potential, driver='GPKG')
