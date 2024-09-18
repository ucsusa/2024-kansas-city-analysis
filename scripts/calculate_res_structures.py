import geopandas as gpd
import pandas as pd

if __name__ == "__main__":
    
    census_data = gpd.read_file(snakemake.input.census_data)
    building_opts = snakemake.config['building_data_options']
    
    
    # specific to armourdale; selecting only block groups within the community
    community = gpd.read_file(snakemake.input.community)
    community_bg = census_data.sjoin(community, how='inner', predicate='within')
    
    # replace `community_bg` with `census_data` for a non-localized analysis.
    building_types = building_opts['building_types']['residential']
    building_data = community_bg[building_types]\
                    .sum()\
                    .to_frame()\
                    .rename(columns={0:'n_units'})

    building_data.to_csv(snakemake.output.res_structures)
    
    
    