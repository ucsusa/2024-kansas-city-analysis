import geopandas as gpd
import pandas as pd

if __name__ == "__main__":
    
    census_data = gpd.read_file(snakemake.input.census_data)
    building_opts = snakemake.config['building_data_options']
    
    
    # specific to armourdale; selecting only block groups within armourdale
    armourdale = gpd.read_file(snakemake.input.armourdale)
    armourdale_bg = census_data.sjoin(armourdale, how='inner', predicate='within')
    
    # replace `armourdale_bg` with `census_data` for a non-armourdale analysis.
    building_data = armourdale_bg[building_opts['building_types']]\
                    .sum()\
                    .to_frame()\
                    .rename(columns={0:'n_units'})

    building_data.to_csv(snakemake.output.res_structures)
    
    
    