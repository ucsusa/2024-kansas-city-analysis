import requests
import json
import pandas as pd
import numpy as np

if __name__ == "__main__":
    
        
    # gather config options
    state_name = snakemake.config['state']
    state = states.lookup(state_name)
    county = snakemake.config['county']
    building_opts = snakemake.config['building_data_options']
    
    lut = pd.read_csv(snakemake.input.spatial_lut)
    
    
    # get the PUMA ID
    # this method is unstable, since some counties might contain multiple PUMAs
    county_and_puma = lut[((lut['state_abbreviation']==state.abbr)\
                          & (lut['resstock_county_id'] == f"{state.abbr}, {county.capitalize()} County"))]['nhgis_puma_gisjoin'].unique()[0]