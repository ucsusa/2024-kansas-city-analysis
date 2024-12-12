import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from us import states

from pathlib import Path
import zipfile
import io
import requests


if __name__ == "__main__":
    state_name = snakemake.config['state']
    state = states.lookup(state_name)
    state_abbr = state.abbr

    # download the LEAD data
    url = f"https://data.openei.org/files/573/{state_abbr}-2018-LEAD-data.zip"
    r = requests.get(url)

    z = zipfile.ZipFile(io.BytesIO(r.content))
    z.extractall(f"data/spatial_data/{state_abbr}-2018-LEAD-data/")


    # cutout data for local community
    dataset = snakemake.output.lead_data
    
    lead_df = pd.read_csv(dataset)
    
    county_bg = gpd.read_file(snakemake.input.county_blockgroups)
    community = gpd.read_file(snakemake.input.community)
    
    bg = int(county_bg.sjoin(community, predicate='within').GEOID.values[0][:11])
    
    community_lead = lead_df[lead_df['FIP']==bg]
    
    community_lead.to_csv(snakemake.output.lead_community)