import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
from census import Census
from us import states
import os
import sys
import yaml

sys.path.append("utils/")
from api_functions import get_tiger_files, get_county_fips

from api_functions import get_tiger_files, get_county_fips


column_names = {
    "B01003_001E":"total_population",
    "B25024_002E":"single-family_detached",
    "B25024_003E":"single-family_attached",
    "B25024_004E":"2 units",
    "B25024_005E":"3-4_units",
    "B25024_006E":"5-9_units",
    "B25024_007E":"10-19_units",
    "B25024_008E":"20-49_units",
    "B25024_009E":"50plus_units",
    "B25024_010E":"mobile_home",
}


if __name__ == "__main__":
    # gather config options
    state_name = snakemake.config['state']
    state = states.lookup(state_name)
    county = snakemake.config['county']
    county_fips = get_county_fips(state_name, county)
    census_year = int(snakemake.config['census_year'])
    
    
    # get census data
    api_key = os.environ.get('CENSUS_API_KEY')
    c = Census(api_key)
    county_census = c.acs5.state_county_blockgroup(fields=tuple(column_names.keys()),
                               state_fips=state.fips,
                               county_fips=str(county_fips),
                               blockgroup="*",
                               year=census_year)
    county_df = pd.DataFrame(county_census)
    county_df.rename(columns=column_names, inplace=True)
    
    county_df['GEOID'] = county_df['state'] + county_df['county'] + county_df['tract'] + county_df['block group']
    county_df.drop(columns=['state','county','tract','block group'], inplace=True)
    
    
    # get the map of state level block groups
    state_map = get_tiger_files(year=census_year,
                         state_abbr = state.abbr,
                         feature='blockgroup')
    state_map = state_map.to_crs(epsg=int(snakemake.config['geographic_crs']))
    
    county_bg = state_map[state_map['COUNTYFP'] == str(county_fips)]
    county_bg = county_bg.drop(columns = ['STATEFP','COUNTYFP', 'TRACTCE',
                                'BLKGRPCE', 'NAMELSAD', 'MTFCC',
                                'FUNCSTAT', 'ALAND','AWATER',
                                'INTPTLAT', 'INTPTLON'])
    
    county_merge = county_bg.merge(county_df, on='GEOID')
    
    
    # combine structure types by unit; harmonize with NREL resstock
    multi_family = ['2 units','3-4_units']
    many_family = ['5-9_units', '10-19_units', '20-49_units','50plus_units']
    
    county_merge['multi-family_with_2_-_4_units'] = county_merge[multi_family].sum(axis=1)
    county_merge['multi-family_with_5plus_units'] = county_merge[many_family].sum(axis=1)
    county_merge = county_merge.drop(columns=multi_family+many_family)
    
    county_merge.to_file(snakemake.output.census_data, driver="GPKG")
    
    state_map.to_file(snakemake.output.state_blockgroups, driver="GPKG")
    county_bg.to_file(snakemake.output.county_blockgroups, driver="GPKG")


    
    