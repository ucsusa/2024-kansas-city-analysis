import requests
import json
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import zipfile
import io
import glob
from us import states

_TIGER_URL = "https://www2.census.gov/geo/tiger/"

def get_tiger_files(year, state_abbr, feature='tract'):
    """
    This function retrievs a TIGER shapefile from the United States Census
    website.

    Parameters
    ----------
    year : int
        The shapefile year of interest.
    state_abbr : str
        The abbreviation for the state of interest.
    feature : str, optional
        Indicates which TIGER file data feature to extract, by default 'tract'.
    """
    
    
    try:
        state = states.lookup(state_abbr)
        assert state, f"{state_abbr} is not a state in the U.S."
    except AssertionError as error:
        raise error

    _FEATURE_FILE = {'tract':f'TRACT/tl_{year}_{state.fips}_tract.zip',
                     'blockgroup':f'BG/tl_{year}_{state.fips}_bg.zip',
                     'county':f"COUNTY/tl_{year}_us_county.zip"}
    data_route = f"TIGER{year}/{_FEATURE_FILE[feature]}"
    
    geo_df = gpd.read_file(_TIGER_URL+data_route)
    
    return geo_df