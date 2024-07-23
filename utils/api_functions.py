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
_RETAIL_SERVICE_URL = ("https://services1.arcgis.com/Hp6G80Pky0om7QvQ/"
                      "arcgis/rest/services/Retail_Service_Territories/"
                      "FeatureServer/0/query?")
RETAIL_SERVICE_COLUMNS = ["CNTRL_AREA",
                          "PLAN_AREA",
                          "HOLDING_CO",
                          "NET_GEN",
                          "PURCHASED",
                          "RETAIL_MWH",
                          "WSALE_MWH",
                          "TOTAL_MWH",
                          "TRANS_MWH",
                          "CUSTOMERS",
                          "YEAR",
                          "NET_EX",
                          "NAME",
                          "REGULATED",
                          "STATE",
                          "ID",
                          "NAICS_CODE",
                          "NAICS_DESC"]


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





def get_retail_service_area(state_name=None,
                            crs=4326,
                            columns=RETAIL_SERVICE_COLUMNS):
    
    try:
        state = states.lookup(state_name)
        assert state_name, f"{state_name} is not a state in the U.S."
    except AssertionError as error:
        raise error
    
    state_field = f"where=STATE%20%3D%20'{state.abbr}'" if state_name else ""
    crs_field = f"outSR={crs}"
    format_field = f"f=json"