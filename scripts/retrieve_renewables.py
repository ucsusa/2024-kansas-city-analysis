import numpy as np
import sys
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from unyt import m, s, MW, W, kg, g

sys.path.append("functions")

from nrel_data_api import parameters, make_csv_url

# model_years = np.array(snakemake.config['model_years']).astype('int')
model_years = np.array([2018]).astype('int')

def handle_datetime(dataframe):
    """
    Combines time columns into a single timestamp column.
    Expects columns ['year','month','day','hour'].

    Parameters
    ----------
    dataframe : :class:`pd.DataFrame`
        A pandas dataframe.
    """
    frame = dataframe.copy()
    timestamps = pd.DatetimeIndex([])
    for year in model_years:
        period = pd.date_range(start=f"{year}-01-01",
                                        freq=f"1h",
                                        periods=8760)
        timestamps = timestamps.append(period)
    
    frame.index = timestamps
    try:
        frame.set_index(timestamps,inplace=True)
        frame.drop(columns=['Year','Month','Day','Hour','Minute'],inplace=True)
    except:
        raise ValueError
    
    return frame
    

def retrieve_solar_timeseries(region, save_years=True):
    """
    Retrieves data from NREL's national solar radiation database (NSRDB).

    Parameters
    ----------
    region : :class:`gpd.GeoDataFrame`
        A geopandas dataframe containing modeled bus regions.
    """
   
    parameters['attr_list'] = ['ghi']
    outer_pbar = tqdm(snakemake.config['solar_years'], position=0, leave=True)
    all_frames = []
    for year in outer_pbar:
        outer_pbar.set_description(f"Processing {year}")
        parameters['year'] = int(year)
        frames = []
        inner_pbar = tqdm(region[['name','x','y']].values, position=1, leave=True)
        for n, i, j in inner_pbar:
            inner_pbar.set_description(f"Processing {n}")
            parameters['lon'] = i
            parameters['lat'] = j
            URL = make_csv_url(parameters=parameters, 
                            kind='solar')
            df = pd.read_csv(URL, skiprows=2)
            df['date'] = pd.to_datetime(df[['Year','Month','Day','Hour','Minute']])
            df = df.drop(columns=['Year','Month','Day','Hour','Minute']).set_index('date')
            frames.append(df)
    
        solar_df = pd.concat(frames, axis=1)
        if save_years:
            solar_df.to_csv(f"data/timeseries/solar_{year}.csv")    
        all_frames.append(solar_df)

    full_df = pd.concat(all_frames, axis=0)
    # full_df = handle_datetime(full_df)
    
    return full_df


def process_solar_timeseries(df, normalize=True):
    """
    Converts solar radiation timeseries to a
    hypothetical power production.

    Parameters
    ----------
    df : _type_
        _description_
    normalize : bool, optional
        Whether the data should be normalized. Default is true.
    """
    frame = df.copy()
    if normalize:
            frame = frame.divide(frame.max(axis=0), axis=1)
        
    return frame 

if __name__ == "__main__":
    
    regions = gpd.read_file(snakemake.input.supply_regions)
    
    # solar data
    df = retrieve_solar_timeseries(regions)
    df = process_solar_timeseries(df)
    df.to_csv(snakemake.output.solar)