import pandas as pd
import geopandas as gpd
import datetime as dt


if __name__ == "__main__":
    
    if snakemake.config['usrdb_start_date'].lower() == 'today':
        start_date = pd.to_datetime(dt.date.today())
    else:
        start_date = pd.to_datetime(snakemake.config['usrdb_start_date'])
        
    future_date = pd.to_datetime("2099-01-01")
    
    URL = "https://apps.openei.org/USURDB/download/usurdb.csv.gz"
    usrdb = pd.read_csv(URL, low_memory=False, parse_dates=True)
    
    # filter by date
    usrdb.loc[:,'enddate'] = pd.to_datetime(usrdb['enddate'].fillna(future_date))
    usrdb.loc[:,'startdate'] = pd.to_datetime(usrdb['startdate'])
    usrdb = usrdb[(usrdb['enddate'] > start_date)]
    
    # get utility info
    utility_service = gpd.read_file(snakemake.input.utility)
    utility_ids = utility_service.loc[:, 'ID'].values.astype(int)

    sectors = [sector.capitalize() 
               for sector 
               in snakemake.config['energy_sectors']]
    
    # filter by utility and sector
    usrdb = usrdb.loc[usrdb['eiaid'].isin(utility_ids)]
    usrdb = usrdb.loc[usrdb['sector'].isin(sectors)]
    
    # filter: is default?
    usrdb = usrdb[usrdb['is_default']].dropna(how='all',axis=1)
    
    usrdb.to_csv(snakemake.output.rates)