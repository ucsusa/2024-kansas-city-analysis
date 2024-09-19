import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
from census import Census
from us import states

if __name__ == "__main__":
    kc_wards_url = "https://maps.wycokck.org/gisdata/shp/ward_prec_py.zip"
    kck_wards = gpd.read_file(kc_wards_url).to_crs(
        epsg=int(snakemake.config['geographic_crs']))

    armourdale_ward = '06'
    armourdale = kck_wards[kck_wards['WARD'] == armourdale_ward].dissolve(
        "CITY").reset_index(drop=False)

    armourdale = armourdale[['CITY', 'geometry', 'WARD']]

    armourdale.to_file(snakemake.output.armourdale, driver="GPKG")
