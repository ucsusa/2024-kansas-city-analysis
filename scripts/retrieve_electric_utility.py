import pandas as pd
import geopandas as gpd
import sys

sys.path.append("utils")
from api_functions import get_retail_service_area

if __name__ == "__main__":
    state_name = snakemake.config['state']

    cutout = gpd.read_file(snakemake.input.cutout)

    url = get_retail_service_area(state_name=state_name)

    service_gdf = gpd.read_file(url)

    service_gdf = service_gdf.sjoin(cutout)

    service_gdf.to_file(snakemake.output.utility)
