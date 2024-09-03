import pandas as pd

resstock_opts = snakemake.config['building_data_options']

BASE_URL = (f"https://oedi-data-lake.s3.amazonaws.com/nrel-pds-building-stock"
            f"/end-use-load-profiles-for-us-building-stock/")

URL = BASE_URL + (
    f"{resstock_opts['resstock_year']}"
    f"/resstock_{resstock_opts['weather_version']}_release_{resstock_opts['release_version']}"
    f"/geographic_information"
    f"/spatial_tract_lookup_table.csv")


if __name__ == "__main__":

    df = pd.read_csv(URL)
    df.to_csv(snakemake.output.spatial_lut)
