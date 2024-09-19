import pandas as pd
import geopandas as gpd
from us import states


def create_resstock_url(state_abbr,
                        puma_id,
                        building_type,
                        year=2021,
                        product='resstock',
                        weather_version='tmy3',
                        release=1,
                        ):

    BASE_URL = (
        "https://oedi-data-lake.s3.amazonaws.com/nrel-pds-building-stock"
        "/end-use-load-profiles-for-us-building-stock")

    data_route = (f"/{year}"
                  f"/{product}_{weather_version}_release_{release}"
                  "/timeseries_aggregates/by_puma"
                  f"/state={state_abbr}/")

    file = f"{puma_id.lower()}-{building_type}.csv"

    return BASE_URL + data_route + file


if __name__ == "__main__":

    columns = ['timestamp', 'out.electricity.total.energy_consumption']

    # gather config options
    state_name = snakemake.config['state']
    state = states.lookup(state_name)
    county = snakemake.config['county']
    building_opts = snakemake.config['building_data_options']
    sector_buildings = building_opts['building_types']

    # load spatial lut
    lut = pd.read_csv(snakemake.input.spatial_lut)

    # get the PUMA ID
    # this method is unstable, since some counties might contain multiple PUMAs
    county_and_puma = lut[((lut['state_abbreviation'] == state.abbr) & (
        lut['resstock_county_id'] == f"{state.abbr}, {county.capitalize()} County"))]['nhgis_puma_gisjoin'].unique()[0]

    # puma_id = county_and_puma.split(',')[-1].replace(' ','')

    # for sector in list(sectors_buildings.keys()):
    for sector in ['residential']:
        building_types = sector_buildings[sector]
        sector_frames = []
        for bldg_type in building_types:
            # for bldg_type in ["single-family_attached"]:
            bldg_url = create_resstock_url(state_abbr=state.abbr,
                                           puma_id=county_and_puma,
                                           building_type=bldg_type)
            bldg_df = pd.read_csv(bldg_url,
                                  parse_dates=True,
                                  index_col='timestamp',
                                  usecols=columns)

            bldg_df.rename(
                columns={
                    'out.electricity.total.energy_consumption': bldg_type},
                inplace=True)

            sector_frames.append(bldg_df)

        sector_ts = pd.concat(sector_frames, axis=1)
        sector_ts.to_csv(f"data/timeseries/{sector}_load.csv")
