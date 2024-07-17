configfile: "config.yml"

from pathlib import Path
env_file = Path("./.env").resolve()
from dotenv import load_dotenv
load_dotenv(str(env_file))


rule retrieve_spatial_lut:
    output: 
        spatial_lut = "data/spatial_data/spatial_lut.csv"
    script: "scripts/retrieve_lut.py"

rule retrieve_census_data:
    output:
        census_data = "data/spatial_data/county_census_data.gpkg"
    script: "scripts/retrieve_census_data.py"

# a bespoke step to make this analysis specific to armourdale
rule retrieve_armourdale_shape:
    output: 
        armourdale = "data/spatial_data/armourdale_shape.gpkg"
    script: "scripts/retrieve_armourdale.py"

rule calculate_res_structures:
    input: 
        census_data = "data/spatial_data/county_census_data.gpkg",
        armourdale = "data/spatial_data/armourdale_shape.gpkg"
    output: 
        res_structures = "data/residential_buildings.csv"
    script: "scripts/calculate_res_structures.py"

rule retrieve_res_load:
    input:
        spatial_lut = "data/spatial_data/spatial_lut.csv",
        res_structures = "data/residential_buildings.csv"
    output: 
        sfa = "data/timeseries/single-family_attached_load.csv"
    script: "scripts/retrieve_res_load.py"
    