configfile: "config.yml"

from us import states
from pathlib import Path
from dotenv import load_dotenv

state = config['state']
state_abbr = states.lookup(state).abbr

env_file = Path("./.env").resolve()
load_dotenv(str(env_file))

rule targets:
    input: 
        armourdale = "data/spatial_data/armourdale_shape.gpkg",
        census_data = "data/spatial_data/county_census_data.gpkg",
        state_blockgroups = f"data/spatial_data/{state.lower()}_blockgroups.gpkg",
        county_blockgroups = f"data/spatial_data/{config['county'].lower()}_blockgroups.gpkg",
        elec_load = "data/timeseries/residential_elec_load.csv",
        heat_load = "data/timeseries/residential_heat_load.csv",
        weather = "data/timeseries/weather_year.csv",
        res_structures = "data/residential_buildings.csv",
        rates = "data/usrdb_rates.csv",
        project_sunroof = f"data/spatial_data/project-sunroof-census_tract.csv",
        utility="data/spatial_data/electric_utility.gpkg",
        lead_data = f"data/spatial_data/{state_abbr}-2018-LEAD-data/{state_abbr} AMI Census Tracts 2018.csv",
        dag = "dag.png"

rule retrieve_spatial_lut:
    output: 
        spatial_lut = "data/spatial_data/spatial_lut.csv"
    script: "scripts/retrieve_lut.py"

rule retrieve_census_data:
    output:
        census_data = "data/spatial_data/county_census_data.gpkg",
        state_blockgroups = f"data/spatial_data/{state.lower()}_blockgroups.gpkg",
        county_blockgroups = f"data/spatial_data/{config['county'].lower()}_blockgroups.gpkg"
    script: "scripts/retrieve_census_data.py"

rule retrieve_project_sunroof:
    input: 
        blockgroups = f"data/spatial_data/{state.lower()}_blockgroups.gpkg"
    output: 
        project_sunroof = "data/spatial_data/project-sunroof-census_tract.csv",
        local_potential = f"data/spatial_data/{state.lower()}_rooftop_potential.gpkg"
    script: "scripts/retrieve_project_sunroof.py"

# a bespoke step to make this analysis specific to armourdale
rule retrieve_armourdale_shape:
    output: 
        armourdale = "data/spatial_data/armourdale_shape.gpkg"
    script: "scripts/retrieve_armourdale.py"

rule retrieve_electric_utility:
    input: 
      cutout="data/spatial_data/armourdale_shape.gpkg"
    output:
      utility="data/spatial_data/electric_utility.gpkg"
    script: "scripts/retrieve_electric_utility.py"

rule retrieve_usrdb:
    input: 
        utility="data/spatial_data/electric_utility.gpkg"
    output: 
        rates="data/usrdb_rates.csv"
    script: "scripts/retrieve_usrdb.py"

rule calculate_res_structures:
    input: 
        census_data = "data/spatial_data/county_census_data.gpkg",
        armourdale = "data/spatial_data/armourdale_shape.gpkg"
    output: 
        res_structures = "data/residential_buildings.csv"
    script: "scripts/calculate_res_structures.py"

rule retrieve_res_load:
    input:
        spatial_lut = "data/spatial_data/spatial_lut.csv"
    output: 
        elec_load = "data/timeseries/residential_elec_load.csv",
        heat_load = "data/timeseries/residential_heat_load.csv",
        weather = "data/timeseries/weather_year.csv"
    script: "scripts/retrieve_res_load.py"

rule retrieve_lead_data:
    input: 
        community = "data/spatial_data/armourdale_shape.gpkg",
        county_blockgroups = f"data/spatial_data/{config['county'].lower()}_blockgroups.gpkg"
    output: 
        lead_data = f"data/spatial_data/{state_abbr}-2018-LEAD-data/{state_abbr} AMI Census Tracts 2018.csv",
        lead_community = "data/spatial_data/armourdale_lead.csv"
    script: "scripts/retrieve_lead_data.py"
    
rule build_dag:
    input: "Snakefile"
    output:
        "dag.png"
    shell:
        "snakemake --dag | dot -Tpng > {output}"