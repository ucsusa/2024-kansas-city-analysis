configfile: "config.yml"

from pathlib import Path
env_file = Path("./.env").resolve()
from dotenv import load_dotenv
load_dotenv(str(env_file))

rule targets:
    input: 
        sfa = "data/timeseries/residential_load.csv",
        res_structures = "data/residential_buildings.csv",
        rates = "data/usrdb_rates.csv",
        project_sunroof = f"data/spatial_data/project-sunroof-census_tract.csv",
        dag = "dag.png"

rule retrieve_spatial_lut:
    output: 
        spatial_lut = "data/spatial_data/spatial_lut.csv"
    script: "scripts/retrieve_lut.py"

rule retrieve_census_data:
    output:
        census_data = "data/spatial_data/county_census_data.gpkg",
        state_blockgroups = f"data/spatial_data/{config['state'].lower()}_blockgroups.gpkg",
        county_blockgroups = f"data/spatial_data/{config['county'].lower()}_blockgroups.gpkg"
    script: "scripts/retrieve_census_data.py"

rule retrieve_project_sunroof:
    input: 
        blockgroups = f"data/spatial_data/{config['state'].lower()}_blockgroups.gpkg"
    output: 
        project_sunroof = "data/spatial_data/project-sunroof-census_tract.csv",
        local_potential = f"data/spatial_data/{config['state'].lower()}_rooftop_potential.gpkg"
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
        heat_load = "data/timeseries/residential_heat_load.csv"
    script: "scripts/retrieve_res_load.py"
    
rule build_dag:
    input: "Snakefile"
    output:
        "dag.png"
    shell:
        "snakemake --dag | dot -Tpng > {output}"