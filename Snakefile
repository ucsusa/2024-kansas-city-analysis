configfile: "config.yml"

from us import states
from pathlib import Path
from dotenv import load_dotenv

state = config['state']
state_abbr = states.lookup(state).abbr

county_name = config['county']
community_name = config['community_name']

env_file = Path("./.env").resolve()
load_dotenv(str(env_file))

rule targets:
    input: 
        community = f"data/spatial_data/{community_name.lower()}_shape.gpkg",
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
        res_energy_expenses = f"data/{community_name.lower()}_energy_expenses.csv",
        zoning_data = f"data/spatial_data/{community_name.lower()}/zoning.gpkg",
        rescaled_elec_load = "data/timeseries/residential_elec_load_rescaled.csv",
        costs = "data/technology_costs.csv",
        residential_model = f"data/networks/{community_name}_residential.nc",
        dag = "dag.png"

rule build_residential_model:
    input:
        costs = "data/technology_costs.csv",
        utility_costs = "data/utility_technology_costs.csv",
        rescaled_elec_load = "data/timeseries/residential_elec_load_rescaled.csv",
        weather = "data/timeseries/weather_year.csv",
    output:
        residential_model = f"data/networks/{community_name}_residential.nc"
    script:
        "scripts/build_residential_model.py"

rule solve_residential_model:
    input:
        residential_model = f"data/networks/{community_name}_residential.nc"
    output:
        residential_model_solved = f"data/networks/{community_name}_residential_solved.nc"
    script:
        "scripts/solve_residential_model.py"

rule run_sensitivity:
    input:
        residential_model = f"data/networks/{community_name}_residential.nc"
    output:
        sensitivity_results = f"results/full_sensitivity_data.csv"
    script:
        "scripts/run_sensitivity.py"
        

rule retrieve_supply_regions:
    input:
        script = "scripts/retrieve_supply_regions.py",
        community = f"data/spatial_data/{community_name.lower()}_shape.gpkg"
    output: 
        supply_regions = "data/spatial_data/supply_regions.shp"
    script: "scripts/retrieve_supply_regions.py"

rule retrieve_spatial_lut:
    output: 
        spatial_lut = "data/spatial_data/spatial_lut.csv"
    script: "scripts/retrieve_lut.py"

rule retrieve_outage_data:
    input:
        "scripts/retrieve_outage_data.py"
    output: 
        outages = "data/timeseries/outages.csv",
        county_outages = f"data/timeseries/{county_name.lower()}_outages.csv"
    script: f"{input}"

rule retrieve_census_data:
    output:
        census_data = "data/spatial_data/county_census_data.gpkg",
        state_blockgroups = f"data/spatial_data/{state.lower()}_blockgroups.gpkg",
        county_blockgroups = f"data/spatial_data/{county_name.lower()}_blockgroups.gpkg"
    script: "scripts/retrieve_census_data.py"

rule retrieve_project_sunroof:
    input: 
        blockgroups = f"data/spatial_data/{state.lower()}_blockgroups.gpkg",
        community = f"data/spatial_data/{community_name.lower()}_shape.gpkg"
    output: 
        project_sunroof = "data/spatial_data/project-sunroof-census_tract.csv",
        local_potential = f"data/spatial_data/{community_name.lower()}_rooftop_potential.gpkg"
    script: "scripts/retrieve_project_sunroof.py"

# a bespoke step to make this analysis specific to community
rule retrieve_community_shape:
    output: 
        community = f"data/spatial_data/{community_name.lower()}_shape.gpkg"
    script: "scripts/retrieve_community_cutout.py"

rule retrieve_electric_utility:
    input: 
      cutout=f"data/spatial_data/{community_name.lower()}_shape.gpkg"
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
        community = f"data/spatial_data/{community_name.lower()}_shape.gpkg"
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
        community = f"data/spatial_data/{community_name.lower()}_shape.gpkg",
        county_blockgroups = f"data/spatial_data/{config['county'].lower()}_blockgroups.gpkg"
    output: 
        lead_data = f"data/spatial_data/{state_abbr}-2018-LEAD-data/{state_abbr} AMI Census Tracts 2018.csv",
        lead_community = f"data/spatial_data/{community_name.lower()}_lead.csv"
    script: "scripts/retrieve_lead_data.py"

rule retrieve_nrel_costs:
    output: 
        costs = "data/technology_costs.csv"
    script: "scripts/retrieve_nrel_costs.py"

rule retrieve_renewable_profiles:
    input:
        supply_regions = "data/spatial_data/supply_regions.shp"
    output:
        solar = "data/timeseries/solar.csv"
    script: "scripts/retrieve_renewables.py"

rule calculate_historical_expenses:
    input:
        lead_community = f"data/spatial_data/{community_name.lower()}_lead.csv" 
    output: 
        res_energy_expenses = f"data/{community_name.lower()}_energy_expenses.csv"
    script: "scripts/calculate_historical_expenses.py"

rule retrieve_community_spatial_data:
    input: 
        community = f"data/spatial_data/{community_name.lower()}_shape.gpkg"
    output: 
        zoning_data = f"data/spatial_data/{community_name.lower()}/zoning.gpkg"
    script: "scripts/retrieve_shapefiles.py"

rule calculate_rescaled_load:
    input: 
        res_energy_expenses = f"data/{community_name.lower()}_energy_expenses.csv",
        elec_load = "data/timeseries/residential_elec_load.csv",
        heat_load = "data/timeseries/residential_heat_load.csv",
        res_structures = "data/residential_buildings.csv"
    output: 
        rescaled_elec_load = "data/timeseries/residential_elec_load_rescaled.csv",
        rescaled_heat_load = "data/timeseries/residential_heat_load_rescaled.csv",
    script: "scripts/calculate_residential_load.py"
    
rule build_dag:
    input: "Snakefile"
    output:
        "dag.png"
    shell:
        "snakemake --dag | dot -Tpng > {output}"