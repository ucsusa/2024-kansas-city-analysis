# 2024 Kansas City Analysis
This repository holds analysis for the energy system in Kansas City, Kansas. Located in Wyandotte County, Kansas.

## Workflow

The flow of data through the modeling process is shown in the graph below.

![DAG](dag.png)

There are a few categories of steps:
* **Retrieve**: In a `retrieve` step, data are primarily downloaded and lightly processed (e.g., ensuring good formatting and data types).
* **Calculate**: In a `calculate` step, data are transformed through some calculation.
* *place holder for future additions*


## Steps

### `retrieve_census_data`
In this step, data from the U.S. Census Bureau are queried. The datasets gathered, here, are:
* Total population and
* the number and types of residential building units.

### `retrieve_armourdale_shape`
In this step, the "shape" of the community of interest is retrieved. This shape can be used as a cut-out
to subset other geospatial data later.

> [!NOTE]
> This data is specific to the particular community of Armourdale in Kansas City, Kansas. If you 
> wish to model a different community, should omit this step or replace it with a different shape.
> For example, by specifying a few census tracts.

### `retrieve_spatial_lut`
This step downloads the spatial lookup table (LUT) for NREL's ResStock datasets. The spatial LUT 
cross references census tracts, counties, and states with public use microdata areas (PUMAs). As
well as how the data are stored within NREL's models.

### `retreive_res_load`
Simulated building load data is collected from NREL's ResStock database in this step. Currently, 
the data collected are aggregated building data for the building types defined in the `config.yml` file.
Future versions may include an option to specify individual buildings.