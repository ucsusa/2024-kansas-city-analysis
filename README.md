# 2024 Kansas City Analysis
This repository holds analysis for the energy system in Kansas City, Kansas. Located in Wyandotte County, Kansas.


# Installation

## Requirements

* `git` - version control software
    * [Windows Installation instructions](https://git-scm.com/download/win)
    * [MacOS Installation instructions](https://git-scm.com/download/mac)
    * [Linux Installation instructions](https://git-scm.com/download/linux)
* Python installed with either `conda` or `mamba`(recommended)
    * Download `mamba` installer [here](https://github.com/conda-forge/miniforge).
    * 'anaconda' ('conda') installation instructions [here](https://docs.anaconda.com/anaconda/install/windows/).

> [!NOTE]
> Make sure you add Python to PATH during installation.

## Installation Steps
0. Open command prompt or terminal window. Copy and paste the following commands.

1. Clone the repository

```bash
git clone https://github.com/ucsusa/2024-kansas-city-analysis.git
```

2. Set up the environment

```bash
cd 2024-kansas-city-analysis
mamba env create  # mamba and conda may be used interchangeably, here
mamba activate kansas-city
```

3. Creating the `.env` file

Users should copy the `.env.template` file into a new file simply called `.env`.
This file contains "secret" information, such as API keys, emails, and other data
that should remain local. In order to run the current model, users must have API keys
from the following organizations:

* [U.S. Census API](https://api.census.gov/data/key_signup.html)

These keys may be added directly to the `.env` file.    

## Running the model 

This project uses the workflow management tool, `snakemake`, to create a reproducible data pipeline.
Running the command

```bash
snakemake --cores=1
```

will run the workflow illustrated in the directed acyclic graph (DAG) shown below.