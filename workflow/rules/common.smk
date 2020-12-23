from snakemake.utils import validate
import pandas as pd
import yaml
import pprint

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

SAMPLESDF = pd.read_csv(config["samples"],sep="\t",header=0,index_col="sampleName")
SAMPLES = list(SAMPLESDF.index)
#now path to fast5 folder for sampleA will be sampledf["path_to_fast5_parent_folder"]["sampleA"]
# validate(SAMPLESDF, schema="../schemas/samples.schema.yaml")

## Load tools from TSV file
# tools = pd.read_csv(config["tools"], sep="\t",header=0,index_col="tool")
# validate(tools, schema="../schemas/tools.schema.yaml")
#now envmodule for loading guppy will be tools["version"]["guppy"]

## Load tools from YAML file
with open(config["tools"]) as f:
	TOOLS = yaml.safe_load(f)
# pprint.pprint(tools)


## Load project metadata from TSV file
# project = pd.read_csv(config["project"], sep="\t",header=0,index_col="key")
# workdir = project["value"]["workdir"]
## Project related metadata is directly moved into config.yaml file....no need for project.tsv
