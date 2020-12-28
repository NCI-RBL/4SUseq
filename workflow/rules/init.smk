from snakemake.utils import validate
import pandas as pd
import yaml
import pprint

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
# validate(config, schema="../schemas/config.schema.yaml")

SAMPLESDF = pd.read_csv(config["samples"],sep="\t",header=0,index_col="sampleName")
SAMPLES = list(SAMPLESDF.index)
#now path to r1 fastq file for sampleA will be sampledf["path_to_R1_fastq"]["sampleA"]
# validate(SAMPLESDF, schema="../schemas/samples.schema.yaml")

## Load tools from YAML file
with open(config["tools"]) as f:
	TOOLS = yaml.safe_load(f)
# pprint.pprint(tools)

