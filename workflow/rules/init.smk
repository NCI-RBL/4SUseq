from snakemake.utils import validate
import pandas as pd
import yaml
import pprint
import os
from os.path import join

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

# configfile: "config.yaml"
# validate(config, schema="../schemas/config.schema.yaml")

# set memory
MEMORY="200"

WORKDIR = config["workdir"]
RESOURCESDIR = config["resourcesdir"]
SCRIPTSDIR = config["scriptsdir"]
# define subfolders
TRIMDIR = join(WORKDIR,"trim")
FASTUNIQDIR = join(WORKDIR,"fastuniq")
RAWFASTQDIR = join(WORKDIR,"raw_fastq")
if not os.path.exists(RAWFASTQDIR):
	os.mkdir(RAWFASTQDIR)

SAMPLESDF = pd.read_csv(config["samples"],sep="\t",header=0,index_col="sampleName")

SAMPLES = list(SAMPLESDF.index)
#now path to r1 fastq file for sampleA will be sampledf["path_to_R1_fastq"]["sampleA"]
# validate(SAMPLESDF, schema="../schemas/samples.schema.yaml")
for sample in SAMPLES:
	f1=os.path.basename(SAMPLESDF["path_to_R1_fastq"][sample])
	f2=os.path.basename(SAMPLESDF["path_to_R2_fastq"][sample])
	if not os.path.exists(join(RAWFASTQDIR,f1)):
		os.symlink(SAMPLESDF["path_to_R1_fastq"][sample],join(RAWFASTQDIR,f1))
	if not os.path.exists(join(RAWFASTQDIR,f2)):
		os.symlink(SAMPLESDF["path_to_R2_fastq"][sample],join(RAWFASTQDIR,f2))
# print(SAMPLESDF)

## Load tools from YAML file
with open(config["tools"]) as f:
	TOOLS = yaml.safe_load(f)
# pprint.pprint(tools)


if not os.path.exists(TRIMDIR):
	os.mkdir(TRIMDIR)
if not os.path.exists(FASTUNIQDIR):
	os.mkdir(FASTUNIQDIR)

GENOME=config["genome"]
GTF=config["gtf"][GENOME]
