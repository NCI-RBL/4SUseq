from snakemake.utils import validate
import pandas as pd
import yaml
import pprint
import os
from os.path import join
import json

def check_existence(filename):
    """Checks if file exists on filesystem
    :param filename <str>: Name of file to check
    """
    filename=filename.strip()
    if not os.path.exists(filename):
        sys.exit("File: {} does not exists!".format(filename))
    return True


def check_readaccess(filename):
    """Checks permissions to see if user can read a file
    :param filename <str>: Name of file to check
    """
    filename=filename.strip()
    check_existence(filename)
    if not os.access(filename,os.R_OK):
        sys.exit("File: {} exists, but user cannot read from file due to permissions!".format(filename))
    return True


def check_writeaccess(filename):
    """Checks permissions to see if user can write to a file
    :param filename <str>: Name of file to check
    """
    filename=filename.strip()
    check_existence(filename)
    if not os.access(filename,os.W_OK):
        sys.exit("File: {} exists, but user cannot write to file due to permissions!".format(filename))
    return True

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

# configfile: "config/config.yaml"
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
check_readaccess(GTF)
print("# ANNOTATION FILE :",GTF)
REFFA=config["reffa"][GENOME]
check_readaccess(REFFA)
print("# REFERENCE FASTA :",REFFA)
RSEMINDEXDIR=join(WORKDIR,"rsem_index")

#########################################################
# READ CLUSTER PER-RULE REQUIREMENTS
#########################################################

## Load cluster.json
try:
    CLUSTERJSON = config["clusterjson"]
except KeyError:
    CLUSTERJSON = join(WORKDIR,"cluster.json")
check_readaccess(CLUSTERJSON)
with open(CLUSTERJSON) as json_file:
    CLUSTER = json.load(json_file)

## Create lambda functions to allow a way to insert read-in values
## as rule directives
getthreads=lambda rname:int(CLUSTER[rname]["threads"]) if rname in CLUSTER and "threads" in CLUSTER[rname] else int(CLUSTER["__default__"]["threads"])
getmemg=lambda rname:CLUSTER[rname]["mem"] if rname in CLUSTER else CLUSTER["__default__"]["mem"]
getmemG=lambda rname:getmemg(rname).replace("g","G")
#########################################################