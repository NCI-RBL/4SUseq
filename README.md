# 4SUseq
Analysis of 4SU(4-thiouridine) labeled (and unlabeled) total RNA samples. This snakemake workflow is designed specifically for analysis of samples generated under project **RBL-4**

## RBL-4
### Description: 

To investigate newly synthesized transcripts, RNAs in mouse stem cells (WT and Ro-/-) are labeled with 4SU, extracted and treated to generate T to C mutations in the final library following slam-seq protocol. We would like to compare the labeled RNA levels in WT and Ro-/-.

### Pipeline:

Snakemake is used to orchestrate a workflow illustrated by the following schematic:
 ![img](https://github.com/RBL-NCI/4SUseq/blob/main/4sUSeq_pipeline.png)

> **NOTE**: The pink portion of the above figure is not implemented, i.e, there is no distinction between labeled and unlabeled samples. T-to-C mutations are called for all samples and there is no correction of labeled samples by subtracting away mutations seen in unlabeled samples.

### Location:

Latest version of the pipeline has been checked out at `/data/RBL_NCI/Pipelines/4sUSeq`

### How to run?

To run an already checked out version of the pipeline, one may run the following command (replace `v0.3.0` with the version of interest):

```bash
%  bash /data/RBL_NCI/Pipelines/4sUSeq/v0.3.0/find_mutations.sh
Pipeline Dir: /gpfs/gsfs11/users/RBL_NCI/Pipelines/4sUSeq/v0.3.0
Git Commit/Tag: 6fd5086f25bc66c406f574014754ecac4c6974af	v0.3.0
find_mutations.sh: find T-to-C mutations
USAGE:
  bash find_mutations.sh <MODE> <path_to_workdir>
Required Positional Argument:
  MODE: [Type: Str] Valid options:
    a) init <path_to_workdir> : initialize workdir
    b) run <path_to_workdir>: run with slurm
    c) reset <path_to_workdir> : DELETE workdir dir and re-init it
    e) dryrun <path_to_workdir> : dry run snakemake to generate DAG
    f) unlock <path_to_workdir> : unlock workdir if locked by snakemake
    g) runlocal <path_to_workdir>: run without submitting to sbatch
```

The command has 2 arguments:

* **mode**: This defines if you want to initialize, reset, dryrun, unlock, runlocal(debug) or run the pipeline
* **workdir**: This defines the absolute path to the output folder

The steps involved in firing up the workflow will be:

1. **init**: initialize the workdir. This sets up the output directory and gets it ready to run.

```bash
%  bash /data/RBL_NCI/Pipelines/4sUSeq/v0.3.0/find_mutations.sh init /path/to/workdir
Pipeline Dir: /gpfs/gsfs11/users/RBL_NCI/Pipelines/4sUSeq/v0.3.0
Git Commit/Tag: 6fd5086f25bc66c406f574014754ecac4c6974af	v0.3.0
Working Dir: /path/to/workdir
/gpfs/gsfs11/users/RBL_NCI/Pipelines/4sUSeq/v0.3.0
/path/to/workdir
Logs Dir: /path/to/workdir/logs
Stats Dir: /path/to/workdir/stats
Done Initializing /path/to/workdir. You can now edit /path/to/workdir/config.yaml and /path/to/workdir/samples.tsv
```

2. **samples.tsv**. Editing this file in the output directory is required to define the sample manifest.  This file requires a header and 3 tab-delimited columns like:

```bash
sampleName	path_to_R1_fastq	path_to_R2_fastq
KO1_4SU_labeled	/pathto/KO1_4SU_labeled.R1.fastq.gz	/pathto/KO1_4SU_labeled.R2.fastq.gz
```

3. **dryrun**: Once the sample manifest is set, we can now dryrun the pipeline. This ensures that
   * samples.tsv defined sample datasets are "readable"
   * outputfolder is "writable"
   * create a DAG to illustrate flow-of-jobs with

```bash
%  bash /data/RBL_NCI/Pipelines/4sUSeq/v0.3.0/find_mutations.sh dryrun /path/to/workdir
```

4. **run**: Once we have confirmed that the DAG looks good, we can now submit the pipeline to the cluster:

```bash
%  bash /data/RBL_NCI/Pipelines/4sUSeq/v0.3.0/find_mutations.sh run /path/to/workdir
```

### Resources:
Genome and annotation versions supported:

| Genome | Gencode Version |
| ------ | --------------- |
| hg38   | Release 36      |
| mm10   | Release 25      |

Required resources(hisat2 indices etc.) are created using the following scripts and stored on biowulf at `/data/RBL_NCI/Wolin/mESC_slam_analysis/resources`

__hg38__
```
# download files
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.primary_assembly.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz
gzip -d gencode.v36.primary_assembly.annotation.gtf.gz
gzip -d GRCh38.primary_assembly.genome.fa.gz
# rename files
ln GRCh38.primary_assembly.genome.fa hg38.fa
ln gencode.v36.primary_assembly.annotation.gtf hg38.v36.gtf
# build index
swarm -f build_index.sh -t 2 -g 100 --partition=ccr,norm --time=24:00:00
## build_index.sh:
## module load hisat/2.2.1.0 && hisat2-build hg38.fa hg38 > hisat2-build.log 2>&
# convert gtf to splicesites
python /usr/local/apps/hisat/2.2.1.0/hisat2_extract_splice_sites.py hg38.v36.gtf > hg38.v36.splicesites.txt
```

__mm10__
```
# download files
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
gzip -d gencode.vM25.primary_assembly.annotation.gtf.gz
gzip -d GRCm38.primary_assembly.genome.fa.gz
# rename files
ln GRCm38.primary_assembly.genome.fa mm10.fa
ln gencode.vM25.primary_assembly.annotation.gtf mm10.v25.gtf
# build index
swarm -f build_index.sh -t 2 -g 100 --partition=ccr,norm --time=24:00:00
## build_index.sh:
## module load hisat/2.2.1.0 && hisat2-build mm10.fa mm10 > hisat2-build.log 2>&1
# convert gtf to splicesites
python /usr/local/apps/hisat/2.2.1.0/hisat2_extract_splice_sites.py mm10.v25.gtf > mm10.v25.splicesites.txt
```

The `splicesites.txt` files are moved to the `resources` folder of this pipeline.

For adding annotations to the vcf files a `tab` file with gene-level annotations is generated using `create_tab.sh`:
```
python gtf2tab.py mm10.v25.gtf > mm10.genes.v25.tab
bgzip mm10.genes.v25.tab
tabix -s1 -b2 -e3 mm10.genes.v25.tab.gz
```
and the `gtf2tab.py` script:
```bash
import sys
def get_feature(feature,col9):
	col9=col9.split()
	for i,x in enumerate(col9):
		if x == feature:
			j=i+1
			break
	value=col9[j]
	value=value.replace("\"","")
	value=value.replace(";","")
	value=value.strip()
	return value

hdr="#CHROM FROM TO strand gene_id gene_type gene_name"
hdr=hdr.replace(" ","\t")
print(hdr)
fh=open(sys.argv[1],'r')
lines=fh.readlines()
fh.close()
for l in lines:
	if l.startswith("#"):
		continue
	l=l.strip().split("\t")
	if l[2] != "gene":
		continue
	chrom=l[0]
	start=l[3]
	end=l[4]
	strand=l[6]
	gene_id=get_feature("gene_id",l[8])
	gene_type=get_feature("gene_type",l[8])
	gene_name=get_feature("gene_name",l[8])
	x=[chrom,start,end,strand,gene_id,gene_type,gene_name]
	print("\t".join(x))
```

> **NOTE**: 
>
> "Exon" and "Gene" level counts are reported for each sample. What does this mean?
>
> * "Exon" level counts actually means counts at the gene level but only including reads aligning within annotated exons, i.e., excluding reads aligned in the intronic regions. These can be thought of as counts to spliced mRNA.
> * "Gene" level counts also depict gene level counts, but including reads aligned in the intronic region along with reads aligned in the exons. These can be thought of as counts to pre-spliced mRNA.
