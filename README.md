# 4SUseq
Analysis of 4SU(4-thiouridine) labeled (and unlabeled) total RNA samples. This snakemake workflow is designed specifically for analysis of samples generated under project **RBL-4**

## RBL-4
### Description: 

To investigate newly synthesized transcripts, RNAs in mouse stem cells (WT and Ro-/-) are labeled with 4SU, extracted and treated to generate T to C mutations in the final library following slam-seq protocol. We would like to compare the labeled RNA levels in WT and Ro-/-.

### Pipeline:

Snakemake is used to orchestrate a workflow illustrated by the following schematic:
 ![img](https://github.com/RBL-NCI/4SUseq/blob/main/4sUSeq_pipeline.png)

### Resources:
Genome supported:

* hg38
* mm10

Required resources(hisat2 indices etc.) are created using the following scripts and stored on biowulf at `/data/RBL_NCI/Wolin/mESC_slam_analysis/resource`

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
```
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
