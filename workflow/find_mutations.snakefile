# import stuff
from snakemake.utils import R
from os.path import join
from os import listdir
import os, sys

# init.smk already defines
# SAMPLES
# SAMPLESDF
# TOOLS

include: "rules/init.smk"
include: "rules/trim.smk"
include: "rules/find_mutations.smk"

report: "report/workflow.rst"

rule all:
	input:
		# trimmed fastqs
		expand(join(TRIMDIR,"{sample}.R1.trim.fastq.gz"), sample=SAMPLES),
		expand(join(TRIMDIR,"{sample}.R2.trim.fastq.gz"), sample=SAMPLES),
		# fastuniq fastqs
		expand(join(FASTUNIQDIR,"{sample}.R1.trim.fastuniq.fastq.gz"), sample=SAMPLES),
		expand(join(FASTUNIQDIR,"{sample}.R2.trim.fastuniq.fastq.gz"), sample=SAMPLES),
		# # lanes.txt
		join(RAWFASTQDIR,"lanes.txt"),
		join(TRIMDIR,"lanes.txt"),
		join(FASTUNIQDIR,"lanes.txt"),
		# # qc-->fastqscreen
		expand(join(WORKDIR,"qc","fastqscreen","{sample}.R1.trim_screen.txt"), sample=SAMPLES),
		expand(join(WORKDIR,"qc","fastqscreen","{sample}.R2.trim_screen.txt"), sample=SAMPLES),
		# # hisat on fastuniq reads
		expand(join(WORKDIR,"hisat2","{sample}.bam"), sample=SAMPLES),
		expand(join(WORKDIR,"hisat2","{sample}.plus.bam"), sample=SAMPLES),
		expand(join(WORKDIR,"hisat2","{sample}.minus.bam"), sample=SAMPLES),
		# # call snps with bcftools
		expand(join(WORKDIR,"vcf","{sample}.plus.vcf.gz"), sample=SAMPLES),
		expand(join(WORKDIR,"vcf","{sample}.minus.vcf.gz"), sample=SAMPLES),
		# # split bam to mutated and unmutated BAMs
		expand(join(WORKDIR,"bams","{sample}.mutated.bam"), sample=SAMPLES),
		expand(join(WORKDIR,"bams","{sample}.unmutated.bam"), sample=SAMPLES),	
		# nfragments qc table
		# join(WORKDIR,"qc","nfragments","nfragments.tsv")