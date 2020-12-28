def get_fastqs(wildcards):
	d=dict()
	d["R1"]=SAMPLESDF["path_to_R1_fastq"][wildcards.sample]
	d["R2"]=SAMPLESDF["path_to_R2_fastq"][wildcards.sample]
	return d


rule cutadapt:
	input:
		unpack(get_fastqs)
	output:
		of1=join(TRIMDIR,"{sample}.R1.trim.fastq.gz"),
		of2=join(TRIMDIR,"{sample}.R2.trim.fastq.gz")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		adapters=join(RESOURCES_DIR,"TruSeq_and_nextera_adapters.consolidated.fa")
	envmodules: TOOLS["cutadapt"]["version"]
	threads: 56
	shell:"""
cutadapt --pair-filter=any \
--nextseq-trim=2 \
--trim-n \
-n 5 -O 5 \
-q 10,10 -m 35:35 \
-b file:{params.adapters} \
-B file:{params.adapters} \
-j {threads} \
-o {output.of1} -p {output.of2} \
{input.R1} {input.R2}
"""

# rule fastqc:
# 	input:
# 		expand(join(WORKDIR,"fastqs","{sample}.fastq.gz"),sample=SAMPLES)
# 	output:
# 		expand(join(WORKDIR,"qc","fastqc","{sample}_fastqc.zip"),sample=SAMPLES)
# 	params:
# 		outdir=join(WORKDIR,"qc","fastqc")
# 	threads: 16
# 	# envmodules: TOOLS["fastqc"]["version"]
# 	container: "nciccbr/ccbr_fastqc_v0.11.9:latest"
# 	log: join(WORKDIR,"logs","fastqc.log")
# 	shell:"""
# fastqc {input} -t {threads} -o {params.outdir}
# """

# rule pycoqc:
# 	input:
# 		# sequencing_summary=join(WORKDIR,"fastqs","{sample}.sequencing_summary.txt")
# 		sequencing_summary=rules.guppy.output.sequencing_summary
# 	output:
# 		pycoQChtml=join(WORKDIR,"qc","pycoQC","{sample}.pycoQC.html"),
# 		pycoQCjson=join(WORKDIR,"qc","pycoQC","{sample}.pycoQC.json"),
# 	params:
# 		outdir=join(WORKDIR,"qc","fastqc")
# 	# conda: "envs/pycoqc.yaml"
# 	container: "docker://nciccbr/ccbr_pycoqc_v2.5.0.23:latest"
# 	log: join(WORKDIR,"logs","{sample}.pycoQC.log")
# 	shell:"""
# pycoQC -f {input.sequencing_summary} -o {output.pycoQChtml} -j {output.pycoQCjson} 2>&1 |tee -a {log}
# """