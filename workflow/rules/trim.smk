def get_fastqs(wildcards):
	d=dict()
	d["R1"]=SAMPLESDF["path_to_R1_fastq"][wildcards.sample]
	d["R2"]=SAMPLESDF["path_to_R2_fastq"][wildcards.sample]
#	print("in get_fastqs")
#	print(d)
	return d


rule cutadapt:
	input:
		unpack(get_fastqs)
	output:
		of1=join(TRIMDIR,"{sample}.R1.trim.fastq.gz"),
		of2=join(TRIMDIR,"{sample}.R2.trim.fastq.gz")
	params:
		sample="{sample}",
		adapters=join(RESOURCESDIR,"TruSeq_and_nextera_adapters.consolidated.fa")
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