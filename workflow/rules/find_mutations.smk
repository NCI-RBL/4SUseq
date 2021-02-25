rule get_fastuniq_readids:
	input:
		if1=rules.cutadapt.output.of1,
		if2=rules.cutadapt.output.of2
	output:
		of1=join(FASTUNIQDIR,"{sample}.R1.trim.fastuniq.fastq.gz"),
		of2=join(FASTUNIQDIR,"{sample}.R2.trim.fastuniq.fastq.gz"),
		readids=join(FASTUNIQDIR,"{sample}.fastuniq.readids"),
	params:
		sample="{sample}",
		workdir=WORKDIR,
		outdir=FASTUNIQDIR,
		fastuniq=join(RESOURCESDIR,"fastuniq")
	envmodules: TOOLS["pigz"]["version"]
	threads: 56
	shell:"""
zcat {input.if1} > /lscratch/${{SLURM_JOBID}}/{params.sample}.R1.fastq
zcat {input.if2} > /lscratch/${{SLURM_JOBID}}/{params.sample}.R2.fastq
echo "/lscratch/${{SLURM_JOBID}}/{params.sample}.R1.fastq" > /lscratch/${{SLURM_JOBID}}/{params.sample}.list
echo "/lscratch/${{SLURM_JOBID}}/{params.sample}.R2.fastq" >> /lscratch/${{SLURM_JOBID}}/{params.sample}.list
{params.fastuniq} -i /lscratch/${{SLURM_JOBID}}/{params.sample}.list -t q -o /lscratch/${{SLURM_JOBID}}/{params.sample}.R1.trim.fastuniq.fastq -p /lscratch/${{SLURM_JOBID}}/{params.sample}.R2.trim.fastuniq.fastq -c 0
awk '{{if (! ((FNR + 3) % 4)) {{ print(substr($1,2))}}}}' /lscratch/${{SLURM_JOBID}}/{params.sample}.R1.trim.fastuniq.fastq > {output.readids}
pigz -p4 /lscratch/${{SLURM_JOBID}}/{params.sample}.R1.trim.fastuniq.fastq && mv /lscratch/${{SLURM_JOBID}}/{params.sample}.R1.trim.fastuniq.fastq.gz {output.of1}
pigz -p4 /lscratch/${{SLURM_JOBID}}/{params.sample}.R2.trim.fastuniq.fastq && mv /lscratch/${{SLURM_JOBID}}/{params.sample}.R2.trim.fastuniq.fastq.gz {output.of2}

echo "DONE!"
"""

rule get_fastq_nreads:
	input:
		f1=RAWFASTQDIR,
		f2=TRIMDIR,
		f3=FASTUNIQDIR,
		x=expand(join(FASTUNIQDIR,"{sample}.R1.trim.fastuniq.fastq.gz"),sample=SAMPLES)
	output:
		o1=join(RAWFASTQDIR,"lanes.txt"),
		o2=join(TRIMDIR,"lanes.txt"),
		o3=join(FASTUNIQDIR,"lanes.txt")
	params:
		workdir=WORKDIR,
		outdir=FASTUNIQDIR,
		bashscript=join(SCRIPTSDIR,"get_lanes.sh"),
		scriptsdir=SCRIPTSDIR
	threads: 56
	shell:"""
for f in {input.f1} {input.f2} {input.f3};do
cd $f
bash {params.bashscript}
done
"""

rule fastq_screen:
	input:
		if1=rules.cutadapt.output.of1,
		if2=rules.cutadapt.output.of2
	output:
		out1=join(WORKDIR,"qc","fastqscreen","{sample}.R1.trim_screen.txt"),
		out2=join(WORKDIR,"qc","fastqscreen","{sample}.R2.trim_screen.txt")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		outdir=join(WORKDIR,"qc","fastqscreen"),
		conf=join(RESOURCESDIR,"fastq_screen.conf")
	threads: 56
	envmodules: TOOLS["fastq_screen"]["version"], TOOLS["bowtie"]["version"]
	shell:"""
fastq_screen --conf {params.conf} \
 --outdir "{params.outdir}" \
 --threads {threads} --subset 1000000 \
 --aligner bowtie2 --force \
 {input.if1} \
 {input.if2}
"""	

rule hisat:
	input:
		if1=rules.cutadapt.output.of1,
		if2=rules.cutadapt.output.of2
	output:
		bam=join(WORKDIR,"hisat2","{sample}.bam"),
		plusbam=join(WORKDIR,"hisat2","{sample}.plus.bam"),
		minusbam=join(WORKDIR,"hisat2","{sample}.minus.bam"),
		flagstat1=join(WORKDIR,"hisat2","{sample}.post_secondary_supplementary_filter.bam.flagstat"),
		flagstat2=join(WORKDIR,"hisat2","{sample}.post_insertion_filter.bam.flagstat"),
		flagstat3=join(WORKDIR,"hisat2","{sample}.bam.flagstat"),
		flagstat4=join(WORKDIR,"hisat2","{sample}.plus.bam.flagstat"),
		flagstat5=join(WORKDIR,"hisat2","{sample}.minus.bam.flagstat"),
	params:
		sample="{sample}",
		mem=MEMORY,
		workdir=WORKDIR,
		outdir=join(WORKDIR,"hisat2"),
		genome=config["genome"],
		hisatindex=config["hisatindexdir"],
		hisat_rna_strandness=config["hisat_rna_strandness"],
		splicesites=join(RESOURCESDIR,config["genome"]+".splicesites.txt"),
		filter_script=join(SCRIPTSDIR,"filter_bam.py"),
		mapqfilter=config['mapqfilter'],
		ninsertionfilter=config['ninsertionfilter']
	threads: 56
	envmodules: TOOLS["hisat"]["version"], TOOLS["samtools"]["version"], TOOLS["sambamba"]["version"],  TOOLS["picard"]["version"],  TOOLS["bbtools"]["version"]
	shell:"""
# Align with HiSAT2 and remove secondary/supplementary alignments
hisat2 \
 -x {params.hisatindex}/{params.genome}/{params.genome} \
 -1 {input.if1} \
 -2 {input.if2} \
 --rna-strandness {params.hisat_rna_strandness} \
 --known-splicesite-infile {params.splicesites} \
 --summary-file {params.outdir}/{params.sample}.hisat2.summary.log \
 --threads {threads} \
 --rg-id {params.sample} --rg SM:{params.sample} \
 --un-gz {params.outdir}/{params.sample}.unmapped.fastq.gz \
 --mp 4,2 \
 | samtools view -@{threads} -bS -F 4 -F 256 - > /dev/shm/{params.sample}.post_secondary_supplementary_filter.tmp.bam

# Sort the "post_secondary_supplementary_filter" BAM and collect some stats
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out=/dev/shm/{params.sample}.post_secondary_supplementary_filter.bam /dev/shm/{params.sample}.post_secondary_supplementary_filter.tmp.bam && rm -f /dev/shm/{params.sample}.post_secondary_supplementary_filter.tmp.bam
sambamba flagstat --nthreads={threads} /dev/shm/{params.sample}.post_secondary_supplementary_filter.bam > {output.flagstat1}

# Apply the "number of insertions in read" filter
python {params.filter_script} -i /dev/shm/{params.sample}.post_secondary_supplementary_filter.bam -o /dev/shm/{params.sample}.post_insertion_filter.tmp.bam -q 0 -n {params.ninsertionfilter}

# Sort, FixMateInfo, collect stats
sambamba sort -n --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out=/dev/shm/{params.sample}.post_insertion_filter.bam /dev/shm/{params.sample}.post_insertion_filter.tmp.bam && rm -f /dev/shm/{params.sample}.post_insertion_filter.tmp.bam
java -Xmx{params.mem}g -jar $PICARDJARPATH/picard.jar FixMateInformation I=/dev/shm/{params.sample}.post_insertion_filter.bam O=/dev/stdout ASSUME_SORTED=true QUIET=true \
 | sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out=/dev/shm/{params.sample}.postsambambasort.bam /dev/stdin
sambamba flagstat --nthreads={threads} /dev/shm/{params.sample}.postsambambasort.bam > {output.flagstat2}

# Apply the MAPQ filter and remove "widowed" reads
python {params.filter_script} -i /dev/shm/{params.sample}.postsambambasort.bam -o {output.bam} -q {params.mapqfilter} -n 1000

# Index and collect stats
sambamba index --nthreads={threads} {output.bam}
sambamba flagstat --nthreads={threads} {output.bam} > {output.flagstat3}

# Split reads into :
# 1. Reads originating from + strand fragments --> these will go for T-to-C mutation detection
#         R2
#       ------>
#     5'-------------------------------------------------3'
#     3'-------------------------------------------------5'
#                                                 <------
#                                                    R1
# 2. Reads originating from - strand fragments --> these will go for A-to-G mutation detection
#         R1
#       ------>
#     5'-------------------------------------------------3'
#     3'-------------------------------------------------5'
#                                                 <------
#                                                    R2

samtools view -@ {threads} -b -f 128 -F 16 {output.bam} > /dev/shm/{params.sample}.F1.bam
samtools view -@ {threads} -b -f 80 {output.bam} > /dev/shm/{params.sample}.F2.bam
samtools merge -c -f -p -@ {threads} /dev/shm/{params.sample}.F.bam /dev/shm/{params.sample}.F1.bam /dev/shm/{params.sample}.F2.bam && rm -f /dev/shm/{params.sample}.F1.bam /dev/shm/{params.sample}.F2.bam
samtools view -@ {threads} -b -f 144 {output.bam} > /dev/shm/{params.sample}.R1.bam
samtools view -@ {threads} -b -f 64 -F 16 {output.bam} > /dev/shm/{params.sample}.R2.bam
samtools merge -c -f -p -@ {threads} /dev/shm/{params.sample}.R.bam /dev/shm/{params.sample}.R1.bam /dev/shm/{params.sample}.R2.bam && rm -f /dev/shm/{params.sample}.R1.bam /dev/shm/{params.sample}.R2.bam

# Sort and collect stats for plus strand BAM
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.plusbam} /dev/shm/{params.sample}.F.bam && rm -f /dev/shm/{params.sample}.F.bam
sambamba flagstat --nthreads={threads} {output.plusbam} > {output.flagstat4}

# Sort and collect stats for minus strand BAM
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.minusbam} /dev/shm/{params.sample}.R.bam && rm -f /dev/shm/{params.sample}.R.bam
sambamba flagstat --nthreads={threads} {output.minusbam} > {output.flagstat5}

"""

rule create_toSNPcalling_BAM:
	input:
		plusbam=rules.hisat.output.plusbam,
		minusbam=rules.hisat.output.minusbam,
		readids=rules.get_fastuniq_readids.output.readids
	output:
		plusbam=join(WORKDIR,"hisat2","{sample}.plus.toSNPcalling.bam"),
		minusbam=join(WORKDIR,"hisat2","{sample}.minus.toSNPcalling.bam"),
	params:
		sample="{sample}",
		mem=MEMORY,
		workdir=WORKDIR,
		script=join(SCRIPTSDIR,"filter_bam_by_readids.py")
	envmodules: TOOLS["sambamba"]["version"]
	threads: 4
	shell:"""
python {params.script} -i {input.plusbam} -o /dev/shm/{params.sample}.plus.bam --readids {input.readids}
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.plusbam} /dev/shm/{params.sample}.plus.bam && rm -f /dev/shm/{params.sample}.plus.bam
python {params.script} -i {input.minusbam} -o /dev/shm/{params.sample}.minus.bam --readids {input.readids}
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.minusbam} /dev/shm/{params.sample}.minus.bam && rm -f /dev/shm/{params.sample}.minus.bam
"""



rule call_mutations:
	input:
		plusbam=rules.create_toSNPcalling_BAM.output.plusbam,
		minusbam=rules.create_toSNPcalling_BAM.output.minusbam
	output:
		plusvcf=join(WORKDIR,"vcf","{sample}.plus.vcf.gz"),
		minusvcf=join(WORKDIR,"vcf","{sample}.minus.vcf.gz"),
		vcf=join(WORKDIR,"vcf","{sample}.vcf.gz")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		outdir=join(WORKDIR,"vcf"),
		genome=config["genome"],
		hisatindex=config["hisatindexdir"],
		tab=join(RESOURCESDIR,config["genome"]+".genes.tab.gz"),
		hdr=join(RESOURCESDIR,"hdr.txt"),
		filter_script=join(SCRIPTSDIR,"filter_bam.py")
	threads: 56
	envmodules: TOOLS["bcftools"]["version"], TOOLS["samtools"]["version"]
	shell:"""
# call mutations with minimum 3 read-support
# plus strand
bcftools mpileup -f {params.hisatindex}/{params.genome}/{params.genome}.fa -a AD,ADF,ADR {input.plusbam} | \
bcftools call -mv -Ob --threads {threads} | \
bcftools filter -i '%QUAL>45' -g3 -Ob --threads {threads} - | \
bcftools view -i '(REF=="T" & ALT="C")' --threads {threads} - > /dev/shm/{params.sample}.TtoC.bcf
bcftools sort -T /dev/shm /dev/shm/{params.sample}.TtoC.bcf | bgzip > {output.plusvcf}
tabix -p vcf {output.plusvcf}
# minus strand
bcftools mpileup -f {params.hisatindex}/{params.genome}/{params.genome}.fa -a AD,ADF,ADR {input.minusbam} | \
bcftools call -mv -Ob --threads {threads} | \
bcftools filter -i '%QUAL>45' -g3 -Ob --threads {threads} - | \
bcftools view -i '(REF=="A" & ALT="G")' --threads {threads} - > /dev/shm/{params.sample}.AtoG.bcf
bcftools sort -T /dev/shm /dev/shm/{params.sample}.AtoG.bcf | bgzip > {output.minusvcf}
tabix -p vcf {output.minusvcf}
# merge mutations
bcftools merge --force-samples -O z -o {output.vcf} {output.plusvcf} {output.minusvcf}
tabix -p vcf {output.vcf}
"""


rule split_bam_by_mutation:
	input:
		plusbam=rules.hisat.output.plusbam,
		minusbam=rules.hisat.output.minusbam,
		plusvcf=rules.call_mutations.output.plusvcf,
		minusvcf=rules.call_mutations.output.minusvcf,
	output:
		mutatedbam=join(WORKDIR,"bams","{sample}.mutated.bam"),
		unmutatedbam=join(WORKDIR,"bams","{sample}.unmutated.bam"),
	params:
		sample="{sample}",
		mem=MEMORY,
		workdir=WORKDIR,
		outdir=join(WORKDIR,"filtered_reads"),
		genome=config["genome"],
		hisatindex=config["hisatindexdir"],
		tsv2readidspy=join(SCRIPTSDIR,"tsv2readids.py"),
		filterbyreadidspy=join(SCRIPTSDIR,"filter_bam_by_readids.py"),
		sam2tsvjar=join(RESOURCESDIR,"sam2tsv.jar")
	threads: 4
	envmodules: TOOLS["java"]["version"], TOOLS["sambamba"]["version"], TOOLS["samtools"]["version"]
	shell:"""
# using sam2tsv from jvarkit --> https://lindenb.github.io/jvarkit/Sam2Tsv.html

plusvcf=$(basename {input.plusvcf})
zcat {input.plusvcf} > /dev/shm/${{plusvcf%.*}}
java -Xmx{params.mem}g -jar {params.sam2tsvjar} \
 --reference {params.hisatindex}/{params.genome}/{params.genome}.fa \
 --skip-N \
 {input.plusbam} | \
python {params.tsv2readidspy} /dev/shm/${{plusvcf%.*}} "T" "C" | \
sort | uniq > /dev/shm/{sample}.plus.readids

minusvcf=$(basename {input.minusvcf})
zcat {input.minusvcf} > /dev/shm/${{minusvcf%.*}}
java -Xmx{params.mem}g -jar {params.sam2tsvjar} \
 --reference {params.hisatindex}/{params.genome}/{params.genome}.fa \
 --skip-N \
 {input.minusbam} | \
python {params.tsv2readidspy} /dev/shm/${{minusvcf%.*}} "A" "G" | \
sort | uniq > /dev/shm/{sample}.minus.readids

python {params.filterbyreadidspy} -i {input.plusbam} -o /dev/shm/{params.sample}.mutated.plus.bam --readids /dev/shm/{sample}.plus.readids -o2 /dev/shm/{params.sample}.unmutated.plus.bam
python {params.filterbyreadidspy} -i {input.minusbam} -o /dev/shm/{params.sample}.mutated.minus.bam --readids /dev/shm/{sample}.minus.readids -o2 /dev/shm/{params.sample}.unmutated.minus.bam

samtools merge -c -f -p -@ {threads} /dev/shm/{params.sample}.mutated.bam /dev/shm/{params.sample}.mutated.plus.bam /dev/shm/{params.sample}.mutated.minus.bam && rm -f /dev/shm/{params.sample}.mutated.plus.bam /dev/shm/{params.sample}.mutated.minus.bam
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.mutatedbam} /dev/shm/{params.sample}.mutated.bam && rm -f /dev/shm/{params.sample}.mutated.bam

samtools merge -c -f -p -@ {threads} /dev/shm/{params.sample}.unmutated.bam /dev/shm/{params.sample}.unmutated.plus.bam /dev/shm/{params.sample}.unmutated.minus.bam && rm -f /dev/shm/{params.sample}.unmutated.plus.bam /dev/shm/{params.sample}.unmutated.minus.bam
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.unmutatedbam} /dev/shm/{params.sample}.unmutated.bam && rm -f /dev/shm/{params.sample}.unmutated.bam

"""

# rule get_nfragments_json:
# 	input:
# 		lanes1=rules.get_fastq_nreads.output.o1,
# 		lanes2=rules.get_fastq_nreads.output.o2,
# 		lanes3=rules.get_fastq_nreads.output.o3,
# 		flagstat1=rules.hisat_on_fastuniq.output.flagstat1,
# 		flagstat2=rules.hisat_on_fastuniq.output.flagstat2,
# 		flagstat3=rules.hisat_on_fastuniq.output.flagstat3,
# 		vcf=rules.call_mutations.output.TtoCvcf
# 	output:
# 		json=join(WORKDIR,"qc","nfragments","{sample}.json")
# 	params:
# 		sample="{sample}",
# 		workdir=WORKDIR,
# 		outdir=join(WORKDIR,"qc","nfragments"),
# 		pyscript=join(SCRIPTSDIR,"get_per_sample_nfragments.py")
# 	shell:"""
# cd {params.workdir}
# python {params.pyscript} {params.sample} {output.json}
# """

# rule create_nfragments_table:
# 	input:
# 		expand(join(WORKDIR,"qc","nfragments","{sample}.json"),sample= SAMPLES)
# 	output:
# 		table=join(WORKDIR,"qc","nfragments","nfragments.tsv")
# 	params:
# 		workdir=WORKDIR,
# 		outdir=join(WORKDIR,"qc","nfragments"),
# 		pyscript=join(SCRIPTSDIR,"nfragments_json2table.py")
# 	shell:"""
# cd {params.outdir}
# python {params.pyscript} {output.table}
# """


