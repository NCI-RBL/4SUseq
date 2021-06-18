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

rule star:
	input:
		if1=rules.cutadapt.output.of1,
		if2=rules.cutadapt.output.of2
	output:
		starbam=join(WORKDIR,"star","{sample}.Aligned.out.bam"),
		bam=join(WORKDIR,"star","{sample}.bam"),
		bamflagstat=join(WORKDIR,"star","{sample}.bam.flagstat"),
		tbam=join(WORKDIR,"star","{sample}.Aligned.toTranscriptome.out.bam"),
		plusbam=join(WORKDIR,"star","{sample}.plus.bam"),
		minusbam=join(WORKDIR,"star","{sample}.minus.bam"),
		postsecondarysupplementaryfilterbamflagstat=join(WORKDIR,"star","{sample}.post_secondary_supplementary_filter.bam.flagstat"),
		postinsertionfilterbamflagstat=join(WORKDIR,"star","{sample}.post_insertion_filter.bam.flagstat"),
		plusbamflagstat=join(WORKDIR,"star","{sample}.plus.bam.flagstat"),
		minusbamflagstat=join(WORKDIR,"star","{sample}.minus.bam.flagstat"),
	params:
		sample="{sample}",
		mem=MEMORY,
		workdir=WORKDIR,
		outdir=join(WORKDIR,"star"),
		genome=config["genome"],
		starindexdir=config["starindexdir"],
		gtf=GTF,
		filter_script=join(SCRIPTSDIR,"filter_bam.py"),
		mapqfilter=config['mapqfilter'],
		ninsertionfilter=config['ninsertionfilter'],
		outFilterIntronMotifs=config['starparams']['outFilterIntronMotifs'],
		outFilterType=config['starparams']['outFilterType'],
	threads: 56
	envmodules: TOOLS["star"]["version"], TOOLS["samtools"]["version"],  TOOLS["picard"]["version"], 
	shell:"""
tmpdir=/lscratch/$SLURM_JOB_ID
cd {params.outdir}
# Align with STAR and remove secondary/supplementary alignments
readlength=$(
	zcat {input.if1} {input.if2} | \
	awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; \
	END {{print maxlen-1}}'
)
STAR --genomeDir {params.starindexdir} \
	--outFilterIntronMotifs {params.outFilterIntronMotifs} \
	--outFilterType {params.outFilterType} \
	--readFilesIn {input.if1} {input.if2} \
	--readFilesCommand zcat \
	--runThreadN {threads} \
	--outFileNamePrefix {params.sample}. \
	--twopassMode Basic \
	--sjdbGTFfile {params.gtf} \
	--quantMode TranscriptomeSAM GeneCounts \
	--outSAMtype BAM Unsorted \
	--alignEndsProtrude 10 ConcordantPair \
	--peOverlapNbasesMin 10 \
	--outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp_{params.sample} \
	--sjdbOverhang ${{readlength}}
rm -rf {params.sample}._STAR*
samtools flagstat {output.starbam} > {output.bamflagstat}
# The above STAR command will generate {sample}.Aligned.out.bam and {sample}.Aligned.toTranscriptome.out.bam
samtools view -@{threads} -b -F 4 -F 256 {output.starbam} > ${{tmpdir}}/{params.sample}.post_secondary_supplementary_filter.tmp.bam

# Sort the "post_secondary_supplementary_filter" BAM and collect some stats
samtools sort -@{threads} --output-fmt BAM -o ${{tmpdir}}/{params.sample}.post_secondary_supplementary_filter.bam ${{tmpdir}}/{params.sample}.post_secondary_supplementary_filter.tmp.bam && rm -f ${{tmpdir}}/{params.sample}.post_secondary_supplementary_filter.tmp.bam
samtools flagstat -@{threads} ${{tmpdir}}/{params.sample}.post_secondary_supplementary_filter.bam > {output.postsecondarysupplementaryfilterbamflagstat}
rm -f ${{tmpdir}}/{params.sample}.post_secondary_supplementary_filter.tmp.bam
samtools index -@{threads} ${{tmpdir}}/{params.sample}.post_secondary_supplementary_filter.bam

# Apply the "number of insertions in read" filter
python {params.filter_script} -i ${{tmpdir}}/{params.sample}.post_secondary_supplementary_filter.bam -o ${{tmpdir}}/{params.sample}.post_insertion_filter.bam -q 0 -n {params.ninsertionfilter}

# Sort, FixMateInfo, collect stats
java -Xmx{params.mem}g -jar $PICARDJARPATH/picard.jar FixMateInformation \
 I=${{tmpdir}}/{params.sample}.post_insertion_filter.bam \
 O=${{tmpdir}}/{params.sample}.postfixmate.bam \
 ASSUME_SORTED=false \
 CREATE_INDEX=true \
 SORT_ORDER=coordinate
samtools flagstat -@{threads} ${{tmpdir}}/{params.sample}.postfixmate.bam > {output.postinsertionfilterbamflagstat}

# Apply the MAPQ filter and remove "widowed" reads
python {params.filter_script} -i ${{tmpdir}}/{params.sample}.postfixmate.bam -o {output.bam} -q {params.mapqfilter} -n 1000

# Index and collect stats
samtools index -@{threads} {output.bam}
samtools flagstat -@{threads} {output.bam} > ${{tmpdir}}/{params.sample}

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

samtools view -@{threads} -b -f 128 -F 16 {output.bam} > ${{tmpdir}}/{params.sample}.F1.bam
samtools view -@{threads} -b -f 80 {output.bam} > ${{tmpdir}}/{params.sample}.F2.bam
samtools merge -c -f -p -@{threads} ${{tmpdir}}/{params.sample}.F.bam ${{tmpdir}}/{params.sample}.F1.bam ${{tmpdir}}/{params.sample}.F2.bam
#rm -f ${{tmpdir}}/{params.sample}.F1.bam ${{tmpdir}}/{params.sample}.F2.bam
samtools view -@{threads} -b -f 144 {output.bam} > ${{tmpdir}}/{params.sample}.R1.bam
samtools view -@{threads} -b -f 64 -F 16 {output.bam} > ${{tmpdir}}/{params.sample}.R2.bam
samtools merge -c -f -p -@{threads} ${{tmpdir}}/{params.sample}.R.bam ${{tmpdir}}/{params.sample}.R1.bam ${{tmpdir}}/{params.sample}.R2.bam
#rm -f ${{tmpdir}}/{params.sample}.R1.bam ${{tmpdir}}/{params.sample}.R2.bam

# Sort and collect stats for plus strand BAM
samtools sort -@56 --output-fmt BAM -o {output.plusbam} ${{tmpdir}}/{params.sample}.F.bam
samtools flagstat -@56 {output.plusbam} > {output.plusbamflagstat}

# Sort and collect stats for minus strand BAM
samtools sort -@56 --output-fmt BAM -o {output.minusbam} ${{tmpdir}}/{params.sample}.R.bam
samtools flagstat -@56 {output.minusbam} > {output.minusbamflagstat}

# cleanup
rm -rf ${{tmpdir}}/{params.sample}.*
rm -rf _STARtmp

"""

rule hisat:
	input:
		if1=rules.cutadapt.output.of1,
		if2=rules.cutadapt.output.of2
	output:
		bam=join(WORKDIR,"hisat2","{sample}.bam"),
		plusbam=join(WORKDIR,"hisat2","{sample}.plus.bam"),
		minusbam=join(WORKDIR,"hisat2","{sample}.minus.bam"),
		postsecondarysupplementaryfilterbamflagstat=join(WORKDIR,"hisat2","{sample}.post_secondary_supplementary_filter.bam.flagstat"),
		postinsertionfilterbamflagstat=join(WORKDIR,"hisat2","{sample}.post_insertion_filter.bam.flagstat"),
		hisat2bamflagstat=join(WORKDIR,"hisat2","{sample}.bam.flagstat"),
		plusbamflagstat=join(WORKDIR,"hisat2","{sample}.plus.bam.flagstat"),
		minusbamflagstat=join(WORKDIR,"hisat2","{sample}.minus.bam.flagstat"),
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
# Align with HISAT and remove secondary/supplementary alignments

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
sambamba flagstat --nthreads={threads} /dev/shm/{params.sample}.post_secondary_supplementary_filter.bam > {output.postsecondarysupplementaryfilterbamflagstat}

# Apply the "number of insertions in read" filter
python {params.filter_script} -i /dev/shm/{params.sample}.post_secondary_supplementary_filter.bam -o /dev/shm/{params.sample}.post_insertion_filter.tmp.bam -q 0 -n {params.ninsertionfilter}

# Sort, FixMateInfo, collect stats
sambamba sort -n --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out=/dev/shm/{params.sample}.post_insertion_filter.bam /dev/shm/{params.sample}.post_insertion_filter.tmp.bam && rm -f /dev/shm/{params.sample}.post_insertion_filter.tmp.bam
java -Xmx{params.mem}g -jar $PICARDJARPATH/picard.jar FixMateInformation I=/dev/shm/{params.sample}.post_insertion_filter.bam O=/dev/stdout ASSUME_SORTED=true QUIET=true \
 | sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out=/dev/shm/{params.sample}.postsambambasort.bam /dev/stdin
sambamba flagstat --nthreads={threads} /dev/shm/{params.sample}.postsambambasort.bam > {output.postinsertionfilterbamflagstat}

# Apply the MAPQ filter and remove "widowed" reads
python {params.filter_script} -i /dev/shm/{params.sample}.postsambambasort.bam -o {output.bam} -q {params.mapqfilter} -n 1000

# Index and collect stats
sambamba index --nthreads={threads} {output.bam}
sambamba flagstat --nthreads={threads} {output.bam} > {output.hisat2bamflagstat}

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
sambamba flagstat --nthreads={threads} {output.plusbam} > {output.plusbamflagstat}

# Sort and collect stats for minus strand BAM
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.minusbam} /dev/shm/{params.sample}.R.bam && rm -f /dev/shm/{params.sample}.R.bam
sambamba flagstat --nthreads={threads} {output.minusbam} > {output.minusbamflagstat}

# cleanup
rm -rf /dev/shm/{params.sample}.*

"""

rule create_toSNPcalling_BAM:
	input:
		plusbam=rules.hisat.output.plusbam,
		minusbam=rules.hisat.output.minusbam,
		readids=rules.get_fastuniq_readids.output.readids
	output:
		plustoSNPcallingbam=join(WORKDIR,"hisat2","{sample}.plus.toSNPcalling.bam"),
		minustoSNPcallingbam=join(WORKDIR,"hisat2","{sample}.minus.toSNPcalling.bam"),
		plustoSNPcallingbamflagstat=join(WORKDIR,"hisat2","{sample}.plus.toSNPcalling.bam.flagstat"),
		minustoSNPcallingbamflagstat=join(WORKDIR,"hisat2","{sample}.minus.toSNPcalling.bam.flagstat"),
	params:
		sample="{sample}",
		mem=MEMORY,
		workdir=WORKDIR,
		script=join(SCRIPTSDIR,"filter_bam_by_readids.py")
	envmodules: TOOLS["sambamba"]["version"]
	threads: 4
	shell:"""
python {params.script} -i {input.plusbam} -o /dev/shm/{params.sample}.plus.bam --readids {input.readids}
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.plustoSNPcallingbam} /dev/shm/{params.sample}.plus.bam && rm -f /dev/shm/{params.sample}.plus.bam
sambamba flagstat --nthreads={threads} {output.plustoSNPcallingbam} > {output.plustoSNPcallingbamflagstat}
python {params.script} -i {input.minusbam} -o /dev/shm/{params.sample}.minus.bam --readids {input.readids}
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.minustoSNPcallingbam} /dev/shm/{params.sample}.minus.bam && rm -f /dev/shm/{params.sample}.minus.bam
sambamba flagstat --nthreads={threads} {output.minustoSNPcallingbam} > {output.minustoSNPcallingbamflagstat}
"""



rule call_mutations:
	input:
		plusbam=rules.create_toSNPcalling_BAM.output.plustoSNPcallingbam,
		minusbam=rules.create_toSNPcalling_BAM.output.minustoSNPcallingbam
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
bcftools view -i '(REF=="T" & ALT="C")' --threads {threads} - > /dev/shm/{params.sample}.TtoC.bcf
bcftools sort -T /dev/shm /dev/shm/{params.sample}.TtoC.bcf | bgzip > {output.plusvcf}
tabix -p vcf {output.plusvcf}
# minus strand
bcftools mpileup -f {params.hisatindex}/{params.genome}/{params.genome}.fa -a AD,ADF,ADR {input.minusbam} | \
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
		mutatedbamflagstat=join(WORKDIR,"bams","{sample}.mutated.bam.flagstat"),
		unmutatedbamflagstat=join(WORKDIR,"bams","{sample}.unmutated.bam.flagstat"),
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
sort | uniq > /dev/shm/{params.sample}.plus.readids

minusvcf=$(basename {input.minusvcf})
zcat {input.minusvcf} > /dev/shm/${{minusvcf%.*}}
java -Xmx{params.mem}g -jar {params.sam2tsvjar} \
 --reference {params.hisatindex}/{params.genome}/{params.genome}.fa \
 --skip-N \
 {input.minusbam} | \
python {params.tsv2readidspy} /dev/shm/${{minusvcf%.*}} "A" "G" | \
sort | uniq > /dev/shm/{params.sample}.minus.readids

python {params.filterbyreadidspy} -i {input.plusbam} -o /dev/shm/{params.sample}.mutated.plus.bam --readids /dev/shm/{params.sample}.plus.readids -o2 /dev/shm/{params.sample}.unmutated.plus.bam
python {params.filterbyreadidspy} -i {input.minusbam} -o /dev/shm/{params.sample}.mutated.minus.bam --readids /dev/shm/{params.sample}.minus.readids -o2 /dev/shm/{params.sample}.unmutated.minus.bam

samtools merge -c -f -p -@ {threads} /dev/shm/{params.sample}.mutated.bam /dev/shm/{params.sample}.mutated.plus.bam /dev/shm/{params.sample}.mutated.minus.bam && rm -f /dev/shm/{params.sample}.mutated.plus.bam /dev/shm/{params.sample}.mutated.minus.bam
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.mutatedbam} /dev/shm/{params.sample}.mutated.bam && rm -f /dev/shm/{params.sample}.mutated.bam
sambamba flagstat --nthreads={threads} {output.mutatedbam} > {output.mutatedbamflagstat}

samtools merge -c -f -p -@ {threads} /dev/shm/{params.sample}.unmutated.bam /dev/shm/{params.sample}.unmutated.plus.bam /dev/shm/{params.sample}.unmutated.minus.bam && rm -f /dev/shm/{params.sample}.unmutated.plus.bam /dev/shm/{params.sample}.unmutated.minus.bam
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.unmutatedbam} /dev/shm/{params.sample}.unmutated.bam && rm -f /dev/shm/{params.sample}.unmutated.bam
sambamba flagstat --nthreads={threads} {output.unmutatedbam} > {output.unmutatedbamflagstat}
"""

rule gtf2bed12:
	input:
		gtf=GTF
	output:
		bed12=join(WORKDIR,GENOME+".bed12")
	envmodules: TOOLS["ucsc"]["version"]
	shell:"""
bn=$(basename {input.gtf})
gtfToGenePred {input.gtf} /dev/shm/${{bn%.*}}.genepred
genePredToBed /dev/shm/${{bn%.*}}.genepred {output.bed12}
"""

rule htseq:
	input:
		mutatedbam=rules.split_bam_by_mutation.output.mutatedbam,
		unmutatedbam=rules.split_bam_by_mutation.output.unmutatedbam,
		bed12=rules.gtf2bed12.output.bed12
	output:
		mutatedcounts_exon=join(WORKDIR,"htseq","{sample}.mutated.exon.counts"),
		mutatedcounts_gene=join(WORKDIR,"htseq","{sample}.mutated.gene.counts"),
		unmutatedcounts_exon=join(WORKDIR,"htseq","{sample}.unmutated.exon.counts"),
		unmutatedcounts_gene=join(WORKDIR,"htseq","{sample}.unmutated.gene.counts"),
	params:
		gtf=GTF,
		sample="{sample}",
		htseqparams=config["htseqparams"]
	envmodules: TOOLS["htseq"]["version"], TOOLS["rseqc"]["version"]
	shell:"""
infer_experiment.py -i {input.mutatedbam} -r {input.bed12} -s 500000 > /dev/shm/{params.sample}.strand.info
strandedness=$(tail -n1 /dev/shm/{params.sample}.strand.info|awk '{{if ( $NF > 0.75 ) {{print "reverse"}} else if ( $NF < 0.25 ) {{print "yes"}} else {{print "no"}}}}')
htseq-count -t "exon" -f bam -r pos --stranded $strandedness --mode intersection-nonempty --additional-attr=gene_name {params.htseqparams} {input.mutatedbam} {params.gtf} > {output.mutatedcounts_exon}
htseq-count -t "exon" -f bam -r pos --stranded $strandedness --mode intersection-nonempty --additional-attr=gene_name {params.htseqparams} {input.unmutatedbam} {params.gtf} > {output.unmutatedcounts_exon}
htseq-count -t "gene" -f bam -r pos --stranded $strandedness --mode intersection-nonempty --additional-attr=gene_name {params.htseqparams} {input.mutatedbam} {params.gtf} > {output.mutatedcounts_gene}
htseq-count -t "gene" -f bam -r pos --stranded $strandedness --mode intersection-nonempty --additional-attr=gene_name {params.htseqparams} {input.unmutatedbam} {params.gtf} > {output.unmutatedcounts_gene}
"""

rule get_split_reads:
	input:
		unpack(get_fastqs),
		bam=join(WORKDIR,"bams","{sample}.{mutated}.bam"),
	output:
		R1=join(WORKDIR,"filtered_reads","{sample}.{mutated}.R1.fastq.gz"),
		R2=join(WORKDIR,"filtered_reads","{sample}.{mutated}.R2.fastq.gz"),
	params:
		sample="{sample}",
		mutated="{mutated}",
		pyscript=join(SCRIPTSDIR,"filter_fastq_by_readids_highmem.py"),
		fastq_pair=join(RESOURCESDIR,"fastq_pair")
	envmodules: TOOLS["samtools"]["version"]
	shell:"""
samtools view {input.bam}|cut -f1|sort|uniq > /dev/shm/{params.sample}.{params.mutated}.readids
R1fastq=$(basename {output.R1})
R1fastq="/dev/shm/${{R1fastq%.*}}"
R2fastq=$(basename {output.R2})
R2fastq="/dev/shm/${{R2fastq%.*}}"
python {params.pyscript} \
	--infq {input.R1} \
	--outfq $R1fastq \
	--readids /dev/shm/{params.sample}.{params.mutated}.readids
python {params.pyscript} \
	--infq {input.R2} \
	--outfq $R2fastq \
	--readids /dev/shm/{params.sample}.{params.mutated}.readids
cd /dev/shm
{params.fastq_pair} $R1fastq $R2fastq
pigz -p4 ${{R1fastq}}.paired.fq
pigz -p4 ${{R2fastq}}.paired.fq
rsync -az --progress ${{R1fastq}}.paired.fq.gz {output.R1}
rsync -az --progress ${{R2fastq}}.paired.fq.gz {output.R2}
rm -rf /dev/shm/*.fastq /dev/shm/*.fq /dev/shm/*.fq.gz /dev/shm/*.fastq.gz /dev/shm/*.readids
"""

rule get_nfragments_json:
	input:
		rules.get_fastq_nreads.output.o1,
		rules.get_fastq_nreads.output.o2,
		rules.get_fastq_nreads.output.o3,
		rules.hisat.output.postsecondarysupplementaryfilterbamflagstat,
		rules.hisat.output.postinsertionfilterbamflagstat,
		rules.hisat.output.hisat2bamflagstat,
		rules.hisat.output.plusbamflagstat,
		rules.hisat.output.minusbamflagstat,
		rules.split_bam_by_mutation.output.mutatedbamflagstat,
		rules.split_bam_by_mutation.output.unmutatedbamflagstat,
		rules.call_mutations.output.vcf,
		rules.call_mutations.output.plusvcf,
		rules.call_mutations.output.minusvcf,
	output:
		json=join(WORKDIR,"qc","nfragments","{sample}.json")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		outdir=join(WORKDIR,"qc","nfragments"),
		pyscript=join(SCRIPTSDIR,"get_per_sample_nfragments.py")
	shell:"""
cd {params.workdir}
python {params.pyscript} {params.sample} {output.json}
"""

rule create_nfragments_table:
	input:
		expand(join(WORKDIR,"qc","nfragments","{sample}.json"),sample= SAMPLES)
	output:
		table=join(WORKDIR,"qc","nfragments","nfragments.tsv")
	params:
		workdir=WORKDIR,
		outdir=join(WORKDIR,"qc","nfragments"),
		pyscript=join(SCRIPTSDIR,"nfragments_json2table.py")
	shell:"""
cd {params.outdir}
python {params.pyscript} {output.table}
"""


