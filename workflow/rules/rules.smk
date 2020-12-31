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

rule fastuniq:
	input:
		if1=rules.cutadapt.output.of1,
		if2=rules.cutadapt.output.of2
	output:
		of1=join(FASTUNIQDIR,"{sample}.R1.trim.fastuniq.fastq.gz"),
		of2=join(FASTUNIQDIR,"{sample}.R2.trim.fastuniq.fastq.gz")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		outdir=FASTUNIQDIR,
		fastuniq=join(RESOURCESDIR,"fastuniq")
	envmodules: TOOLS["pigz"]["version"]
	threads: 56
	shell:"""
zcat {input.if1} > /dev/shm/{params.sample}.R1.fastq
zcat {input.if2} > /dev/shm/{params.sample}.R2.fastq
echo -ne "/dev/shm/{params.sample}.R1.fastq\n/dev/shm/{params.sample}.R2.fastq\n" > /dev/shm/{params.sample}.list
o1={output.of1}
o2={output.of2}
{params.fastuniq} -i /dev/shm/{params.sample}.list -t q -o ${{o1%.*}} -p ${{o2%.*}} -c 0
pigz -p4 -f ${{o1%.*}}
pigz -p4 -f ${{o2%.*}}
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
bash {params.bashscript} {params.scriptsdir}
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

rule hisat_on_fastuniq:
	input:
		if1=rules.fastuniq.output.of1,
		if2=rules.fastuniq.output.of2
	output:
		bam=join(WORKDIR,"hisat2","{sample}.toSNPcalling.bam"),
		flagstat1=join(WORKDIR,"hisat2","{sample}.post_secondary_supplementary_filter.bam.flagstat"),
		flagstat2=join(WORKDIR,"hisat2","{sample}.post_insertion_filter.bam.flagstat"),
		flagstat3=join(WORKDIR,"hisat2","{sample}.toSNPcalling.bam.flagstat")
	params:
		sample="{sample}",
		mem=MEMORY,
		workdir=WORKDIR,
		outdir=join(WORKDIR,"hisat2"),
		genome=config["genome"],
		hisatindex=config["hisatindex"],
		splicesites=join(RESOURCESDIR,config["genome"]+".splicesites.txt"),
		filter_script=join(SCRIPTSDIR,"filter_bam.py"),
		mapqfilter=config['mapqfilter'],
		ninsertionfilter=config['ninsertionfilter']
	threads: 56
	envmodules: TOOLS["hisat"]["version"], TOOLS["samtools"]["version"], TOOLS["sambamba"]["version"],  TOOLS["picard"]["version"],  TOOLS["bbtools"]["version"]
	shell:"""
hisat2 \
 -x {params.hisatindex}/{params.genome}/{params.genome} \
 -1 {input.if1} \
 -2 {input.if2} \
 --rna-strandness RF \
 --known-splicesite-infile {params.splicesites} \
 --summary-file {params.outdir}/{params.sample}.hisat2.summary.log \
 --threads {threads} \
 --rg-id {params.sample} --rg SM:{params.sample} \
 --un-gz {params.outdir}/{params.sample}.unmapped.fastq.gz \
 --mp 4,2 \
 | samtools view -@{threads} -bS -F 4 -F 256 - > /dev/shm/{params.sample}.post_secondary_supplementary_filter.tmp.bam
sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out=/dev/shm/{params.sample}.post_secondary_supplementary_filter.bam /dev/shm/{params.sample}.post_secondary_supplementary_filter.tmp.bam && rm -f /dev/shm/{params.sample}.post_secondary_supplementary_filter.tmp.bam
sambamba flagstat --nthreads={threads} /dev/shm/{params.sample}.post_secondary_supplementary_filter.bam > {output.flagstat1}
python {params.filter_script} -i /dev/shm/{params.sample}.post_secondary_supplementary_filter.bam -o /dev/shm/{params.sample}.post_insertion_filter.tmp.bam -q 0 -n {params.ninsertionfilter}
sambamba sort -n --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out=/dev/shm/{params.sample}.post_insertion_filter.bam /dev/shm/{params.sample}.post_insertion_filter.tmp.bam && rm -f /dev/shm/{params.sample}.post_insertion_filter.tmp.bam
java -Xmx{params.mem}g -jar $PICARDJARPATH/picard.jar FixMateInformation I=/dev/shm/{params.sample}.post_insertion_filter.bam O=/dev/stdout ASSUME_SORTED=true QUIET=true \
 | sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out=/dev/shm/{params.sample}.postsambambasort.bam /dev/stdin
sambamba flagstat --nthreads={threads} /dev/shm/{params.sample}.postsambambasort.bam > {output.flagstat2}
python {params.filter_script} -i /dev/shm/{params.sample}.postsambambasort.bam -o {output.bam} -q {params.mapqfilter} -n 1000
sambamba index --nthreads={threads} {output.bam}
sambamba flagstat --nthreads={threads} {output.bam} > {output.flagstat3}
"""

rule call_mutations:
	input:
		bam=rules.hisat_on_fastuniq.output.bam
	output:
		genicbcf=join(WORKDIR,"vcf","{sample}.genic.bcf"),
		TtoCvcf=join(WORKDIR,"vcf","{sample}.TtoC.vcf")
	params:
		sample="{sample}",
		workdir=WORKDIR,
		outdir=join(WORKDIR,"vcf"),
		genome=config["genome"],
		hisatindex=config["hisatindex"],
		tab=join(RESOURCESDIR,config["genome"]+".genes.tab.gz"),
		hdr=join(RESOURCESDIR,"hdr.txt"),
		filter_script=join(SCRIPTSDIR,"filter_bam.py")
	threads: 56
	envmodules: TOOLS["bcftools"]["version"]
	shell:"""
bcftools mpileup -f {params.hisatindex}/{params.genome}/{params.genome}.fa -a AD,ADF,ADR {input.bam} | \
bcftools call -mv -Ob --threads {threads} | \
bcftools filter -i '%QUAL>45' -g3 -Ob --threads {threads} - | \
bcftools annotate -a {params.tab} -h {params.hdr} -c CHROM,FROM,TO,strand,gene_id,gene_type,gene_name -Ob --threads {threads} - | \
bcftools view -i 'INFO/strand=="-" | INFO/strand=="+"' -Ob --threads {threads} - > {output.genicbcf}
bcftools view -i '(REF=="A" & ALT="G") | (REF=="T" & ALT=="C")' --threads {threads} {output.genicbcf} > {output.TtoCvcf}
"""

rule hisat_on_cutadapt:
	input:
		if1=rules.cutadapt.output.of1,
		if2=rules.cutadapt.output.of2
	output:
		bam=join(WORKDIR,"hisat2_allreads","{sample}.bam"),
		flagstat=join(WORKDIR,"hisat2_allreads","{sample}.bam.flagstat")
	params:
		sample="{sample}",
		mem=MEMORY,
		workdir=WORKDIR,
		outdir=join(WORKDIR,"hisat2_allreads"),
		genome=config["genome"],
		hisatindex=config["hisatindex"],
		splicesites=join(RESOURCESDIR,config["genome"]+".splicesites.txt"),
		filter_script=join(SCRIPTSDIR,"filter_bam.py")
	threads: 56
	envmodules: TOOLS["hisat"]["version"], TOOLS["samtools"]["version"], TOOLS["sambamba"]["version"],  TOOLS["picard"]["version"],  TOOLS["bbtools"]["version"]
	shell:"""
hisat2 \
 -x {params.hisatindex}/{params.genome}/{params.genome} \
 -1 {input.if1} \
 -2 {input.if2} \
 --rna-strandness RF \
 --known-splicesite-infile {params.splicesites} \
 --summary-file {params.outdir}/{params.sample}.hisat2.summary.log \
 --threads {threads} \
 --rg-id {params.sample} --rg SM:{params.sample} \
 --un-gz {params.outdir}/{params.sample}.unmapped.fastq.gz \
 --mp 4,2 \
 | samtools view -@{threads} -bS -F 4 -F 256 - \
 | sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.bam} /dev/stdin
sambamba flagstat --nthreads={threads} {output.bam} > {output.flagstat}
"""

rule get_nonmutant_reads:
	input:
		bam=rules.hisat_on_cutadapt.output.bam,
		vcf=rules.call_mutations.output.TtoCvcf,
		infq1=rules.cutadapt.output.of1,
		infq2=rules.cutadapt.output.of2
	output:
		outfq1=join(WORKDIR,"filtered_reads","{sample}.filtered.R1.fastq.gz"),
		outfq2=join(WORKDIR,"filtered_reads","{sample}.filtered.R2.fastq.gz"),
	params:
		sample="{sample}",
		mem=MEMORY,
		workdir=WORKDIR,
		outdir=join(WORKDIR,"filtered_reads"),
		genome=config["genome"],
		hisatindex=config["hisatindex"],
		tsv2readidspy=join(SCRIPTSDIR,"tsv2readids.py"),
		filterbyreadidspy=join(SCRIPTSDIR,"filter_fastq_by_readids.py"),
		sam2tsvjar=join(RESOURCESDIR,"sam2tsv.jar")
	envmodules: TOOLS["java"]["version"]
	shell:"""
java -Xmx{params.mem}g -jar {params.sam2tsvjar} \
 --reference {params.hisatindex}/{params.genome}/{params.genome}.fa \
 --skip-N \
 {input.bam} | \
python {params.tsv2readidspy} {input.vcf} | \
sort | uniq > /dev/shm/{sample}.readids
o1={output.outfq1}
o2={output.outfq2}
python {params.filterbyreadidspy} --infq {input.infq1} --outfq ${{o1%.*}} --readids /dev/shm/{sample}.readids --complement
python {params.filterbyreadidspy} --infq {input.infq2} --outfq ${{o2%.*}} --readids /dev/shm/{sample}.readids --complement
pigz -p4 ${{o1%.*}}
pigz -p4 ${{o2%.*}}
"""

rule get_nfragments_json:
	input:
		lanes1=rules.get_fastq_nreads.output.o1,
		lanes2=rules.get_fastq_nreads.output.o2,
		lanes3=rules.get_fastq_nreads.output.o3,
		flagstat1=rules.hisat_on_fastuniq.output.flagstat1,
		flagstat2=rules.hisat_on_fastuniq.output.flagstat2,
		flagstat3=rules.hisat_on_fastuniq.output.flagstat3,
		vcf=rules.call_mutations.output.TtoCvcf
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


