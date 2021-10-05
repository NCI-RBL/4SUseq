
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
    threads: getthreads("fastq_screen")
    envmodules: TOOLS["fastq_screen"]["version"], TOOLS["bowtie"]["version"]
    shell:"""
set -e -x -o pipefail
fastq_screen --conf {params.conf} \
 --outdir "{params.outdir}" \
 --threads {threads} --subset 1000000 \
 --aligner bowtie2 --force \
 {input.if1} \
 {input.if2}
"""	

rule bbtools:
    input:
        if1=rules.cutadapt.output.of1,
        if2=rules.cutadapt.output.of2
    output:
        hist=join(WORKDIR,"qc","insert_size_hist","{sample}.bbtools.hist.txt")
    params:
        sample="{sample}",
        workdir=WORKDIR,
        outdir=join(WORKDIR,"qc","fastqscreen"),
        mem=getmemg("bbtools")
    threads: getthreads("fastq_screen")
    envmodules: TOOLS["bbtools"]["version"],
    shell:"""
set -e -x -o pipefail
bbtools bbmerge-auto \
in1={input.if1} \
in2={input.if2} \
ihist={output.hist} \
k=62 extend2=200 rem ecct -Xmx{params.mem}
"""

rule get_nfragments_json:
    input:
        rules.get_fastq_nreads.output.o1,
        rules.get_fastq_nreads.output.o2,
        rules.get_fastq_nreads.output.o3,
        rules.star.output.postsecondarysupplementaryfilterbamflagstat,
        rules.star.output.postinsertionfilterbamflagstat,
        rules.star.output.bamflagstat,
        rules.star.output.plusbamflagstat,
        rules.star.output.minusbamflagstat,
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
set -e -x -o pipefail
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
set -e -x -o pipefail
cd {params.outdir}
python {params.pyscript} {output.table}
"""

rule qualimap:
    input:
        bam=rules.star.output.bam
        strand_info=rules.strand_info.output.strand_info
    output:
        html=join(WORKDIR,"qc","qualimap","{sample}","report.html")
    params:
        sample="{sample}",
        outdir=join(WORKDIR,"qc","qualimap","{sample}"),
        mem=getmemG("qualimap")
    threads: getthreads("qualimap")
    shell:"""
set -e -x -o
qualimap bamqc \
 --java-mem-size={params.mem} \
 -bam {input.bam} \
 -c -nw 2800 -ip -gd mm10 -hm 3 \
 -outdir {params.outdir} 
"""