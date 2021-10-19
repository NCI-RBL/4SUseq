rule create_rsem_index:
    input:
        reffa=REFFA,
        gtf=GTF
    output:
        ump=join(RSEMINDEXDIR,"ref.transcripts.ump")
    params:
        rsemindexdir=RSEMINDEXDIR,
    envmodules: TOOLS['rsem']['version']
    threads: getthreads("create_rsem_index")
    shell:"""
set -euf -o pipefail
cd {params.rsemindexdir}
rsem-prepare-reference -p {threads} --gtf {input.gtf} {input.reffa} ref
rsem-generate-ngvector ref.transcripts.fa ref.transcripts
"""

rule create_bed12:
    input:
        gtf=GTF
    output:
        bed12=join(RSEMINDEXDIR,"ref.bed12")
    params:
        rsemindexdir=RSEMINDEXDIR
    envmodules:
        TOOLS['ucsc']['version'], 
    shell:"""
set -euf -o pipefail
cd {params.rsemindexdir}
genesgtf={input.gtf}
bn=$(basename $genesgtf)
gtfToGenePred -genePredExt -geneNameAsName2 $genesgtf ${{bn}}.tmp
awk -v OFS="\\t" '{{print $2,$4,$5,$1,"0",$3,$6,$7,"0",$8,$9,$10}}' ${{bn}}.tmp > {output.bed12}
rm -f ${{bn}}.tmp
"""

rule get_rsem_counts:
    input:
        starbam=rules.star.output.starbam, # used only for strandinfo
        tbam=rules.star.output.tbam,
        mutatedtbam=rules.split_tbam_by_mutation.output.mutatedtbam,
        unmutatedtbam=rules.split_tbam_by_mutation.output.unmutatedtbam,      
    output:
        strandinfo=join(WORKDIR,"rsem","{sample}","{sample}.strandinfo"),
        gcounts=join(WORKDIR,"rsem","{sample}","{sample}.RSEM.genes.results"),
        tcounts=join(WORKDIR,"rsem","{sample}","{sample}.RSEM.isoforms.results"),
        mutatedgcounts=join(WORKDIR,"rsem","{sample}","{sample}_mutated.RSEM.genes.results"),
        mutatedtcounts=join(WORKDIR,"rsem","{sample}","{sample}_mutated.RSEM.isoforms.results"),
        unmutatedgcounts=join(WORKDIR,"rsem","{sample}","{sample}_unmutated.RSEM.genes.results"),
        unmutatedtcounts=join(WORKDIR,"rsem","{sample}","{sample}_unmutated.RSEM.isoforms.results"),
    envmodules: TOOLS['rseqc']['version'], TOOLS['rsem']['version'], TOOLS["samtools"]["version"]
    threads: getthreads("get_rsem_counts")
    params:
        sample="{sample}",
        rsemdir=join(WORKDIR,"rsem")
    shell:"""
set -exuf -o pipefail
cd {params.rsemdir}
# inter strandedness
samtools index -@{threads} {input.starbam}
infer_experiment.py -r {input.bed12} -i {input.starbam} -s 1000000 > {output.strandinfo}
# Get strandedness to calculate Forward Probability
fp=$(tail -n1 {output.strandinfo} | awk '{{if($NF > 0.75) print "0.0"; else if ($NF<0.25) print "1.0"; else print "0.5";}}')
echo "Forward Probability Passed to RSEM: $fp"
rsemindex=$(echo {input.bed12}|sed "s@.bed12@@g")

rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  \
        --bam --paired-end -p {threads}  {input.tbam} $rsemindex {params.sample} --time \
        --temporary-folder /lscratch/$SLURM_JOBID --keep-intermediate-files --forward-prob=${{fp}} --estimate-rspd
mv {params.sample}.genes.results {output.gcounts}
mv {params.sample}.isoforms.results {output.tcounts}

rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  \
        --bam --paired-end -p {threads}  {input.mutatedtbam} $rsemindex {params.sample}_mutated --time \
        --temporary-folder /lscratch/$SLURM_JOBID --keep-intermediate-files --forward-prob=${{fp}} --estimate-rspd
mv {params.sample}_mutated.genes.results {output.mutatedgcounts}
mv {params.sample}_mutated.isoforms.results {output.mutatedtcounts}

rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  \
        --bam --paired-end -p {threads}  {input.unmutatedtbam} $rsemindex {params.sample}_unmutated --time \
        --temporary-folder /lscratch/$SLURM_JOBID --keep-intermediate-files --forward-prob=${{fp}} --estimate-rspd
mv {params.sample}_unmutated.genes.results {output.unmutatedgcounts}
mv {params.sample}_unmutated.isoforms.results {output.unmutatedtcounts}

"""