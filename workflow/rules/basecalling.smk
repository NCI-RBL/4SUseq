def get_fast5_path(wildcards):
    return SAMPLESDF["path_to_fast5_parent_folder"][wildcards.sample]


rule guppy:
    input:
        fast5path=get_fast5_path
    output:
        outfastq=join(WORKDIR,"fastqs","{sample}.fastq.gz"),
        sequencing_summary=join(WORKDIR,"fastqs","{sample}.sequencing_summary.txt.gz")
    params:
        flowcell = config["flowcell"],
        kit = config["kit"],
        sample = "{sample}"
    envmodules: TOOLS["guppy"]["version"], TOOLS["pigz"]["version"]
    log: join(WORKDIR,"logs","{sample}.guppy.log")
    shell:"""
guppy_basecaller --print_workflows 2>&1 | tee -a {log}
guppy_basecaller \
 --input_path {input.fast5path} \
 --recursive \
 --flowcell {params.flowcell} \
 --kit {params.kit} \
 -x cuda:all \
 --records_per_fastq 0 \
 --qscore_filtering \
 --save_path /lscratch/$SLURM_JOBID/{params.sample} 2>&1 |tee -a {log}
find /lscratch/$SLURM_JOBID/{params.sample} -name "*.fastq" -exec cat {{}} \; \
 | gzip -n -> {output.outfastq} 2>&1 |tee -a {log}
pigz -p $(nproc) /lscratch/$SLURM_JOBID/{params.sample}/sequencing_summary.txt && cp /lscratch/$SLURM_JOBID/{params.sample}/sequencing_summary.txt.gz {output.sequencing_summary} 2>&1 |tee -a {log}
"""
## Files created by guppy look like this:
# -rw-r--r-- 1 kopardevn CCBR 3.0K Nov 19 14:20 sequencing_summary.txt
# -rw-r--r-- 1 kopardevn CCBR 3.7K Nov 19 14:20 fastq_runid_d531634aaf7cba4fd8f7a1fba1d8ed9f0f81be2a_0_0.fastq
# -rw-r--r-- 1 kopardevn CCBR 3.3K Nov 19 14:20 fastq_runid_84d34f66eed213a95bd9b6aff2d24aac498555ff_0_0.fastq
# -rw-r--r-- 1 kopardevn CCBR  11K Nov 19 14:20 fastq_runid_013ea2ec6aebadbd80826ad673b152e04460f452_0_0.fastq
# -rw-r--r-- 1 kopardevn CCBR  64K Nov 19 14:20 sequencing_telemetry.js
# -rw-r--r-- 1 kopardevn CCBR 7.4K Nov 19 14:20 guppy_basecaller_log-2020-11-19_14-20-41.log

rule fastqc:
    input:
        expand(join(WORKDIR,"fastqs","{sample}.fastq.gz"),sample=SAMPLES)
    output:
        expand(join(WORKDIR,"qc","fastqc","{sample}_fastqc.zip"),sample=SAMPLES)
    params:
        outdir=join(WORKDIR,"qc","fastqc")
    threads: 16
    # envmodules: TOOLS["fastqc"]["version"]
    container: "nciccbr/ccbr_fastqc_v0.11.9:latest"
    log: join(WORKDIR,"logs","fastqc.log")
    shell:"""
fastqc {input} -t {threads} -o {params.outdir}
"""

rule pycoqc:
    input:
        # sequencing_summary=join(WORKDIR,"fastqs","{sample}.sequencing_summary.txt")
        sequencing_summary=rules.guppy.output.sequencing_summary
    output:
        pycoQChtml=join(WORKDIR,"qc","pycoQC","{sample}.pycoQC.html"),
        pycoQCjson=join(WORKDIR,"qc","pycoQC","{sample}.pycoQC.json"),
    params:
        outdir=join(WORKDIR,"qc","fastqc")
    # conda: "envs/pycoqc.yaml"
    container: "docker://nciccbr/ccbr_pycoqc_v2.5.0.23:latest"
    log: join(WORKDIR,"logs","{sample}.pycoQC.log")
    shell:"""
pycoQC -f {input.sequencing_summary} -o {output.pycoQChtml} -j {output.pycoQCjson} 2>&1 |tee -a {log}
"""