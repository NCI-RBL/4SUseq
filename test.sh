#!/usr/bin/env bash
set -euo pipefail

# ## setting PIPELINE_HOME
# ## clone the pipeline to a folder
# ## git clone https://github.com/RBL-NCI/4SUseq.git
# ## and set that as the pipeline home
PIPELINE_HOME="/data/RBL_NCI/Wolin/mESC_slam_analysis/workflow/4SUseq"
echo "Pipeline Dir: $PIPELINE_HOME"
# hg38 and mm10 hisat2 prebuilt indices
HISAT2INDEXDIR="/data/RBL_NCI/Wolin/mESC_slam_analysis/resource"
# ## make current folder as working directory
a=$(readlink -f $0)
# echo $a
for i in `seq 1 20`;do
b=$(echo $a|sed "s/\/gpfs\/gsfs${i}\/users/\/data/g")
a=$b
done
# echo $a
WORKDIR=$(dirname $a)

function usage() { cat << EOF
test.sh: test the workflow.
USAGE:
  bash test.sh <MODE> 
Required Positional Argument:
  MODE: [Type: Str] Valid options:
    a) init: initial workdir
    b) inittest: initial workdir for testing 
    c) reset: cleanup followed by init
    d) cleanup: delete folders to get ready for re-init
    e) unlock: snakemake --unlock
    f) dryrun: snakemake --dry-run
    g) run: run test data
    h) runlocal: run without submitting to sbatch
EOF
}

function err() { cat <<< "$@" 1>&2; }

function init() {
echo "Working Dir: $WORKDIR"
cd $WORKDIR
WORKDIR_esc=${WORKDIR////\\\/}
PIPELINE_HOME_esc=${PIPELINE_HOME////\\\/}
RESOURCESDIR="${WORKDIR}/resources"
RESOURCESDIR_esc=${RESOURCESDIR////\\\/}
HISAT2INDEXDIR_esc=${HISAT2INDEXDIR////\\\/}
# make sure that config folder exists in the current folder
if [ ! -d $WORKDIR/config ]; then cp -r $PIPELINE_HOME/config $WORKDIR/;echo "Config Dir: $WORKDIR/config";fi
# change the working directory to current directory
# sed -i "s/WORKDIR/\/data\/RBL_NCI\/Wolin\/mESC_slam_analysis\/test_workdir/g" $WORKDIR/config/config.yaml
sed -i "s/WORKDIR/$WORKDIR_esc/g" $WORKDIR/config/config.yaml
sed -i "s/WORKDIR/$WORKDIR_esc/g" $WORKDIR/config/config_test.yaml
sed -i "s/HISAT2INDEXDIR/$HISAT2INDEXDIR_esc/g" $WORKDIR/config/config.yaml
sed -i "s/HISAT2INDEXDIR/$HISAT2INDEXDIR_esc/g" $WORKDIR/config/config_test.yaml
sed -i "s/RESOURCESDIR/$RESOURCESDIR_esc/g" $WORKDIR/config/config.yaml
sed -i "s/RESOURCESDIR/$RESOURCESDIR_esc/g" $WORKDIR/config/config_test.yaml
sed -i "s/PIPELINE_HOME/$PIPELINE_HOME_esc/g" $WORKDIR/config/samples_test.tsv
if [ ! -d $WORKDIR/scripts ]; then cp -r $PIPELINE_HOME/workflow/scripts $WORKDIR/;echo "Scripts Dir: $WORKDIR/scripts";fi
if [ ! -d $WORKDIR/resources ]; then cp -r $PIPELINE_HOME/resources $WORKDIR/;echo "Resources Dir: $WORKDIR/resources";fi

#create log and stats folders
if [ ! -d $WORKDIR/logs ]; then mkdir -p $WORKDIR/logs;echo "Logs Dir: $WORKDIR/logs";fi
if [ ! -d $WORKDIR/stats ];then mkdir -p $WORKDIR/stats;echo "Stats Dir: $WORKDIR/stats";fi

}

function inittest {
  init
  if [ -f $WORKDIR/config/config.yaml.original ];then rm -f $WORKDIR/config/config.yaml.original;fi
  if [ -f $WORKDIR/config/cluster.json.original ];then rm -f $WORKDIR/config/cluster.json.original;fi
  cp $WORKDIR/config/config.yaml $WORKDIR/config/config.yaml.original
  cp $WORKDIR/config/cluster.json $WORKDIR/config/cluster.json.original
  cp $WORKDIR/config/config_test.yaml $WORKDIR/config/config.yaml
  cp $WORKDIR/config/cluster_test.json $WORKDIR/config/cluster.json
}

function run () {
  ## initialize if not already done
  echo "Working Dir: $WORKDIR"
  if [ ! -d $WORKDIR/config ]; then err "Error: config folder not found ... initialize first!";usage && exit 1;fi
  for f in config.yaml tools.yaml cluster.json samples.tsv; do
    if [ ! -f $WORKDIR/config/$f ]; then err "Error: '${f}' file not found in config folder ... initialize first!";usage && exit 1;fi
  done
  ## Archive previous run files
  if [ -f ${WORKDIR}/snakemake.log ];then 
    modtime=$(stat ${WORKDIR}/snakemake.log |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
    mv ${WORKDIR}/snakemake.log ${WORKDIR}/stats/snakemake.log.${modtime} && gzip -n ${WORKDIR}/stats/snakemake.log.${modtime}
  fi
  if [ -f ${WORKDIR}/snakemake.stats ];then 
    modtime=$(stat ${WORKDIR}/snakemake.stats |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
    mv ${WORKDIR}/snakemake.stats ${WORKDIR}/stats/snakemake.stats.${modtime} && gzip -n ${WORKDIR}/stats/snakemake.stats.${modtime}
  fi


  module load python/3.7
  module load snakemake/5.24.1
  # module load singularity

  # --use-singularity \
  # --singularity-args "-B ${WORKDIR}" \


  if [ "$1" == "local" ];then

  snakemake -s ${PIPELINE_HOME}/workflow/Snakefile \
  --directory $WORKDIR \
  --use-envmodules \
  --printshellcmds \
  --latency-wait 120 \
  --configfile $WORKDIR/config/config.yaml \
  --cores all \
  --stats ${WORKDIR}/snakemake.stats \
  2>&1|tee ${WORKDIR}/snakemake.log

  else

  snakemake $1 -s ${PIPELINE_HOME}/workflow/Snakefile \
  --directory $WORKDIR \
  --use-envmodules \
  --printshellcmds \
  --latency-wait 120 \
  --configfile $WORKDIR/config/config.yaml \
  --cluster-config $WORKDIR/config/cluster.json \
  --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}" \
  -j 500 \
  --rerun-incomplete \
  --keep-going \
  --stats ${WORKDIR}/snakemake.stats \
  2>&1|tee ${WORKDIR}/snakemake.log

  fi

  mv slurm*out stats && for a in $(ls stats/slurm*out);do gzip -n $a;done
}

function cleanup() {
  rm -rf ${WORKDIR}/config
  rm -rf ${WORKDIR}/logs
  rm -rf ${WORKDIR}/stats
  rm -rf ${WORKDIR}/scripts
  rm -rf ${WORKDIR}/resources
  rm -rf *snakemake*
}


function main(){

  if [ $# -eq 0 ]; then usage; exit 1; fi

  case $1 in
    init) init && exit 0;;
    inittest) inittest && exit 0;;
    dryrun) run --dry-run && exit 0;;
    unlock) run --unlock && exit 0;;
    run) run "" && exit 0;;
    runlocal) run local && exit 0;;
    cleanup) cleanup && exit 0;;
    reset) cleanup && init && exit 0;;
    -h | --help | help) usage && exit 0;;
    -* | --*) err "Error: Failed to provide mode: <init|run>."; usage && exit 1;;
    *) err "Error: Failed to provide mode: <init|run>. '${1}' is not supported."; usage && exit 1;;
  esac
}

main "$@"





