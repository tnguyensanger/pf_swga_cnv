#!/bin/bash

source activate pf_swga_cnv_analysis_py35_linux

CWD=$(dirname "$(readlink -f "$0")")

snakemake --snakefile $CWD/xhmm.snakefile --cluster-config $CWD/cluster_config.yml --cluster "bsub -n {cluster.ncpus} -P {cluster.project} -q {cluster.queue} -J {cluster.jobname} -M {cluster.memory} -R '{cluster.resources}' -e {cluster.error} -o {cluster.output}" --jobs 1
