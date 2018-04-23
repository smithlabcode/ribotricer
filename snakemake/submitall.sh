#!/bin/bash
snakemake --snakefile Snakefile\
    --config config_path=$1\
    --js $PWD/jobscript.sh\
    --printshellcmds\
    --cluster-config $PWD/cluster.yaml\
    --jobname "{rulename}.{jobid}.$1"\
    --keep-going\
    --stats $PWD/$1.riboraptor.stats\
    --timestamp\
    --rerun-incomplete\
    -j 100\
    --cluster 'sbatch --partition={cluster.partition} --ntasks={cluster.cores} --mem={cluster.mem} --time={cluster.time} -o {cluster.logout} -e {cluster.logerror}'
