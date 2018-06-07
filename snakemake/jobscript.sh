#!/usr/bin/env bash
DIR=$(dirname "$(readlink -f "$0")")

echo "=============================="
cd $SLURM_SUBMIT_DIR
echo "=============================="
echo $SLURM_SUBMIT_DIR
echo "=============================="
echo "=============================="
echo $PATH
echo "=============================="
echo $SLURM_JOB_ID
echo "=============================="

{exec_job}

# Report resource consumption because it's not reported by default
echo "------------------------------"
scontrol show job $SLURM_JOB_ID

# if the job succeeds, snakemake 
# touches jobfinished, thus if it exists cat succeeds. if cat fails, the error code indicates job failure
# an error code of 100 is needed since UGER only prevents execution of dependent jobs if the preceding
# job exits with error code 100

cat $1 &>/dev/null && exit 0 || exit 100
