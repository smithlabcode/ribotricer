#PBS -S /bin/bash
#PBS -q cmb
#PBS -e /home/cmb-06/as/wenzhenl/logs
#PBS -o /home/cmb-06/as/wenzhenl/logs
#PBS -l vmem=100G
#PBS -l pmem=100G
#PBS -l mem=100G
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -N jupyter
cd /staging/as/wenzhenl/jupyter
jupyter notebook --port=11111 --no-browser --NotebookApp.iopub_data_rate_limit=10000000000
