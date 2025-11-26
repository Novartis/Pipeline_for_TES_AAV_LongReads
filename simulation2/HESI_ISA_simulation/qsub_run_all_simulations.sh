#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N hesi_sim
#$ -l h_rt=100000
#$ -l m_mem_free=250G
#$ -pe smp 40
#$ -o qsub.out
#$ -e qsub.err

./run_all_simulations.sh

