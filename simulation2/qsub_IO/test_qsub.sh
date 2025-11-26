#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N qsub_test
#$ -l h_rt=100000
#$ -l m_mem_free=100G
#$ -pe smp 10
#$ -o test.out
#$ -e test.err

cd ..
snakemake all --configfile config/test_config.yaml


