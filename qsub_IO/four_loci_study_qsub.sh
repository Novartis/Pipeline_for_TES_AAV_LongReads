#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N four_loci_study
#$ -l h_rt=100000
#$ -l m_mem_free=100G
#$ -pe smp 10
#$ -o four_loci_study.out
#$ -e four_loci_study.err


snakemake all --rerun-incomplete --keep-going --keep-incomplete --printshellcmds --jobs 40 --snakefile ../workflow/Snakefile.smk --configfile ../config/four_loci_study_config.yaml



