#!/bin/bash
#Torque Configuration
#PBS -l walltime=6:00:00
#PBS -l mem=16gb
#PBS -l nodes=1:ppn=16
#PBS -q batch
#PBS -N SlamDunk_D1515T46.R1
#PBS -j oe
#PBS -o /data/tmp/mvromman/01_EV_track/20231019_rc_run/D1515T46.R1/out_SlamDunk_D1515T46.R1.txt

singularity exec docker://tobneu/slamdunk slamdunk all \
	--reference /data/tmp/mvromman/01_EV_track/references/GRCh38.p14.genome.fa \
	--bed /data/tmp/mvromman/01_EV_track/references/gencode.v44.chr_patch_hapl_scaff.annotation.bed \
	--outputDir /data/tmp/mvromman/01_EV_track/20231019_rc_run/D1515T46.R1 \
	--trim-5p 0 \
	--topn 100 \
	--threads 16 \
	/data/tmp/mvromman/01_EV_track/data_rc/D1515T46/D1515T46.R1_rc.fastq.gz
