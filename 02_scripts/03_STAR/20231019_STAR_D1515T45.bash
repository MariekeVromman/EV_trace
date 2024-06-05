#!/bin/bash
#Torque Configuration
#PBS -l walltime=6:00:00
#PBS -l mem=64gb
#PBS -l nodes=1:ppn=6
#PBS -q batch
#PBS -N STAR_D1515T45
#PBS -j oe
#PBS -o /data/tmp/mvromman/01_EV_track/20231018_STAR/D1515T45/out_STAR_D1515T45.txt

/bioinfo/local/build/Centos/STAR/STAR-2.7.5a/bin/Linux_x86_64/STAR \
	--runThreadN 6 \
	--readFilesIn /data/tmp/mvromman/01_EV_track/data/D1515T45/D1515T45.R1.fastq.gz /data/tmp/mvromman/01_EV_track/data/D1515T45/D1515T45.R2.fastq.gz \
	--genomeDir /data/tmp/mvromman/01_EV_track/references/STAR_index \
	--outFileNamePrefix /data/tmp/mvromman/01_EV_track/20231018_STAR/D1515T45/D1515T45_ \
	--readFilesCommand zcat

