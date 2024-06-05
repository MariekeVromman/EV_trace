#!/bin/bash
#Torque Configuration
#PBS -l walltime=1:00:00
#PBS -l mem=16gb
#PBS -l nodes=1:ppn=1
#PBS -q batch
#PBS -N FC_D1515T43
#PBS -j oe
#PBS -o /data/tmp/mvromman/01_EV_track/data_rc/rec_comp_out.txt


#/data/users/mvromman/seqkit seq -r -p -t DNA /data/tmp/mvromman/01_EV_track/data/D1515T43/D1515T43.R1.fastq.gz | gzip -c > /data/tmp/mvromman/01_EV_track/data_rc/D1515T43/D1515T43.R1_rc.fastq.gz
#/data/users/mvromman/seqkit seq -r -p -t DNA /data/tmp/mvromman/01_EV_track/data/D1515T44/D1515T44.R1.fastq.gz | gzip -c > /data/tmp/mvromman/01_EV_track/data_rc/D1515T44/D1515T44.R1_rc.fastq.gz
#/data/users/mvromman/seqkit seq -r -p -t DNA /data/tmp/mvromman/01_EV_track/data/D1515T45/D1515T45.R1.fastq.gz | gzip -c > /data/tmp/mvromman/01_EV_track/data_rc/D1515T45/D1515T45.R1_rc.fastq.gz
#/data/users/mvromman/seqkit seq -r -p -t DNA /data/tmp/mvromman/01_EV_track/data/D1515T46/D1515T46.R1.fastq.gz | gzip -c > /data/tmp/mvromman/01_EV_track/data_rc/D1515T46/D1515T46.R1_rc.fastq.gz
#/data/users/mvromman/seqkit seq -r -p -t DNA /data/tmp/mvromman/01_EV_track/data/D1515T47/D1515T47.R1.fastq.gz | gzip -c > /data/tmp/mvromman/01_EV_track/data_rc/D1515T47/D1515T47.R1_rc.fastq.gz
#/data/users/mvromman/seqkit seq -r -p -t DNA /data/tmp/mvromman/01_EV_track/data/D1515T48/D1515T48.R1.fastq.gz | gzip -c > /data/tmp/mvromman/01_EV_track/data_rc/D1515T48/D1515T48.R1_rc.fastq.gz
