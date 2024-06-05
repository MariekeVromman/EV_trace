
#!/bin/bash
#Torque Configuration
#PBS -l walltime=6:00:00
#PBS -l mem=16gb
#PBS -l nodes=1:ppn=6
#PBS -q batch
#PBS -N bam_D1515T44
#PBS -j oe
#PBS -o /data/tmp/mvromman/01_EV_track/20231019_rc_run/sort_bam_out_D1515T44.txt



# sort original BAM file
#/bioinfo/local/build/Centos/samtools/samtools-1.9/bin/samtools sort -@ 6 /data/tmp/mvromman/01_EV_track/20231019_rc_run/D1515T44.R1/map/D1515T44.R1*.bam -o /data/tmp/mvromman/01_EV_track/20231019_rc_run/D1515T44.R1/map/D1515T44.R1_sorted.bam

# generate index for original BAM file
/bioinfo/local/build/Centos/samtools/samtools-1.9/bin/samtools index /data/tmp/mvromman/01_EV_track/20231019_rc_run/D1515T44.R1/map/D1515T44.R1_sorted.bam

# get region of interest SAM file
#mkdir /data/tmp/mvromman/01_EV_track/20231019_rc_run/D1515T44.R1/spec_transcript
/bioinfo/local/build/Centos/samtools/samtools-1.9/bin/samtools view -h /data/tmp/mvromman/01_EV_track/20231019_rc_run/D1515T44.R1/map/D1515T44.R1_sorted.bam "chr1:825138-849592" > /data/tmp/mvromman/01_EV_track/20231019_rc_run/D1515T44.R1/spec_transcript/D1515T44.R1.sam

# generate region of interest BAM file
/bioinfo/local/build/Centos/samtools/samtools-1.9/bin/samtools view -bS /data/tmp/mvromman/01_EV_track/20231019_rc_run/D1515T44.R1/spec_transcript/D1515T44.R1.sam > /data/tmp/mvromman/01_EV_track/20231019_rc_run/D1515T44.R1/spec_transcript/D1515T44.R1.bam

# generate index for region of interst BAM file
/bioinfo/local/build/Centos/samtools/samtools-1.9/bin/samtools index /data/tmp/mvromman/01_EV_track/20231019_rc_run/D1515T44.R1/spec_transcript/D1515T44.R1.bam
