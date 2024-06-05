#!/usr/bin/python3

#ml load Python/3.6.6-intel-2018b

import os
import argparse
import subprocess


parser = argparse.ArgumentParser(description='test')

parser.add_argument('-b', nargs=1, required=True, help='The base directory where the sample subdirectories are located', metavar='base_dir')
parser.add_argument('-o', nargs=1, required=True, help='The directory where output should be created', metavar='output_dir')


args = parser.parse_args()
base_dir = args.b[0]
output_dir_base = args.o[0]


def sbatch(job_name, command, outp_dir):
	
	printlines = [
		"#!/bin/bash",
		"#Torque Configuration",
		"#PBS -l walltime=6:00:00",
		"#PBS -l mem=16gb",
		"#PBS -l nodes=1:ppn=16",
		"#PBS -q batch",
		"#PBS -N {0}".format(job_name),
		"#PBS -j oe",
		"#PBS -o {0}/out_{1}.txt\n".format(outp_dir, job_name)
		]


	printlines.extend(command.split("; "))

	printjobfile('{0}/20231019_{1}.bash'.format(outp_dir, job_name),printlines)

	sbatch_command = '{0}/20231019_{1}.bash'.format(outp_dir, job_name)
	sbatch_response = subprocess.getoutput(sbatch_command)
	print(sbatch_response)
	job_id = sbatch_response.split('.')[0].strip()
	return job_id

def printjobfile(filename, printlines):
	with open(filename, 'w') as the_file:
		for line in printlines:
			the_file.write(line + "\n")

def slamdunk(sampleID, out_dir, read):

	command = "singularity exec docker://tobneu/slamdunk slamdunk all \\; "
	command = command + "\t--reference /data/tmp/mvromman/01_EV_track/references/GRCh38.p14.genome.fa \\; "
	command = command + "\t--bed /data/tmp/mvromman/01_EV_track/references/gencode.v44.chr_patch_hapl_scaff.annotation.bed \\; "
	command = command + "\t--outputDir /data/tmp/mvromman/01_EV_track/20231019_rc_run/{0}.{1} \\; ".format(sampleID, read)
	command = command + "\t--trim-5p 0 \\; "
	command = command + "\t--topn 100 \\; "
	command = command + "\t--threads 16 \\; "
	#command = command + "\t--multimap \\; "

	if read == 'R1':
		command = command + "\t{0}_rc/{1}/{1}.{2}_rc.fastq.gz".format(base_dir, sampleID, read)

	elif read == 'R2':
		command = command + "\t{0}/{1}/{1}.{2}.fastq.gz".format(base_dir, sampleID, read)
	
	command = command + ""

	job_id = sbatch('SlamDunk_' + sampleID + '.' + read, command, out_dir)
	
	return job_id

for samplename in os.listdir(base_dir):
	if os.path.isdir(os.path.join(base_dir,samplename)):

		# for R1
		read_nr = 'R1'
		output_dir = output_dir_base + "/" + samplename + '.R1'
		os.makedirs(output_dir, exist_ok=True)
		job = slamdunk(samplename, output_dir, read_nr)

		# for R2
		read_nr = 'R2'
		output_dir = output_dir_base + "/" + samplename + '.R2'
		os.makedirs(output_dir, exist_ok=True)
		job = slamdunk(samplename, output_dir, read_nr)