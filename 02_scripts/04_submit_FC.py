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
		"#PBS -l walltime=1:00:00",
		"#PBS -l mem=16gb",
		"#PBS -l nodes=1:ppn=6",
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

def FC(sampleID, out_dir):

	command = "/bioinfo/local/build/Centos/subread/subread-1.5.1-Linux-x86_64/bin/featureCounts \\; "
	command = command + "\t-p -O  \\; "
	command = command + "\t-T 6 \\; "
	command = command + "\t-a /data/tmp/mvromman/01_EV_track/references/gencode.v44.chr_patch_hapl_scaff.annotation.gtf \\; "
	command = command + "\t-o /data/tmp/mvromman/01_EV_track/20231019_FC/{0}/{0}_FC.txt \\; ".format(sampleID)
	command = command + "\t/data/tmp/mvromman/01_EV_track/20231018_STAR/{0}/{0}_Aligned.out.sam; ".format(sampleID)
	command = command + ""

	job_id = sbatch('FC_' + sampleID, command, out_dir)
	
	return job_id

for samplename in os.listdir(base_dir):
	if os.path.isdir(os.path.join(base_dir,samplename)):

		# for R1
		output_dir = output_dir_base + "/" + samplename
		os.makedirs(output_dir, exist_ok=True)
		job = FC(samplename, output_dir)
