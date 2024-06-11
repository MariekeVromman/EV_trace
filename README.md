# EV_trace
This repo consist of all the scripts performed for the EV project in collaboration with Panagiotis Papoutsoglou.

For this project, the [SLAMseq protocol](https://www.lexogen.com/wp-content/uploads/2017/11/059UG142V0102_SLAMseq_User-Guide.pdf) was used to label RNA. This protocol induces a T>C conversion in the RNA strands, which can be used as a label. Regular RNA sequencing is performed and the tagged RNA strands can be identified as 'mutations' using [SlamDunk](https://t-neumann.github.io/slamdunk/), a NextFlow pipeline. See also the presentation in Teams `Marieke/20240404_student_meeting.pptx` for more information about the technique itself.

I had a couple of question on the pipeline, and the developers replied very fast (see mail_SlamDunk.txt), so do not hesitate to contact them.

The samples analysed come from one sequencing run of labeled total cell line RNA (KDI project: EVRNATRACE, sequencing run: D1515). There are 6 samples in total: 2 replicates x 3 different concentrations:
- D1515T43	250_uM_S4U_1
- D1515T44	250_uM_S4U_2
- D1515T45	125_uM_S4U_1
- D1515T46	125_uM_S4U_2
- D1515T47	0_uM_S4U_2
- D1515T48	0_uM_S4U_3

> [!WARNING]  
> The currenct SlamDunk pipeline only detects T>C conversions in SE data. For PE data, the pipeline was run twice: once for R1 and once for R2. For R1, first, the reverse compelent should be taken of the reads (indicated as R1_rc in this repository). This can be done with [SeqKit](https://github.com/shenwei356/seqkit) (example script: ...).

This repository contains 4 folders
1. **data**  
This folders contains the output data from mapping the fastq files with STAR and generating counts with FeatureCounts, and the output data from running the SlamDunk pipeline. As this is a big folder, it is not included in the github repo itself, but it is present on the hard disk. 

2. **scripts**  
This folder contains al scripts used on the cluster to generate the data in the `data` folder. For most tools, a Python script is used to generate a Bash script for each sample. Next, the Bash script are submitted to the cluster.
  
    - first, the `01_seqkit_revcomp.bash` script was used to take the reverse complement of the R1 fastqs (see warning above).
    - then, `02_submit_SlamDunk.py` can be used to generate a bash script for each sample (stored in 02_SlamDunk)
    - also, `STAR` and `FeatureCounts` were run to generate counts for each sample

3. **data-analysis**  
This folder contains the R scripts used to analyse the data further, and to generate the figures.

4. **figures**  
This folder contains all the generated figures. See also the presentation on Teams in `EV lncRNA team/Marieke/EV_labeling_project_updates.pptx`.
