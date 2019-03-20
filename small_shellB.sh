#!/bin/bash
#$ -l mem=15G,time=55::
#$ -cwd -S /bin/bash -N myJob
Rscript 7b_projection_clusterB_PA.R $1

