#!/bin/bash
#$ -l mem=15G,time=200::
#$ -cwd -S /bin/bash -N myJob
Rscript 10_projection_cluster_uncertainty_GAM.R $1 $2
