#!/bin/bash

#Submits an array of jobs to LOTUS
#To run: bsub < run_python_script_bsub.sh 

#BSUB -oo tot_pcp-%J-%I.o
#BSUB -eo tot_pcp-%J-%I.e

#BSUB -W 12:00 #estimated run time in HH:MM



#for lt in 1 2 3 4 5 6 7 ; do
python2.7 compute_ukmo_nwp_fcst_TOTAL_precip_composites.py 7 #total_pcp

