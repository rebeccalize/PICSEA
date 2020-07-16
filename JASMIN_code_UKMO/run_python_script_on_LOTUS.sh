#!/bin/bash

#Submits an array of jobs to LOTUS
#To run: bsub < run_python_script_bsub.sh 

#BSUB -oo mjo_stats-%J-%I.o
#BSUB -eo mjo_stats-%J-%I.e

#BSUB -W 05:00 #estimated run time in HH:MM



#python2.7 compute_UKMO_deterministic_track_intensity_statistics_per_MJO_phase_all_years_with_confidence_intervals.py 8 1

#python2.7 compute_UKMO_control_ensmean_track_intensity_statistics_per_MJO_phase_with_confidence_intervals_CORRECT.py 8 1

#python2.7 compute_UKMO_ensemble_members_track_intensity_statistics_per_MJO_phase_with_confidence_intervals_CORRECT.py 8 1

#python2.7 compute_UKMO_deterministic_intensity_values_per_MJO_phase_all_years_with_confidence_intervals.py 8 1


python2.7 extract_GPM-IMERG_TC-related_precip_for_ECMWF_analysis.py
	
