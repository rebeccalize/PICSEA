#!/bin/ksh

#After running TRACK for a year of data (UKMO NWP), this script checks whether any of the jobs failed,
#by checking that the final output file (...tcident.new) exists for every date that was run

#year=$1

#iterate over months, days and 00/12 for every day of this year
for year in 2019 ; do
	for month in $(seq 1 12); do    #leading zero
   		python2.7 reformat_analysis_track_files_yearmonth.py ${year} ${month}
	done
done

for year in 2020 ; do
	for month in $(seq 1 6); do    #leading zero
   		python2.7 reformat_analysis_track_files_yearmonth.py ${year} ${month}
	done
done
