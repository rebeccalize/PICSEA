#!/bin/ksh

#Reformat the track files that have been copied across for use in stats and plotting

y1=$1
y2=$2

cd /perm/mo/more/TIGGE/Y${y1}${y2}/reformatted_track_files_per_storm_correct/SIO_storms/

#for 2017-2018, the 2018 tracks do actually have MSLP and UV10 already, but at the time of running this, the 2017 tracks don't, yet
#the copying across (script 1) and this script will need re-running when the MSLP and UV10 have been added to the tracks
#for 2014, 2015, 2016 and 2017, but not 2018 - the 2018 files will just need moving to a more appropriately named directory

for dirname in ./tr0* ; do
#for trnum in 0003 0004 0005 0006 0007 0008 0009 0010 0011 0012 0013 0014 0015 0016 0017 0018 0019 0020 0021 ; do #

    echo ${dirname}
    trnum=${dirname//[^0-9]/}
    echo ${trnum}
    
    
    python2.7 /home/mo/more/PICSEA/reformat_ecmwf_track_files.py ${y1} ${y2} ${trnum}
	
    
done
