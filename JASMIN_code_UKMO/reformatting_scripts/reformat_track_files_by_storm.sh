#!/bin/ksh

#Reformat the track files that have been copied across for use in stats and plotting

y1=$1
y2=$2

cd /gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/${y1}_${y2}/

for dirname in ./tr0* ; do

    echo ${dirname}
    trnum=${dirname//[^0-9]/}
    echo ${trnum}
    
    python2.7 /home/users/emerton/analysis_scripts/reformatting_scripts/reformat_ukmo_deterministic_ibtracs_analysis_tracks_for_one_storm.py ${y1} ${y2} ${trnum}
    
done

