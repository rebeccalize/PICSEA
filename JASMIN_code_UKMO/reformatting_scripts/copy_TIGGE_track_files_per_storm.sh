#!/bin/ksh

y1=$1
y2=$2

for month in 7 8 9 10 11 12 1 2 3 4 5 6 ; do #

	month=`printf "%02d\n" $month`
	
	nodays=`cal ${month} ${year} | awk 'NF {DAYS = $NF}; END {print DAYS}'`
	
	
	if [[ ${month} -ge 7 ]] ; then
		
		year=${y1}
		
	elif [[ ${month} -le 6 ]] ; then
	
		year=${y2}
		
	fi
	
	echo ${y1} ${y2}
	echo ${month}
	
	
	for d in $(seq 1 $nodays); do
	
		d=`printf "%02d\n" $d`
		
		for h in 0 12 ; do
		
			h=`printf "%02d\n" $h`
			
			cd /gws/nopw/j04/klingaman/emerton/TIGGE/UKMO/Y${y1}${y2}/${year}${month}${d}${h}/MATCH-ECMWF-ANALYSIS-IBT-3WAY/
			
			for filename in ./trmatch_cntl_tr* ; do
			
				echo ${filename}
			
				trnum=${filename:15}
			
				echo ${trnum}
			
				savedir=/gws/nopw/j04/klingaman/emerton/ukmo_TIGGE_analysis/reformatted_track_files_per_storm/Y${y1}${y2}
			
				if [ ! -d ${savedir} ] ; then
					mkdir ${savedir}
				fi
			
				cp ${filename} ${savedir}/${year}${month}${d}${h}_${trnum}.txt
				
			done
			
		done
		
	done
	
done
