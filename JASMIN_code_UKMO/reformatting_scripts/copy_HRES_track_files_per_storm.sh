#!/bin/ksh

#Copy the track files across to the directory for reformatting them to use for stats and plotting

year=$1

#FOR 2018, THIS HAS BEEN RUN FOR JAN - JUNE
#BUT NOT YET JULY - DEC, AS WE DIDN'T YET HAVE IBTRACS TO RUN MATCHING!


#iterate over months, days and 00/12 for every day of this year
for month in $(seq 1 5); do    #leading zero
    month=`printf "%02d\n" $month`
    
    #find the total number of days in this month of this year
    nodays=`cal ${month} ${year} | awk 'NF {DAYS = $NF}; END {print DAYS}'`
    
    
    #copying them into directories from july of one year to end of june the following year, 
    #so specify years to use for each month
    if [[ ${month} -ge 7 ]] ; then
    	
	    y1=${year}
	    ((y2=${year}+1))
	
    elif [[ ${month} -le 6 ]] ; then
    
    	((y1=${year}-1))
	    y2=$year
	    
	fi
	
	echo ${y1} ${y2}
    
    for d in $(seq 1 $nodays); do        
        #leading zero
        d=`printf "%02d\n" $d`
        
        for h in 0 12 ; do        
            #leading zero
            h=`printf "%02d\n" $h`
    
   	
	        cd /gws/nopw/j04/klingaman/emerton/TRACK/UKMO/Y${year}-SH/UKMO_${year}${month}${d}${h}_VOR_VERTAVG_T63_DET/MATCH-ECMWF-ANALYSIS-IBT-3WAY/
	    
	        #iterate over each track file in the matched directory
	        for filename in ./trmatch_cntl_tr* ; do
	    
	    	    echo ${filename} 	
	
	    	    #get the track number from the filename	    	    
	    	    trnum=${filename//[^0-9]/}
	    
	    	    #echo ${trnum}
	    
	    	    savedir=/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/${y1}_${y2}/tr${trnum}
	    
	    	    if [ ! -d ${savedir} ]; then 
	    	        mkdir ${savedir};
	    	    fi
	    
                cp ${filename} ${savedir}/${year}${month}${d}${h}_tr${trnum}.txt
            
            done
            
        done
        
    done
    
done
    
