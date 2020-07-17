#!/bin/ksh

#Copy the track files across to the directory for reformatting them to use for stats and plotting

y1=$1
y2=$2



#iterate over months, days and 00/12 for every day of this year
for month in 11 12 1 2 3 4 5 6 ; do    #7 8 9 10 11 12 1 2 3 4 5 6
    month=`printf "%02d\n" $month`
    
    #find the total number of days in this month of this year
    nodays=`cal ${month} ${year} | awk 'NF {DAYS = $NF}; END {print DAYS}'`
    
    
    #copying them into directories from july of one year to end of june the following year, 
    #so specify years to use for each month
    if [[ ${month} -ge 7 ]] ; then
    	
	year=${y1}
	
    elif [[ ${month} -le 6 ]] ; then
    
    	year=${y2}
	    
    fi
	
    echo ${y1} ${y2}
    
    for d in $(seq 1 $nodays); do       #$nodays 
        #leading zero
        d=`printf "%02d\n" $d`
        
        for h in 0 12 ; do        
	
		#leading zero
		h=`printf "%02d\n" $h`
    
   		#ALWAYS USE THE MATCH-ECMWF-ANALYSIS-IBT-3WAY-ALLVAR DIRECTORY!! anything else either doesn't have intensity data, or isn't 3-way matched
		
		cd /perm/mo/more/TIGGE/Y${y1}${y2}/${year}${month}${d}${h}/MATCH-ECMWF-ANALYSIS-IBT-3WAY-ALLVAR/
	    
	        #iterate over each track file in the matched directory
	        for filename in ./trmatch_cntl_tr* ; do
	    
	    	    echo ${filename} 	
	
	    	    #get the track number from the filename	    	    
	    	    #trnum=${filename//[^0-9]/}
		    
		    trnum1=${filename:15} #removes first 15 digits from the fileame (i.e. "./tr_trmatch_cntl" is removed)
	    
	    	    echo ${trnum}
	    
	    	    savedir=/perm/mo/more/TIGGE/Y${y1}${y2}/reformatted_track_files_per_storm_correct
	    
	    	    if [ ! -d ${savedir} ]; then 
	    	        mkdir ${savedir};
	    	    fi
	    
                    cp ${filename} ${savedir}/${year}${month}${d}${h}_${trnum1}.txt #${year}${month}${d}${h}_tr${trnum}.txt
		
		
		done
		
		#deterministic:
		
		cd /perm/mo/more/TIGGE/Y${y1}${y2}/${year}${month}${d}${h}/MATCH-ECMWF-ANALYSIS-IBT-3WAY-ALLVAR/DET/
	    
	        #iterate over each track file in the matched directory
	        for filename in ./trmatch_cntl_tr* ; do
	    
	    	    echo ${filename} 	
	
	    	    #get the track number from the filename	    	    
	    	    #trnum=${filename//[^0-9]/}
		    
		    trnum2=${filename:15} #removes first 15 digits from the fileame (i.e. "./tr_trmatch_cntl" is removed)
	    
	    	    echo ${trnum2}
	    
	    	    savedir=/perm/mo/more/TIGGE/Y${y1}${y2}/reformatted_track_files_per_storm_correct
	    
	    	    if [ ! -d ${savedir} ]; then 
	    	        mkdir ${savedir};
	    	    fi
	    
                    cp ${filename} ${savedir}/${year}${month}${d}${h}_${trnum2}_det.txt #${year}${month}${d}${h}_tr${trnum}.txt
		
		
            
                done
            
        done
        
    done
    
done
    
