#!/bin/ksh


module load cdo
    
for year in 2018 ; do #   $(seq 1 19);
	  
    cd /vol/floods/more/rebecca_TC_do_not_delete/PICSEA_DATA/PRECIP/monthly_files_nc/${year}/ENS
          
    for month in 1 2 3 4 5 6 7 8 9 10 11 12 ; do        
        #leading zero
        month=`printf "%02d\n" $month`
	
	echo ${year}${month}
	
	gunzip ECMWF_PRECIP_${year}${month}_ENS.nc.gz
	
	cdo splitday ECMWF_PRECIP_${year}${month}_ENS.nc ECMWF_PRECIP_${year}${month}_ENS.d
	
	gzip *.nc
	   
    done
        
done
