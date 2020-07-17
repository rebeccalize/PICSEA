#!/bin/ksh

#After donwloading ECMWF precipitation data in .grb format, this script converts the files for one month into netcdf format, storing data as floats
#grib_to_netcdf (eccodes) doesn't work because the grid isn't regular. The cdo command was sent by Ervin and seems to work on test files
#FYI, the full resolution precipitation files are 2.1GB each...

year=$1
month=$2
month=`printf "%02d\n" $month`

#find the total number of days in this month of this year
nodays=`cal ${month} ${year} | awk 'NF {DAYS = $NF}; END {print DAYS}'`

cd /scratch/mo/more/idai_kenneth_ecmwf_forecast_tc_related_precip

module load cdo
    
for d in $(seq 1 19); do #   $(seq 1 19);    
    #leading zero
    d=`printf "%02d\n" $d`
          
    for h in 0 12 ; do        
        #leading zero
        h=`printf "%02d\n" $h`
	
	echo ${year}${month}${d}${h}
	
	
	cdo enssum idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_006.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_012.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_018.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_024.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_0days.nc
	
	cdo enssum idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_030.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_036.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_042.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_048.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_1days.nc
	
	cdo enssum idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_054.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_060.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_066.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_072.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_2days.nc
	
	cdo enssum idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_078.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_084.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_090.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_096.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_3days.nc
	
	cdo enssum idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_102.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_108.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_114.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_120.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_4days.nc
	
	cdo enssum idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_126.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_132.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_138.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_144.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_5days.nc
	
	cdo enssum idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_150.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_156.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_162.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_168.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_6days.nc
	
	cdo enssum idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_174.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_180.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_186.nc idai_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_192.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_ensmean_${year}${month}${d}${h}_lt_7days.nc
	
	
	##################################################
	
	
	#cdo enssum idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_006.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_012.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_018.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_024.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_0days.nc
	
	#cdo enssum idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_030.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_036.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_042.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_048.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_1days.nc
	
	#cdo enssum idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_054.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_060.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_066.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_072.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_2days.nc
	
	#cdo enssum idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_078.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_084.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_090.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_096.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_3days.nc
	
	#cdo enssum idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_102.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_108.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_114.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_120.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_4days.nc
	
	#cdo enssum idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_126.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_132.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_138.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_144.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_5days.nc
	
	#cdo enssum idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_150.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_156.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_162.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_168.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_6days.nc
	
	#cdo enssum idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_174.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_180.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_186.nc idai_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_192.nc 24h_precip/idai_24h_tc_related_precip_ecmwf_det_${year}${month}${d}${h}_lt_7days.nc
            
    done
        
done
