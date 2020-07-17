#!/bin/ksh

#sum the precip forecasts at each lead time, 0-7 days ahead, across the whole storm

#year=2019
#month=3
#month=`printf "%02d\n" $month`

#find the total number of days in this month of this year
#nodays=`cal ${month} ${year} | awk 'NF {DAYS = $NF}; END {print DAYS}'`

cd /scratch/mo/more/idai_kenneth_ecmwf_forecast_tc_related_precip/24h_precip

module load cdo
    
for lt in $(seq 0 2); do #      
	#leading zero
	d=`printf "%02d\n" $d`
        
	
	##################################################
	
	cdo enssum idai_24h_tc_related_precip_ecmwf_det_201903[0-1][0-9]00_lt_${lt}days.nc idai_precip_sum_00runs_det_lt_${lt}days.nc

	cdo enssum idai_24h_tc_related_precip_ecmwf_det_201903[0-1][0-9]12_lt_${lt}days.nc idai_precip_sum_12runs_det_lt_${lt}days.nc
	

	
	cdo enssum idai_precip_sum_00runs_det_lt_${lt}days.nc idai_precip_sum_12runs_det_lt_${lt}days.nc idai_precip_sum_det_DOUBLE_lt_${lt}days.nc
	
	cdo divc,2 idai_precip_sum_det_DOUBLE_lt_${lt}days.nc idai_precip_sum_det_lt_${lt}days.nc
	
	cdo divc,19 idai_precip_sum_det_lt_${lt}days.nc idai_precip_mmday_det_lt_${lt}days.nc
	
	
	##################################################
	
	
    
        
done




for lt in 3 ; do #      
	#leading zero
	d=`printf "%02d\n" $d`
        
	
	##################################################
	
	cdo enssum idai_24h_tc_related_precip_ecmwf_det_201903[0-1][0-8]00_lt_${lt}days.nc idai_precip_sum_00runs_det_lt_${lt}days.nc

	cdo enssum idai_24h_tc_related_precip_ecmwf_det_201903[0-1][0-8]12_lt_${lt}days.nc idai_precip_sum_12runs_det_lt_${lt}days.nc
	
	
	
	cdo enssum idai_precip_sum_00runs_det_lt_${lt}days.nc idai_precip_sum_12runs_det_lt_${lt}days.nc idai_precip_sum_det_DOUBLE_lt_${lt}days.nc
	
	cdo divc,2 idai_precip_sum_det_DOUBLE_lt_${lt}days.nc idai_precip_sum_det_lt_${lt}days.nc
	
	cdo divc,18 idai_precip_sum_det_lt_${lt}days.nc idai_precip_mmday_det_lt_${lt}days.nc
	
	
	##################################################
	
	

        
done

for lt in 4 5 ; do #      
	#leading zero
	d=`printf "%02d\n" $d`
        
	
	##################################################
	
	cdo enssum idai_24h_tc_related_precip_ecmwf_det_201903[0-1][0-7]00_lt_${lt}days.nc idai_precip_sum_00runs_det_lt_${lt}days.nc

	cdo enssum idai_24h_tc_related_precip_ecmwf_det_201903[0-1][0-7]12_lt_${lt}days.nc idai_precip_sum_12runs_det_lt_${lt}days.nc
	
	
	
	cdo enssum idai_precip_sum_00runs_det_lt_${lt}days.nc idai_precip_sum_12runs_det_lt_${lt}days.nc idai_precip_sum_det_DOUBLE_lt_${lt}days.nc
	
	cdo divc,2 idai_precip_sum_det_DOUBLE_lt_${lt}days.nc idai_precip_sum_det_lt_${lt}days.nc
	
	cdo divc,17 idai_precip_sum_det_lt_${lt}days.nc idai_precip_mmday_det_lt_${lt}days.nc
	
	
	##################################################


        
done
