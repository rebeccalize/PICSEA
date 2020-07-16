import picsea_library as pl
from netCDF4 import Dataset
import numpy as np



ecdir = "/gws/nopw/j04/klingaman/emerton/ecmwf/idai_precip_analysis_ecmwf/"
gpmdir = "/gws/nopw/j04/klingaman/emerton/2018-2019_GPM-IMERG_precip_composites/2019/"
# arrays to hold data sum for all years
#pcp_sum, lons_overall, lats_overall = (np.zeros((400, 1440)) for i in range(3))

#cdays_sum = 0

#map the precipitation composite for one year/TC season, in mm/day

for lt in [0,2,3,4,5]: #0,1,2,3,5

	if lt == 5:
		cdays = 18
	
	else:
		cdays = 19


	#ecfile = ecdir + "/idai_precip_mmday_ensmean_lt_"+str(lt)+"days.nc"
	ecfile = ecdir + "/idai_precip_sum_mm_det_lt_"+str(lt)+"days_regridded_to_imerg.nc"
	gpmfile = gpmdir + "idai_analysis.comp_pcp_tc.201903_lt_"+str(lt)+"days.nc"
	outfile_ec = "idai_ecmwf_DET_mmday.comp_pcp_tc."+str(lt)+"-days.png"
	outfile_gpm = "idai_analysis_gpm_imerg_mmday.comp_pcp_tc."+str(lt)+"-days.png"
	outfile_bias = "idai_ecmwf_DET_vs_gpm_imerg_precip_bias_mmday.comp_pcp_tc."+str(lt)+"-days.png"


	ffec = Dataset(ecfile, 'r')
	ecpcp = ffec.variables['tc_precip'][:]
	lons = ffec.variables['lon'][:]
	lats = ffec.variables['lat'][:]

	#print np.shape(ecpcp)

	ffgpm = Dataset(gpmfile, 'r')
	gpmpcp = ffgpm.variables['pcp'][:]
	gpmpcp = np.transpose(gpmpcp,[1,0])
	
	ecpcp_mmday = ecpcp / cdays
	gpmpcp_mmday = gpmpcp / cdays
	
	bias = ecpcp_mmday - gpmpcp_mmday
	
	print "bias max: ", np.nanmax(bias)
	print "bias min: ", np.nanmin(bias)

	#print ecpcp
	#print gpmpcp_mmday
	#print bias

	label1 = " "#""UKMO precipitation (mm/day) "+res+" "+str(lt)+" days lead time"
	label2 = " "#""TRMM precipitation (mm/day) "+res
	label3 = " "#""UKMO - TRMM precipitation bias (mm/day) " + res + " " + str(lt) + " days lead time"
	bounds1 = np.linspace(0, 60, 7)
	bounds3 = np.linspace(-30,30,7)
	pl.map_composite_data(ecpcp_mmday, lats, lons, outfile_ec, "SIO2", 'Purples', label1, bounds1)
	pl.map_composite_data(gpmpcp_mmday, lats, lons, outfile_gpm, "SIO2", 'Purples', label2, bounds1)
	pl.map_composite_data(bias, lats, lons, outfile_bias, "SIO2", 'RdBu', label3, bounds3,'bias')
	
	
	#AND THEN FOR SUM, RATHER THAN MM/DAY:
	
	ecfilesum = ecdir + "/idai_precip_sum_mm_det_lt_"+str(lt)+"days_regridded_to_imerg.nc"
	
	ffecsum = Dataset(ecfilesum, 'r')
	ecpcpsum = ffecsum.variables['tc_precip'][:]
	
	print np.nanmax(ecpcpsum)
	
	biassum = ecpcpsum - gpmpcp
	
	print "biassum max: ", np.nanmax(biassum)
	print "biassum min: ", np.nanmin(biassum)
	
	outfile_ecsum = "idai_ecmwf_DET_mm_total.comp_pcp_tc."+str(lt)+"-days_800mmlegend.png"
	outfile_gpmsum = "idai_analysis_gpm_imerg_mm_total.comp_pcp_tc."+str(lt)+"-days_800mmlegend.png"
	outfile_biassum = "idai_ecmwf_DET_vs_gpm_imerg_precip_bias_mm_total.comp_pcp_tc."+str(lt)+"-days.png"
	
	bounds1 = np.linspace(0, 800, 9)
	#bounds2 = np.linspace(0,2000,11)
	bounds3 = np.linspace(-500,500,11)
	
	pl.map_composite_data(ecpcpsum, lats, lons, outfile_ecsum, "SIO2", 'Purples', label1, bounds1)
	pl.map_composite_data(gpmpcp, lats, lons, outfile_gpmsum, "SIO2", 'Purples', label2, bounds1)
	pl.map_composite_data(biassum, lats, lons, outfile_biassum, "SIO2", 'RdBu', label3, bounds3,'bias')



