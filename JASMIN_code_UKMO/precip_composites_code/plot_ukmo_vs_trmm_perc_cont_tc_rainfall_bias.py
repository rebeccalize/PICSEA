import picsea_library as pl
from netCDF4 import Dataset
import numpy as np


### Script to plot the track density of a particular TC season ###

ukmodir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_precip_composites/"
trmmdir = "/gws/nopw/j04/klingaman/emerton/ibtracs_precip_composites/"

res_list = ['n320', 'n512', 'n768']

# arrays to hold data sum for all years
#pcp_sum, lons_overall, lats_overall = (np.zeros((400, 1440)) for i in range(3))

#cdays_sum = 0

#map the precipitation composite for one year/TC season, in mm/day
for res in res_list:

	for lt in [0,1,2,3,4,5,6]: #

		if res == 'n320':
			ukmofile = ukmodir + str(lt) + "/ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.072006-032010_n320.txt"
			trmmfile = trmmdir + "trmm_ibtracs.perc_cont_tc_rainfall.072006-032010_n320.txt"
			outfile_bias = "ukmo_bias.perc_cont_tc_rainfall.072006-032010.n320."+str(lt)+"-days.png"
			lonlatfile = trmmdir + "trmm.comp_pcp_tc.072006-032010_regridded_n320.nc"
			fflatlon = Dataset(lonlatfile,'r')
			lons = fflatlon.variables['longitude'][:]
			lats = fflatlon.variables['latitude'][:]			
			
		elif res == 'n512':
			ukmofile = ukmodir + str(lt) + "/ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.032010-072014_n512.txt"
			trmmfile = trmmdir + "trmm_ibtracs.perc_cont_tc_rainfall.032010-072014_n512.txt"
			outfile_bias = "ukmo_bias.perc_cont_tc_rainfall.032010-072014.n512."+str(lt)+"-days.png"
			lonlatfile = trmmdir + "trmm.comp_pcp_tc.032010-072014_regridded_n512.nc"
			fflatlon = Dataset(lonlatfile,'r')
			lons = fflatlon.variables['longitude'][:]
			lats = fflatlon.variables['latitude'][:]
			
		elif res == 'n768':
			ukmofile = ukmodir + str(lt) + "/ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.072014-122016_n768.txt"
			trmmfile = trmmdir + "trmm_ibtracs.perc_cont_tc_rainfall.072014-122016_n768.txt"
			outfile_bias = "ukmo_bias.perc_cont_tc_rainfall.072014-122016.n768."+str(lt)+"-days.png"
			lonlatfile = trmmdir + "trmm.comp_pcp_tc.072014-122016_regridded_n768.nc"
			fflatlon = Dataset(lonlatfile,'r')
			lons = fflatlon.variables['longitude'][:]
			lats = fflatlon.variables['latitude'][:]

		ukmo_data=np.genfromtxt(ukmofile)
		trmm_data=np.genfromtxt(trmmfile)

		bias = ukmo_data - trmm_data

		print ukmo_data
		print trmm_data
		print bias

		label = " "#""UKMO - TRMM precipitation bias (mm/day) " + res + " " + str(lt) + " days lead time"
		bounds = np.linspace(-50, 50, 11)
		#pl.map_composite_data(bias, lats, lons, outfile_bias, "SIO2", 'PuOr', label, bounds,'bias')


#here, add the composites for the total across the entire period, once remapping is finished

#here, add the composites for the total across the entire period, once remapping is finished
for lt in [0,1,2,3,4,5,6]:
	ukmofile = ukmodir + str(lt) + "/ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.072006-122016_n320.txt"
	trmmfile = trmmdir + "trmm_ibtracs.perc_cont_tc_rainfall.072006-122016_n320.txt"
	outfile_bias = "ukmo_bias.perc_cont_tc_rainfall.072006-122016.n320." + str(lt) + "-days.png"
	lonlatfile = trmmdir + "trmm.comp_pcp_tc.072006-122016_n320.nc"
	fflatlon = Dataset(lonlatfile,'r')
	lons = fflatlon.variables['longitude'][:]
	lats = fflatlon.variables['latitude'][:]

	ukmo_data=np.genfromtxt(ukmofile)
	trmm_data=np.genfromtxt(trmmfile)

	bias = ukmo_data - trmm_data

	print ukmo_data
	print trmm_data
	print bias

	label = " "#""UKMO - TRMM precipitation bias (mm/day) " + res + " " + str(lt) + " days lead time"
	bounds = np.linspace(-25, 25, 11)
	#pl.map_composite_data(bias, lats, lons, outfile_bias, "SIO2", 'PuOr', label, bounds,'bias')


	for season in ['DJFM','NDJFMA','MJJASO']:
		ukmofile = ukmodir + str(lt) + "/ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.072006-122016"+season+".n320.txt"
		trmmfile = trmmdir + "trmm_ibtracs.perc_cont_tc_rainfall.072006-122016."+season+"_n320.txt"
		outfile_bias = "ukmo_bias.perc_cont_tc_rainfall.072006-122016."+str(lt)+"-days."+season+".n320.png"
		lonlatfile = trmmdir+"trmm.comp_pcp_tc.072006-122016_n320.nc"
		fflatlon = Dataset(lonlatfile, 'r')
		lons = fflatlon.variables['longitude'][:]
		lats = fflatlon.variables['latitude'][:]

		ukmo_data = np.genfromtxt(ukmofile)
		trmm_data = np.genfromtxt(trmmfile)

		bias = ukmo_data - trmm_data

		label = " "  # ""UKMO - TRMM precipitation bias (mm/day) " + res + " " + str(lt) + " days lead time"
		bounds = np.linspace(-25, 25, 11)
		pl.map_composite_data(bias, lats, lons, outfile_bias, "SIO2", 'PuOr', label, bounds, 'bias')
