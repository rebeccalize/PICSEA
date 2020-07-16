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
			ukmofile = ukmodir + str(lt) + "/ukmo_nwp.comp_pcp_tc."+str(lt)+"_days.072006-032010.n320.nc"
			trmmfile = trmmdir + "trmm.comp_pcp_tc.072006-032010_regridded_n320.nc"
			outfile_ukmo = "ukmo_nwp.comp_pcp_tc.072006-032010.n320."+str(lt)+"-days.png"
			outfile_trmm = "trmm.comp_pcp_tc.072006-032010.n320.png"
			outfile_bias = "ukmo_bias.comp_pcp_tc.072006-032010.n320."+str(lt)+"-days.png"
		elif res == 'n512':
			ukmofile = ukmodir + str(lt) + "/ukmo_nwp.comp_pcp_tc."+str(lt)+"_days.032010-072014.n512.nc"
			trmmfile = trmmdir + "trmm.comp_pcp_tc.032010-072014_regridded_n512.nc"
			outfile_ukmo = "ukmo_nwp.comp_pcp_tc.032010-072014.n512."+str(lt)+"-days.png"
			outfile_trmm = "trmm.comp_pcp_tc.032010-072014.n512.png"
			outfile_bias = "ukmo_bias.comp_pcp_tc.032010-072014.n512."+str(lt)+"-days.png"
		elif res == 'n768':
			ukmofile = ukmodir + str(lt) + "/ukmo_nwp.comp_pcp_tc."+str(lt)+"_days.072014-122016.n768.nc"
			trmmfile = trmmdir + "trmm.comp_pcp_tc.072014-122016_regridded_n768.nc"
			outfile_ukmo = "ukmo_nwp.comp_pcp_tc.072014-122016.n768."+str(lt)+"-days.png"
			outfile_trmm = "trmm.comp_pcp_tc.072014-122016.n768.png"
			outfile_bias = "ukmo_bias.comp_pcp_tc.072014-122016.n768."+str(lt)+"-days.png"


		ffukmo = Dataset(ukmofile, 'r')
		ukmopcp = ffukmo.variables['pcp'][:]
		lons = ffukmo.variables['longitude'][:]
		lats = ffukmo.variables['latitude'][:]
		cdaysukmo = ffukmo.getncattr('contributing_days') #contributing days is a global attribute rather than a variable

		ukmopcp_mm = ukmopcp / cdaysukmo

		fftrmm = Dataset(trmmfile, 'r')
		trmmpcp = fftrmm.variables['pcp'][:]
		cdaystrmm = fftrmm.getncattr('contributing_days')

		trmmpcp_mm = trmmpcp / cdaystrmm

		bias = ukmopcp_mm - trmmpcp_mm

		print ukmopcp_mm
		print trmmpcp_mm
		print bias

		label1 = " "#""UKMO precipitation (mm/day) "+res+" "+str(lt)+" days lead time"
		label2 = " "#""TRMM precipitation (mm/day) "+res
		label3 = " "#""UKMO - TRMM precipitation bias (mm/day) " + res + " " + str(lt) + " days lead time"
		bounds1 = np.linspace(0, 2, 9)
		bounds3 = np.linspace(-1.5,1.5,7)
		pl.map_composite_data(ukmopcp_mm, lats, lons, outfile_ukmo, "SIO2", 'Blues', label1, bounds1)
		pl.map_composite_data(trmmpcp_mm, lats, lons, outfile_trmm, "SIO2", 'Blues', label2, bounds1)
		pl.map_composite_data(bias, lats, lons, outfile_bias, "SIO2", 'RdBu', label3, bounds3,'bias')


#here, add the composites for the total across the entire period, once remapping is finished

#here, add the composites for the total across the entire period, once remapping is finished
for lt in [0,1,2,3,4,5,6]:
	ukmofile = ukmodir + str(lt) + "/ukmo_nwp.comp_pcp_tc." + str(lt) + "_days.072006-122016.n320.nc"
	trmmfile = trmmdir + "trmm.comp_pcp_tc.072006-122016_n320.nc"
	outfile_ukmo = "ukmo_nwp.comp_pcp_tc.072006-122016.n320." + str(lt) + "-days.png"
	outfile_trmm = "trmm.comp_pcp_tc.072006-1220160.n320.png"
	outfile_bias = "ukmo_bias.comp_pcp_tc.072006-122016.n320." + str(lt) + "-days.png"

	ffukmo = Dataset(ukmofile, 'r')
	ukmopcp = ffukmo.variables['pcp'][:]
	lons = ffukmo.variables['longitude'][:]
	lats = ffukmo.variables['latitude'][:]
	cdaysukmo = ffukmo.getncattr('contributing_days')  # contributing days is a global attribute rather than a variable

	ukmopcp_mm = ukmopcp / cdaysukmo

	fftrmm = Dataset(trmmfile, 'r')
	trmmpcp = fftrmm.variables['pcp'][:]
	cdaystrmm = fftrmm.getncattr('contributing_days')

	trmmpcp_mm = trmmpcp / cdaystrmm

	bias = ukmopcp_mm - trmmpcp_mm

	print ukmopcp_mm
	print trmmpcp_mm
	print bias

	label1 = " " #""UKMO precipitation (mm/day) 072006 - 122016 " + str(lt) + " days lead time"
	label2 = " "#""TRMM precipitation (mm/day) 072006 - 122016"
	label3 = " "#"UKMO - TRMM precipitation bias (mm/day) 072006 - 122016 " + str(lt) + " days lead time"
	bounds1 = np.linspace(0, 2, 9)
	bounds3 = np.linspace(-1.5,1.5,7)
	pl.map_composite_data(ukmopcp_mm, lats, lons, outfile_ukmo, "SIO2", 'Blues', label1, bounds1)
	pl.map_composite_data(trmmpcp_mm, lats, lons, outfile_trmm, "SIO2", 'Blues', label2, bounds1)
	pl.map_composite_data(bias, lats, lons, outfile_bias, "SIO2", 'RdBu', label3, bounds3, 'bias')