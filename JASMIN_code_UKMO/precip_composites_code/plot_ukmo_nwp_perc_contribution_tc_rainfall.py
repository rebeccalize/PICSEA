import picsea_library as pl
from netCDF4 import Dataset
import numpy as np


### Script to calculate the percentage contribution of TCs to the total rainfall in the UKMO deterministic forecasts ###

datadir_tc = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_precip_composites/"
datadir_tot = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_TOTAL_precip_composites/"

nedit 

#year1s = [2006]
#year2s = [2007]


for lt in [0,1,2,3,4,5]: #0,1,2,3,4,5,

	if lt == 6:
		year1s = [2009, 2010, 2011, 2012, 2013, 2014, 2015] #,   #, 2016, 2017
		year2s = [2010, 2011, 2012, 2013, 2014, 2015, 2016] #,  #, 2017, 2018
	
	else:
		year1s = [2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015] #,   #, 2016, 2017
		year2s = [2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016] #,  #, 2017, 2018

	#map the precipitation composite for one year/TC season, in mm/day, for this lead time
	for y1, y2 in zip(year1s, year2s):

		#Do 2009 and 2014 separately, first, so that we can do it for the two separate resolutions

		#if y1 in [2009, 2014]:
		if y1 == 2009:
			res_list = ['.n320', '.n512']
			key1='unknown' #for some reason the split resolution files have the precip variable named 'uknown' rather than 'pcp'
			key2='unknown'
		elif y1 == 2014:
			key1='unknown'
			key2='pcp'
			res_list = ['.n512', '.n768']
		else:
			res_list = ['']
			key1='pcp'
			key2='pcp'


		for res in res_list:
			infile_tc = datadir_tc + str(lt) + "/ukmo_nwp.comp_pcp_tc."+str(lt)+"-days."+str(y1)+"-"+str(y2)+res+".nc"
			print infile_tc
			infile_tot = datadir_tot + str(lt) + "/ukmo_nwp.comp_pcp_total."+str(lt)+"-days."+str(y1)+"-"+str(y2)+res+".nc"
			print infile_tot
			outfile = "ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days."+str(y1)+"-"+str(y2)+res+".png"

			ff_tc = Dataset(infile_tc, 'r')
			pcp_tc = ff_tc.variables[key1][:]
			lons = ff_tc.variables['longitude'][:]
			lats = ff_tc.variables['latitude'][:]
			#pcp_tc_overall += pcp_tc
			cdays_tc = ff_tc.getncattr('contributing_days') #contributing days is a global attribute rather than a variable

			ff_tot = Dataset(infile_tot, 'r')
			pcp_tot = ff_tot.variables[key2][:]
			#pcp_tot_overall += pcp_tot
			cdays_tot = ff_tot.getncattr('contributing_days') #contributing days is a global attribute rather than a variable
			print cdays_tot

			perc_cont = np.zeros_like(pcp_tot)
			perc_cont[:,:] = np.nan
			print len(perc_cont[:,0])
			print len(perc_cont[0,:])

			print np.max(pcp_tc)
			print np.max(pcp_tot)

			for x in range(len(perc_cont[:,0])):
				for y in range(len(perc_cont[0,:])):
					perc_cont[x,y] = (pcp_tc[x,y]/pcp_tot[x,y])*100.

			print np.nanmax(perc_cont)


			label = " "#""TC contribution to precipitation (%) " + str(y1) + "-" + str(y2)
			bounds = np.linspace(0, 100, 11)
			#pl.map_composite_data(perc_cont, lats, lons, outfile, "SIO2", 'BuPu', label, bounds,'perc')
			lons_overall = lons
			lats_overall = lats
			#cdays_sum += cdays

			outfile_data = "/ukmo_nwp_perc_cont_tc_rainfall."+str(lt)+"-days." + str(y1) + "-" + str(y2) + res+".txt"
			#np.savetxt(datadir_tc+str(lt) + outfile_data, perc_cont, '%.4f')


	#now for each resolution
	for res in ['n320','n512','n768']:

		if res == 'n320':
			infile_tc = datadir_tc + str(lt) + "/ukmo_nwp.comp_pcp_tc."+str(lt)+"_days.072006-032010.n320.nc"
			infile_tot = datadir_tot + str(lt) + "/ukmo_nwp.comp_pcp_tot."+str(lt)+"_days.072006-032010.n320.nc"
			outfile = "ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.072006-032010_n320.png"
			outfile_data = "/ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.072006-032010_n320.txt"

		elif res == 'n512':
			infile_tc = datadir_tc + str(lt) + "/ukmo_nwp.comp_pcp_tc."+str(lt)+"_days.032010-072014.n512.nc"
			infile_tot = datadir_tot + str(lt) + "/ukmo_nwp.comp_pcp_tot."+str(lt)+"_days.032010-072014.n512.nc"
			outfile = "ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.032010-072014_n512.png"
			outfile_data = "/ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.032010-072014_n512.txt"

		elif res == 'n768':
			infile_tc = datadir_tc + str(lt) + "/ukmo_nwp.comp_pcp_tc."+str(lt)+"_days.072014-122016.n768.nc"
			infile_tot = datadir_tot + str(lt) + "/ukmo_nwp.comp_pcp_tot."+str(lt)+"_days.072014-122016.n768.nc"
			outfile = "ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.072014-122016_n768.png"
			outfile_data = "/ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.072014-122016_n768.txt"

		print infile_tc
		print infile_tot

		ff_tc = Dataset(infile_tc, 'r')
		pcp_tc = ff_tc.variables['pcp'][:]
		lons = ff_tc.variables['longitude'][:]
		lats = ff_tc.variables['latitude'][:]
		#pcp_tc_overall += pcp_tc
		cdays_tc = ff_tc.getncattr('contributing_days') #contributing days is a global attribute rather than a variable

		ff_tot = Dataset(infile_tot, 'r')
		pcp_tot = ff_tot.variables['pcp'][:]
		#pcp_tot_overall += pcp_tot
		cdays_tot = ff_tot.getncattr('contributing_days') #contributing days is a global attribute rather than a variable
		print cdays_tot

		print np.max(pcp_tc)
		print np.max(pcp_tot)

		perc_cont = np.zeros_like(pcp_tot)
		perc_cont[:, :] = np.nan

		perc_cont = (pcp_tc/pcp_tot)*100.

		print np.nanmax(perc_cont)


		label = " "#""TC contribution to precipitation (%) " + str(y1) + "-" + str(y2)
		bounds = np.linspace(0, 100, 11)
		#pl.map_composite_data(perc_cont, lats, lons, outfile, "SIO2", 'BuPu', label, bounds,'perc')

		#np.savetxt(datadir_tc+str(lt) + outfile_data, perc_cont, '%.4f')


	infile_tc = datadir_tc + str(lt) + "/ukmo_nwp.comp_pcp_tc."+str(lt)+"_days.072006-122016.n320.nc"
	infile_tot = datadir_tot + str(lt) + "/ukmo_nwp.comp_pcp_tot."+str(lt)+"_days.072006-122016.n320.nc"
	outfile = "ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.072006-122016_n320.png"
	outfile_data = "/ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.072006-122016_n320.txt"



	ff_tc = Dataset(infile_tc, 'r')
	pcp_tc = ff_tc.variables['pcp'][:]
	lons = ff_tc.variables['longitude'][:]
	lats = ff_tc.variables['latitude'][:]
	# pcp_tc_overall += pcp_tc
	cdays_tc = ff_tc.getncattr('contributing_days')  # contributing days is a global attribute rather than a variable

	ff_tot = Dataset(infile_tot, 'r')
	pcp_tot = ff_tot.variables['pcp'][:]
	# pcp_tot_overall += pcp_tot
	cdays_tot = ff_tot.getncattr('contributing_days')  # contributing days is a global attribute rather than a variable
	print cdays_tot

	print np.max(pcp_tc)
	print np.max(pcp_tot)

	perc_cont = np.zeros_like(pcp_tot)
	perc_cont[:, :] = np.nan

	perc_cont = (pcp_tc / pcp_tot) * 100.

	print np.nanmax(perc_cont)

	label = " "  # ""TC contribution to precipitation (%) " + str(y1) + "-" + str(y2)
	bounds = np.linspace(0, 50, 11)
	#pl.map_composite_data(perc_cont, lats, lons, outfile, "SIO2", 'BuPu', label, bounds)

	#np.savetxt(datadir_tc + str(lt) + outfile_data, perc_cont, '%.4f')

	for season in ['DJFM','NDJFMA','MJJASO']:


		infile_tc = datadir_tc + str(lt) + "/ukmo_nwp.comp_pcp_tc."+str(lt)+"_days.072006-122016."+season+".n320.nc"
		infile_tot = datadir_tot + str(lt) + "/ukmo_nwp.comp_pcp_tot."+str(lt)+"_days.072006-122016."+season+".n320.nc"
		outfile = "ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.072006-122016."+season+".n320.png"
		outfile_data = "/ukmo_nwp.perc_cont_tc_rainfall."+str(lt)+"-days.072006-122016"+season+".n320.txt"

		print infile_tc
		print infile_tot

		ff_tc = Dataset(infile_tc, 'r')
		pcp_tc = ff_tc.variables['pcp'][:]
		lons = ff_tc.variables['longitude'][:]
		lats = ff_tc.variables['latitude'][:]
		cdays_tc = ff_tc.getncattr('contributing_days') #contributing days is a global attribute rather than a variable

		ff_tot = Dataset(infile_tot, 'r')
		pcp_tot = ff_tot.variables['pcp'][:]
		cdays_tot = ff_tot.getncattr('contributing_days') #contributing days is a global attribute rather than a variable
		print cdays_tot

		print np.max(pcp_tc)
		print np.max(pcp_tot)

		perc_cont = np.zeros_like(pcp_tot)
		perc_cont[:, :] = np.nan

		perc_cont = (pcp_tc/pcp_tot)*100.

		print np.nanmax(perc_cont)


		label = " "#""TC contribution to precipitation (%) " + str(y1) + "-" + str(y2)
		bounds = np.linspace(0, 100, 11)
		pl.map_composite_data(perc_cont, lats, lons, outfile, "SIO2", 'BuPu', label, bounds,'perc')

		np.savetxt(datadir_tc+str(lt) + outfile_data, perc_cont, '%.4f')

