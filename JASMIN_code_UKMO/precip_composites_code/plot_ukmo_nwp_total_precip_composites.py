import picsea_library as pl
from netCDF4 import Dataset
import numpy as np



datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_TOTAL_precip_composites/"


#map the precipitation composite for one year/TC season, in mm/day
for lt in [6]:

	y1s_res = [2009, 2014]
	y2s_res = [2010, 2015]

	if lt == 6:
		year1s = [2010, 2011, 2012, 2013, 2015]
		year2s = [2011, 2012, 2013, 2014, 2016]
		resolutions = ['n512', 'n768']
	else:
		year1s = [2006, 2007, 2008, 2010, 2011, 2012, 2013, 2015]
		year2s = [2007, 2008, 2009, 2011, 2012, 2013, 2014, 2016]
		resolutions = ['n320', 'n512', 'n768']

	# arrays to hold data sum for all years, for each resolution
	pcp_sum_dict, lons_overall_dict, lats_overall_dict, cdays_sum_dict = ({} for i in range(4))

	for y1, y2 in zip(year1s, year2s):
		if y1 in [2006,2007,2008]:
			xr=481
			yr=640
			res='n320'
		elif y1 in [2010, 2011, 2012, 2013]:
			xr=769
			yr=1024
			res='n512'
		elif y1 in [2015]:
			xr=1152
			yr=1536
			res='n768'


		infile = datadir + str(lt) + "/ukmo_nwp.comp_pcp_total."+str(lt)+"-days."+str(y1)+"-"+str(y2)+".nc"
		print infile
		outfile = "ukmo_nwp.comp_pcp_total."+str(lt)+"-days."+str(y1)+"-"+str(y2)+"_mmday.png"
		ff = Dataset(infile, 'r')
		lons = ff.variables['longitude'][:]
		lats = ff.variables['latitude'][:]
		cdays = ff.getncattr('contributing_days')  # contributing days is a global attribute rather than a variable
		pcp_tot = ff.variables['pcp'][:]
		#print pcp_tot

		if res in pcp_sum_dict.keys():
			pcp_sum_dict[res] += pcp_tot
			cdays_sum_dict[res] += cdays
		else:
			print 'here'
			pcp_sum_dict[res] = pcp_tot
			cdays_sum_dict[res] = cdays
			lons_overall_dict[res] = lons
			lats_overall_dict[res] = lats

		#print pcp_sum_dict
		pcp=np.zeros((xr,yr))
		for x in range(xr):
			for y in range(yr):
				pcp[x,y] = pcp_tot[x,y]/cdays #divide the precip total by the number of days
		print np.max(pcp)

		label = res #""TRMM precipitation (mm/day) " + str(y1) + "-" + str(y2)
		bounds = np.linspace(0, 16, 9)
		pl.map_composite_data(pcp, lats, lons, outfile, "SIO2", 'Blues', label, bounds)

	if lt == 6:
		continue
	else:
		infile = datadir + str(lt) + "/ukmo_nwp.comp_pcp_total."+str(lt)+"-days.2009-2010.n320.nc"
		print infile
		outfile = "ukmo_nwp.comp_pcp_total."+str(lt)+"-days.2009-2010.n320_mmday.png"
		ff = Dataset(infile, 'r')
		lons = ff.variables['longitude'][:]
		lats = ff.variables['latitude'][:]
		cdays = ff.getncattr('contributing_days')  # contributing days is a global attribute rather than a variable
		pcp_tot = ff.variables['unknown'][:]
		print pcp_tot
		pcp_sum_dict['n320'] += pcp_tot
		cdays_sum_dict['n320'] += cdays
		pcp=np.zeros((481,640))
		for x in range(481):
			for y in range(640):
				pcp[x,y] = pcp_tot[x,y]/cdays
		label = 'n320'
		bounds=np.linspace(0, 16, 9)
		pl.map_composite_data(pcp, lats, lons, outfile, "SIO2", 'Blues', label, bounds)


	infile = datadir + str(lt) + "/ukmo_nwp.comp_pcp_total."+str(lt)+"-days.2009-2010.n512.nc"
	outfile = "ukmo_nwp.comp_pcp_total."+str(lt)+"-days.2009-2010.n512_mmday.png"
	ff = Dataset(infile, 'r')
	lons = ff.variables['longitude'][:]
	lats = ff.variables['latitude'][:]
	cdays = ff.getncattr('contributing_days')  # contributing days is a global attribute rather than a variable
	pcp_tot = ff.variables['unknown'][:]
	pcp_sum_dict['n512'] += pcp_tot
	cdays_sum_dict['n512'] += cdays
	pcp=np.zeros((769,1024))
	for x in range(769):
		for y in range(1024):
			pcp[x,y] = pcp_tot[x,y]/cdays
	label = 'n512'
	bounds=np.linspace(0, 16, 9)
	pl.map_composite_data(pcp, lats, lons, outfile, "SIO2", 'Blues', label, bounds)


	infile = datadir + str(lt) + "/ukmo_nwp.comp_pcp_total."+str(lt)+"-days.2014-2015.n512.nc"
	outfile = "ukmo_nwp.comp_pcp_total."+str(lt)+"-days.2014-2015.n512_mmday.png"
	ff = Dataset(infile, 'r')
	lons = ff.variables['longitude'][:]
	lats = ff.variables['latitude'][:]
	cdays = ff.getncattr('contributing_days')  # contributing days is a global attribute rather than a variable
	pcp_tot = ff.variables['pcp'][:]
	pcp_sum_dict['n512'] += pcp_tot
	cdays_sum_dict['n512'] += cdays
	pcp=np.zeros((769,1024))
	for x in range(769):
		for y in range(1024):
			pcp[x,y] = pcp_tot[x,y]/cdays
	label = 'n512'
	bounds=np.linspace(0, 16, 9)
	pl.map_composite_data(pcp, lats, lons, outfile, "SIO2", 'Blues', label, bounds)


	infile = datadir + str(lt) + "/ukmo_nwp.comp_pcp_total."+str(lt)+"-days.2014-2015.n768.nc"
	outfile = "ukmo_nwp.comp_pcp_total."+str(lt)+"-days.2014-2015.n768_mmday.png"
	ff = Dataset(infile, 'r')
	lons = ff.variables['longitude'][:]
	lats = ff.variables['latitude'][:]
	cdays = ff.getncattr('contributing_days')  # contributing days is a global attribute rather than a variable
	pcp_tot = ff.variables['pcp'][:]
	pcp_sum_dict['n768'] += pcp_tot
	cdays_sum_dict['n768'] += cdays
	pcp=np.zeros((1152,1536))
	for x in range(1152):
		for y in range(1536):
			pcp[x,y] = pcp_tot[x,y]/cdays
	label = 'n768'
	bounds=np.linspace(0, 16, 9)
	pl.map_composite_data(pcp, lats, lons, outfile, "SIO2", 'Blues', label, bounds)


	#now do the resolution change years separately
	#check the above works first!

	# map the precipitation composite for all years for each resolution
	#print pcp_sum
	#print lons_overall
	#print lats_overall

	#print pcp_sum_dict
	for res in resolutions:
		print res
		cdays_sum=cdays_sum_dict[res]
		pcp_sum=pcp_sum_dict[res]
		
		print pcp_sum
		print cdays
		
		if res == 'n320':
			xr = 481
			yr = 640
		elif res == 'n512':
			xr=769
			yr=1024
		elif res == 'n768':
			xr=1152
			yr=1536

		for x in range(xr):
			for y in range(yr):
				pcp_sum[x, y] = pcp_sum[x, y] / cdays_sum

		print "keys:", pcp_sum_dict.keys()
		print pcp_sum_dict

		outfile = "ukmo_nwp.comp_pcp_total."+str(lt)+"-days."+res+"_mmday.png"
		label = res#""TRMM precipitation (mm/day) " + str(year1s[0]) + "-" + str(year2s[-1])
		bounds = np.linspace(0,16,9)
		pl.map_composite_data(pcp_sum, lats_overall_dict[res], lons_overall_dict[res], outfile, "SIO2", 'Blues', label,bounds)