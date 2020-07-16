import picsea_library as pl
from netCDF4 import Dataset
import numpy as np


### Script to plot the track density of a particular TC season ###

datadirib = "/gws/nopw/j04/klingaman/emerton/ibtracs_precip_composites/"
datadirtr = "/gws/nopw/j04/klingaman/emerton/total_trmm_precip_composites/"

year1s = [2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017]
year2s = [2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018]

#year1s = [2006]
#year2s = [2007]

# arrays to hold data sum for all years
pcpib_overall, pcptr_overall, lats_overall, lons_overall, perc_overall = (np.zeros((400, 1440)) for i in range(5))

cdays_sum = 0

#map the precipitation composite for one year/TC season, in mm/day
for y1, y2 in zip(year1s, year2s):
	infileib = datadirib + "ibtracs.comp_pcp_tc."+str(y1)+"-"+str(y2)+".nc"
	infiletr = datadirtr + "trmm.comp_pcp_total."+str(y1)+"-"+str(y2)+".nc"
	print infiletr

	outfile = "ibtracs.tc_contrib_precip."+str(y1)+"-"+str(y2)+".png"

	perc_cont = np.zeros((400, 1440))
	perc_cont[:,:] = np.nan

	ffib = Dataset(infileib, 'r')
	pcpib = ffib.variables['pcp'][:]
	lons = ffib.variables['longitude'][:]
	lats = ffib.variables['latitude'][:]
	#pcpib_overall += pcpib
	cdaysib = ffib.getncattr('contributing_days') #contributing days is a global attribute rather than a variable

	fftr = Dataset(infiletr, 'r')
	pcptr = fftr.variables['pcp'][:]
	#pcptr_overall += pcptr
	cdaystr = fftr.getncattr('contributing_days') #contributing days is a global attribute rather than a variable
	print cdaystr

	print np.max(pcpib)
	print np.max(pcptr)

	for x in range(400):
		for y in range(1440):
			perc_cont[x,y] = (pcpib[x,y]/pcptr[x,y])*100.
			pcpib_overall[x,y] += pcpib[x,y]
			pcptr_overall[x,y] += pcptr[x,y]

	print np.nanmax(perc_cont)


	label = " "#""TC contribution to precipitation (%) " + str(y1) + "-" + str(y2)
	bounds = np.linspace(0, 100, 11)
	#pl.map_composite_data(perc_cont, lats, lons, outfile, "SIO2", 'BuPu', label, bounds,'perc')
	lons_overall = lons
	lats_overall = lats
	#cdays_sum += cdays

	#outfile_data = "trmm_ibtracs_perc_cont_tc_rainfall." + str(y1) + "-" + str(y2) + ".txt"
	#np.savetxt(datadirib+outfile_data, perc_cont, '%.4f')

# map the precipitation composite for all years/TC seasons
#print pcp_sum
#print lons_overall
#print lats_overall
#print cdays_sum
for x in range(400):
	for y in range(1440):
		perc_overall[x, y] = (pcpib_overall[x, y] / pcptr_overall[x,y])*100.

#outfile = "ibtracs.tc_contrib_precip."+str(year1s[0])+"-"+str(year2s[-1])+"_v2.png"
#label = " "#""TC contribution to precipitation (%) " + str(year1s[0]) + "-" + str(year2s[-1])
#bounds = np.linspace(0,50,11)
#pl.map_composite_data(perc_overall, lats_overall, lons_overall, outfile, "SIO2", 'BuPu', label,bounds,'perc')
#outfile_data = "trmm_ibtracs_perc_cont_tc_rainfall."+str(year1s[0])+"-"+str(year2s[-1])+".txt"
#np.savetxt(datadirib+outfile_data,perc_overall,'%.4f')

for res in ['n320','n512','n768']:

	if res == 'n320':
		infileib = datadirib + "trmm.comp_pcp_tc.072006-032010_regridded_n320.nc"
		infiletr = datadirtr + "trmm.comp_pcp_tot.072006-032010_regridded_n320.nc"
		outfile = "trmm_ibtracs.perc_cont_tc_rainfall.072006-032010_n320.png"
		outfile_data = "trmm_ibtracs.perc_cont_tc_rainfall.072006-032010_n320.txt"

	elif res == 'n512':
		infileib = datadirib + "trmm.comp_pcp_tc.032010-072014_regridded_n512.nc"
		infiletr = datadirtr + "trmm.comp_pcp_tot.032010-072014_regridded_n512.nc"
		outfile = "trmm_ibtracs.perc_cont_tc_rainfall.032010-072014_n512.png"
		outfile_data = "trmm_ibtracs.perc_cont_tc_rainfall.032010-072014_n512.txt"

	elif res == 'n768':
		infileib = datadirib + "trmm.comp_pcp_tc.072014-122016_regridded_n768.nc"
		infiletr = datadirtr + "trmm.comp_pcp_tot.072014-122016_regridded_n768.nc"
		outfile = "trmm_ibtracs.perc_cont_tc_rainfall.072014-122016_n768.png"
		outfile_data = "trmm_ibtracs.perc_cont_tc_rainfall.072014-122016_n768.txt"

	print infiletr

	ffib = Dataset(infileib, 'r')
	pcpib = ffib.variables['pcp'][:]
	lons = ffib.variables['longitude'][:]
	lats = ffib.variables['latitude'][:]
	#pcpib_overall += pcpib
	cdaysib = ffib.getncattr('contributing_days') #contributing days is a global attribute rather than a variable

	fftr = Dataset(infiletr, 'r')
	pcptr = fftr.variables['pcp'][:]
	#pcptr_overall += pcptr
	cdaystr = fftr.getncattr('contributing_days') #contributing days is a global attribute rather than a variable
	print cdaystr

	print np.max(pcpib)
	print np.max(pcptr)

	perc_cont = np.zeros_like(pcptr)
	perc_cont[:, :] = np.nan

	perc_cont = (pcpib/pcptr)*100.

	print np.nanmax(perc_cont)


	label = " "#""TC contribution to precipitation (%) " + str(y1) + "-" + str(y2)
	bounds = np.linspace(0, 100, 11)
	#pl.map_composite_data(perc_cont, lats, lons, outfile, "SIO2", 'BuPu', label, bounds,'perc')

	#np.savetxt(datadirib+outfile_data, perc_cont, '%.4f')


infileib = datadirib + "trmm.comp_pcp_tc.072006-122016_n320.nc"
infiletr = datadirtr + "trmm.comp_pcp_tot.072006-122016_n320.nc"
outfile = "trmm_ibtracs.perc_cont_tc_rainfall.072006-122016_n320.png"
outfile_data = "trmm_ibtracs.perc_cont_tc_rainfall.072006-122016_n320.txt"

ffib = Dataset(infileib, 'r')
pcpib = ffib.variables['pcp'][:]
lons = ffib.variables['longitude'][:]
lats = ffib.variables['latitude'][:]
# pcpib_overall += pcpib
cdaysib = ffib.getncattr('contributing_days')  # contributing days is a global attribute rather than a variable

fftr = Dataset(infiletr, 'r')
pcptr = fftr.variables['pcp'][:]
# pcptr_overall += pcptr
cdaystr = fftr.getncattr('contributing_days')  # contributing days is a global attribute rather than a variable
print cdaystr

print np.max(pcpib)
print np.max(pcptr)

perc_cont = np.zeros_like(pcptr)
perc_cont[:, :] = np.nan

perc_cont = (pcpib / pcptr) * 100.

print np.nanmax(perc_cont)

label = " "  # ""TC contribution to precipitation (%) " + str(y1) + "-" + str(y2)
bounds = np.linspace(0, 100, 11)
#pl.map_composite_data(perc_cont, lats, lons, outfile, "SIO2", 'BuPu', label, bounds,'perc')

#np.savetxt(datadirib + outfile_data, perc_cont, '%.4f')

for season in ['DJFM','NDJFMA','MJJASO']:

	infileib = datadirib + "trmm.comp_pcp_tc.072006-122016."+season+"_n320.nc"
	infiletr = datadirtr + "trmm.comp_pcp_tot.072006-122016."+season+"_n320.nc"
	outfile = "trmm_ibtracs.perc_cont_tc_rainfall.072006-122016."+season+"_n320.png"
	outfile_data = "trmm_ibtracs.perc_cont_tc_rainfall.072006-122016."+season+"_n320.txt"

	print infiletr

	ffib = Dataset(infileib, 'r')
	pcpib = ffib.variables['pcp'][:]
	lons = ffib.variables['longitude'][:]
	lats = ffib.variables['latitude'][:]
	cdaysib = ffib.getncattr('contributing_days') #contributing days is a global attribute rather than a variable

	fftr = Dataset(infiletr, 'r')
	pcptr = fftr.variables['pcp'][:]
	cdaystr = fftr.getncattr('contributing_days') #contributing days is a global attribute rather than a variable
	print cdaystr

	print np.max(pcpib)
	print np.max(pcptr)

	perc_cont = np.zeros_like(pcptr)
	perc_cont[:, :] = np.nan

	perc_cont = (pcpib/pcptr)*100.

	print np.nanmax(perc_cont)


	label = " "#""TC contribution to precipitation (%) " + str(y1) + "-" + str(y2)
	bounds = np.linspace(0, 100, 11)
	pl.map_composite_data(perc_cont, lats, lons, outfile, "SIO2", 'BuPu', label, bounds,'perc')

	np.savetxt(datadirib+outfile_data, perc_cont, '%.4f')



