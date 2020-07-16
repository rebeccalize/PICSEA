import picsea_library as pl
from netCDF4 import Dataset
import numpy as np



datadir = "/gws/nopw/j04/klingaman/emerton/total_trmm_precip_composites/"

year1s = [2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017]
year2s = [2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018]

#year1s = [2006]
#year2s = [2007]

# arrays to hold data sum for all years
pcp_sum, lons_overall, lats_overall = (np.zeros((400, 1440)) for i in range(3))

cdays_sum = 0

#map the precipitation composite for one year/TC season, in mm/day
for y1, y2 in zip(year1s, year2s):
	infile = datadir + "trmm.comp_pcp_total."+str(y1)+"-"+str(y2)+".nc"
	outfile = "trmm.comp_pcp_total."+str(y1)+"-"+str(y2)+"_mmday.png"
	ff = Dataset(infile, 'r')
	pcp = ff.variables['pcp'][:]
	pcp_sum += pcp
	lons = ff.variables['longitude'][:]
	lats = ff.variables['latitude'][:]
	cdays = ff.getncattr('contributing_days') #contributing days is a global attribute rather than a variable
	print cdays
	for x in range(400):
		for y in range(1440):
			pcp[x,y] = pcp[x,y]/cdays #divide the precip total by the number of days
	print np.max(pcp)
	label = " "#""TRMM precipitation (mm/day) " + str(y1) + "-" + str(y2)
	bounds = np.linspace(0, 16, 9)
	pl.map_composite_data(pcp, lats, lons, outfile, "SIO2", 'Blues', label, bounds)
	lons_overall = lons
	lats_overall = lats
	cdays_sum += cdays

# map the precipitation composite for all years/TC seasons
#print pcp_sum
#print lons_overall
#print lats_overall
print cdays_sum
for x in range(400):
	for y in range(1440):
		pcp_sum[x, y] = pcp_sum[x, y] / cdays_sum

outfile = "trmm.comp_pcp_total."+str(year1s[0])+"-"+str(year2s[-1])+"_mmday.png"
label = " "#""TRMM precipitation (mm/day) " + str(year1s[0]) + "-" + str(year2s[-1])
bounds = np.linspace(0,16,9)
pl.map_composite_data(pcp_sum, lats_overall, lons_overall, outfile, "SIO2", 'Blues', label,bounds)