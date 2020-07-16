import picsea_library as pl
from netCDF4 import Dataset
import numpy as np


### Script to plot the track density of a particular TC season ###

datadir = "/gws/nopw/j04/klingaman/emerton/2018-2019_GPM-IMERG_precip_composites/"

year1s = [2018]
year2s = [2019]

#map the precipitation composite for one year/TC season, in mm/day
for y1, y2 in zip(year1s, year2s):

	infile = datadir + "analysis.comp_pcp_tc."+str(y1)+"-"+str(y2)+"_v2.nc"
	#infile = datadir + "2019/analysis.comp_pcp_tc.201903.nc"
	print infile
	outfile = "analysis.gpm_imerg.comp_pcp_tc."+str(y1)+"-"+str(y2)+"_v2_mmday.png"
	ff = Dataset(infile, 'r')
	pcp = ff.variables['pcp'][:]
	lons = ff.variables['lon'][:]
	lats = ff.variables['lat'][:]
	cdays = ff.getncattr('contributing_days') #contributing days is a global attribute rather than a variable
	print cdays

	pcp_new=np.zeros((1800,3600))

	for x in range(3600):
		for y in range(1800):
			pcp_new[y,x] = pcp[x,y]/cdays #divide the precip total by the number of days

	print pcp_new
	print np.max(pcp_new)
	label = " "#""TC-related precipitation (mm) " + str(y1) + "-" + str(y2)
	bounds = np.linspace(0, 10, 11)
	pl.map_composite_data(pcp_new, lats, lons, outfile, "SIO", 'Blues', label, bounds)


