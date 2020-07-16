import picsea_library as pl
from netCDF4 import Dataset
import numpy as np



datadir = "/gws/nopw/j04/klingaman/emerton/total_gpm_imerg_precip_composites_20182019/"

year1s = [2018]
year2s = [2019]


#map the precipitation composite for one year/TC season, in mm/day
for y1, y2 in zip(year1s, year2s):
	infile = datadir + "gpm_imerg.comp_pcp_total."+str(y1)+"-"+str(y2)+".nc"
	outfile = "gpm_imerg.comp_pcp_total."+str(y1)+"-"+str(y2)+".png"
	ff = Dataset(infile, 'r')
	pcp = ff.variables['pcp'][:]
	lons = ff.variables['lon'][:]
	lats = ff.variables['lat'][:]
	cdays = ff.getncattr('contributing_days') #contributing days is a global attribute rather than a variable
	print cdays
	print np.shape(pcp)

	pcp_new=np.zeros((1800,3600))

	for x in range(3600):
		for y in range(1800):
			pcp_new[y,x] = pcp[x,y]#/cdays #divide the precip total by the number of days

	print np.max(pcp_new)
	label = " "#""TRMM precipitation (mm/day) " + str(y1) + "-" + str(y2)
	bounds = np.linspace(0, 1800, 11)
	print np.shape(pcp_new)
	pl.map_composite_data(pcp_new, lats, lons, outfile, "SIO2", 'Blues', label, bounds)


