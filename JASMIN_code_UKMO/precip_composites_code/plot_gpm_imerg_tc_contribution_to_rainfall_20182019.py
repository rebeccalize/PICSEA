import picsea_library as pl
from netCDF4 import Dataset
import numpy as np


### Script to plot the track density of a particular TC season ###

datadirib = "/gws/nopw/j04/klingaman/emerton/2018-2019_GPM-IMERG_precip_composites/"
datadirtr = "/gws/nopw/j04/klingaman/emerton/total_gpm_imerg_precip_composites_20182019/"

year1s = [2018]
year2s = [2019]


#map the precipitation composite for one year/TC season, in mm/day
for y1, y2 in zip(year1s, year2s):
	infileib = datadirib + "analysis.comp_pcp_tc."+str(y1)+"-"+str(y2)+"_v2.nc"
	infiletr = datadirtr + "gpm_imerg.comp_pcp_total."+str(y1)+"-"+str(y2)+".nc"
	print infiletr

	outfile = "analysis.gpm_imerg.tc_contrib_pcp."+str(y1)+"-"+str(y2)+"_v2.png"

	perc_cont = np.zeros((1800,3600))
	perc_cont[:,:] = np.nan

	ffib = Dataset(infileib, 'r')
	pcpib = ffib.variables['pcp'][:]
	lons = ffib.variables['lon'][:]
	lats = ffib.variables['lat'][:]
	#pcpib_overall += pcpib
	#cdaysib = ffib.getncattr('contributing_days') #contributing days is a global attribute rather than a variable
	print np.shape(pcpib)
	


	fftr = Dataset(infiletr, 'r')
	pcptr = fftr.variables['pcp'][:]
	#pcptr_overall += pcptr
	#cdaystr = fftr.getncattr('contributing_days') #contributing days is a global attribute rather than a variable
	#print cdaystr
	
	print np.shape(pcptr)

	print np.max(pcpib)
	print np.max(pcptr)

	for x in range(3600):
		for y in range(1800):
			perc_cont[y,x] = (pcpib[x,y]/pcptr[x,y])*100.
			if perc_cont[y,x] > 0.:
				print perc_cont[y,x]

	print np.nanmax(perc_cont)


	label = " "#""TC contribution to precipitation (%) " + str(y1) + "-" + str(y2)
	bounds = np.linspace(0, 100, 11)
	pl.map_composite_data(perc_cont, lats, lons, outfile, "SIO", 'BuPu', label, bounds,'perc')


	outfile_data = "gpm_imerg.analysis_tracks.perc_cont_tc_rainfall." + str(y1) + "-" + str(y2) + ".txt"
	np.savetxt(datadirib+outfile_data, perc_cont, '%.4f')




