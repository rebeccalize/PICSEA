import picsea_library as pl
from netCDF4 import Dataset
import numpy as np


### Script to plot the track density of a particular TC season ###

totalfile = "/gws/nopw/j04/klingaman/emerton/total_trmm_precip_composites/trmm.comp_pcp_tot.072006-122016.NDJFMA_n320.nc"
file20182019 = "/gws/nopw/j04/klingaman/emerton/total_gpm_imerg_precip_composites_20182019/gpm_imerg.comp_pcp_total.2018-2019_regridded_n320.nc"

# arrays to hold data sum for all years
#pcp_sum, lons_overall, lats_overall = (np.zeros((400, 1440)) for i in range(3))

#cdays_sum = 0

#map the precipitation composite for one year/TC season, in mm/day
outfile_bias = "bias_20182019_season_compared_to_2006-2016_PERC.png"


fftotal = Dataset(totalfile, 'r')
totalpcp = fftotal.variables['pcp'][:]
lons = fftotal.variables['longitude'][:]
lats = fftotal.variables['latitude'][:]
cdaystotal = fftotal.getncattr('contributing_days') #contributing days is a global attribute rather than a variable

#to get mm/day
totalpcp_mm = totalpcp / cdaystotal

#to get average total per year:
#avgtotpcp = totalpcp / 10

ff20182019 = Dataset(file20182019, 'r')
pcp20182019 = ff20182019.variables['pcp'][:]
print np.shape(pcp20182019)
pcp20182019_new=np.zeros((481,640))
for x in range(640):
	for y in range(481):
		pcp20182019_new[y,x] = pcp20182019[x,y]

cdays20182019 = ff20182019.getncattr('contributing_days')

pcp20182019_mm = pcp20182019_new / cdays20182019

#bias = avgtotpcp - pcp20182019_new

pcp20182019_ma = np.ma.masked_where(pcp20182019_new == 0, pcp20182019_new)
totalpcp_ma = np.ma.masked_where(totalpcp == 0, totalpcp)

bias = pcp20182019_ma - (totalpcp_ma/10)
perc = (bias / (totalpcp_ma/10)) *100

print bias

label3 = " "#""UKMO - TRMM precipitation bias (mm/day) " + res + " " + str(lt) + " days lead time"

bounds3 = np.linspace(-100,100,21)

pl.map_composite_data(perc, lats, lons, outfile_bias, "SIO2", 'RdBu', label3, bounds3,'bias')


