import picsea_library as pl
from netCDF4 import Dataset
import numpy as np


### Script to plot the track density of a particular TC season ###

ukmodir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_precip_composites/"
ukmodensdir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_density/"
#trmmdir = "/gws/nopw/j04/klingaman/emerton/ibtracs_precip_composites/"
#densdir = "/gws/nopw/j04/klingaman/emerton/ibtracs_track_density/"
trmmdir = "/gws/nopw/j04/klingaman/emerton/analysis_trmm_tc_precip_composites/"
densdir = "/gws/nopw/j04/klingaman/emerton/analysis_track_density/"

res_list = ['n320', 'n512', 'n768']

year1s=[2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015]
year2s=[2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016]

#here, need to copy across the code from the precip composite plotting script, for plotting the bias per resolution - not regridded yet, and not calculated the TC
trmmfile = trmmdir + "trmm.comp_pcp_tc.072006-122016_track_density_2.5deg_grid.nc"
fftrmm = Dataset(trmmfile, 'r')
trmmpcp = fftrmm.variables['pcp'][:]

lons_overall = np.zeros((72,144))
lats_overall = np.zeros((72,144))
dens_sum_an = np.zeros((72,144))

for y1, y2 in zip(year1s, year2s):

	infile_ib = densdir+"analysis.density."+str(y1)+str(y2)+".nc"
	
	ffib = Dataset(infile_ib,'r')
	dens_an = ffib.variables['density'][:]
	lons = ffib.variables['longitude'][:]
	lats = ffib.variables['latitude'][:]
		
	lons_overall=lons
	lats_overall=lats
	
	dens_sum_an += dens_an
	


print np.nanmax(trmmpcp)
print np.nanmax(dens_sum_an)

trmm_obs_precip_per_tc_mmper6h = trmmpcp / dens_sum_an
trmm_obs_precip_per_tc_mmph = trmm_obs_precip_per_tc_mmper6h / 6

print np.nanmax(trmm_obs_precip_per_tc_mmper6h)
print np.nanmax(trmm_obs_precip_per_tc_mmph)

bounds1 = np.linspace(0, 50, 11)
label1=' '
outfile = "trmm_pcp_per_tc.analysis_tracks.mmph_2006-2016.png"
pl.map_composite_data(trmm_obs_precip_per_tc_mmph, lats, lons, outfile, "SIO3", 'Purples', label1, bounds1,"contour")


#get the forecast track density for the whole period, for each lead time (saved as yearly files)
fcst_dens_sum = np.zeros((7,72,144)) #lead time, lat, lon
for lt in [0,1,2,3,4,5,6]:
	for y1, y2 in zip(year1s, year2s):
		densfile = ukmodensdir + str(lt) + "/" + str(y1)+str(y2) + "/ukmo_nwp.density.lt"+str(lt)+"_days."+str(y1)+str(y2)+".nc"
		ffdens = Dataset(densfile,'r')
		dens = ffdens.variables['density'][:]
		
		fcst_dens_sum[lt,:,:] += dens[:,:]
		
print fcst_dens_sum


#here, add the composites for the total across the entire period, once remapping is finished
for lt in [0,1,2,3,4,5,6]:
	ukmofile = ukmodir + str(lt) + "/ukmo_nwp.comp_pcp_tc." + str(lt) + "_days.072006-122016.track_density_2.5deg_grid.nc"
	print ukmofile
	#trmmfile = trmmdir + "trmm.comp_pcp_tc.072006-122016_n320.nc"

	outfile_bias = "ukmo_det_pcp_bias.perTC_mmph.vs_trmm_analysis.2006-2016.track_density_2.5deg_grid." + str(lt) + "-days.png"

	ffukmo = Dataset(ukmofile, 'r')
	ukmopcp = ffukmo.variables['pcp'][:]
	lons = ffukmo.variables['longitude'][:]
	lats = ffukmo.variables['latitude'][:]
	#cdaysukmo = ffukmo.getncattr('contributing_days')  # contributing days is a global attribute rather than a variable

	#ukmopcp_mm = ukmopcp / cdaysukmo

	#fftrmm = Dataset(trmmfile, 'r')
	#trmmpcp = fftrmm.variables['pcp'][:]
	#cdaystrmm = fftrmm.getncattr('contributing_days')

	#trmmpcp_mm = trmmpcp / cdaystrmm
	
	ukmo_pcp_perTC_mmper6h = ukmopcp[:,:] / fcst_dens_sum[lt,:,:]

	bias_per_TC_mmper6h = ukmo_pcp_perTC_mmper6h - trmm_obs_precip_per_tc_mmper6h
	
	bias_per_TC_mmph = bias_per_TC_mmper6h / 6
	
	#print dens_sum_ib[dens_sum_ib < 0.5]
	
	for x in range(72):
		for y in range(144):
			if 0 < dens_sum_an[x,y] <= 0.5:
				print dens_sum_an[x,y]
				bias_per_TC_mmph[x,y] = np.ma.masked

	#print ukmopcp_mm
	#print trmmpcp_mm
	#print bias

	#label1 = " " #""UKMO precipitation (mm/day) 072006 - 122016 " + str(lt) + " days lead time"
	#label2 = " "#""TRMM precipitation (mm/day) 072006 - 122016"
	#label3 = " "#"UKMO - TRMM precipitation bias (mm/day) 072006 - 122016 " + str(lt) + " days lead time"
	label = ' '
	#bounds1 = np.linspace(0, 10, 11)
	#bounds3 = np.linspace(-20,20,11)
	bounds3 = np.array([-22,-18,-14,-10,-6,-2,2,6,10,14,18,22])
	#bounds3 = np.array([-20,-16,-12,-8,-4,-2,2,4,8,12,16,20])
	#pl.map_composite_data(ukmopcp_mm, lats, lons, outfile_ukmo, "SIO2", 'Blues', label1, bounds1)
	#pl.map_composite_data(trmmpcp_mm, lats, lons, outfile_trmm, "SIO2", 'Blues', label2, bounds1)
	pl.map_composite_data(bias_per_TC_mmph, lats, lons, outfile_bias, "SIO3", 'RdBu', label, bounds3,"contour",'bias')
	
	
	
	
	
	
