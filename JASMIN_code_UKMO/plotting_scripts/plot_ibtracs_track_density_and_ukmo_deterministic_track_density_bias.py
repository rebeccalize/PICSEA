import picsea_library as pl
from netCDF4 import Dataset
import numpy as np


### Script to plot the track density of a particular TC season ###

datadir_nwp = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_density/"
#datadir_ib = "/gws/nopw/j04/klingaman/emerton/ibtracs_track_density/"
datadir_ib = "/gws/nopw/j04/klingaman/emerton/analysis_track_density/"

year1s=[2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015]
year2s=[2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016]

for lt in range(7): #7
	#arrays to hold data sum for all years
	dens_sum_nwp = np.zeros((72,144))
	dens_peryear_nwp = np.zeros((72,144))
	lons_overall = np.zeros((72,144))
	lats_overall = np.zeros((72,144))
	
	dens_sum_ib = np.zeros((72,144))
	dens_peryear_ib = np.zeros((72,144))
	
	for y1, y2 in zip(year1s, year2s):

		infile_nwp = datadir_nwp+str(lt)+"/"+str(y1)+str(y2)+"/ukmo_nwp.density.lt"+str(lt)+"_days."+str(y1)+str(y2)+".nc"
		outfile = "ukmo_nwp.track_density.lt"+str(lt)+"_days."+str(y1)+str(y2)+".png"
		
		ff = Dataset(infile_nwp, 'r')
		dens = ff.variables['density'][:]
		lons = ff.variables['longitude'][:]
		lats = ff.variables['latitude'][:]
		#label = "Cyclone Track Density "+str(y1)+"-"+str(y2)
		label=' '
		bounds = np.linspace(0, 20, 11)
		#pl.map_composite_data(dens, lats,lons, outfile, "SIO3", 'Greens', label,bounds,"contour")
		
		lons_overall=lons
		lats_overall=lats

		dens_sum_nwp += dens
		
		
		
		infile_ib = datadir_ib+"analysis.density."+str(y1)+str(y2)+".nc"
		outfile_ib = "analysis_track_density."+str(y1)+str(y2)+".png"
		
		ffib = Dataset(infile_ib,'r')
		dens_ib = ffib.variables['density'][:]
		label = ' '
		bounds = np.linspace(0, 20, 11)
		pl.map_composite_data(dens_ib, lats,lons, outfile_ib, "SIO3", 'Greens', label,bounds,"contour") 
		
		dens_sum_ib += dens_ib
		

		outfile_bias = "ukmo_nwp.track_density_bias.vs_analysis.lt"+str(lt)+"_days."+str(y1)+str(y2)+".png"
		bias = dens - dens_ib
				
		#bounds_bias = np.linspace(-5, 5, 11)
		bounds_bias = np.array([-5,-4,-3,-2,-1,1,2,3,4,5])
		pl.map_composite_data(bias, lats,lons, outfile_bias, "SIO3", 'BrBG', label,bounds_bias,"contour",'bias') 

	#map the track density across all the years/TC seasons
	#print dens_sum
	#print lons_overall
	#print lats_overall
	dens_peryear_nwp = dens_sum_nwp / 10
	#print np.max(dens_sum)
	#print np.max(dens_peryear)
	outfile_nwp = "ukmo_nwp.track_density_TCsperyear.lt"+str(lt)+"_days."+str(year1s[0])+str(year2s[-1])+".png"
	#label = "Cyclone Track Density "+str(year1s[0])+"-"+str(year2s[-1])+" (not including 2008-2009 & 2009-2010 seasons)"
	label = ' '
	bounds = np.linspace(0, 10, 11)
	#pl.map_composite_data(dens_peryear_nwp, lats_overall, lons_overall, outfile_nwp, "SIO3", 'Greens', label,bounds,"contour")
	
	
	dens_peryear_ib = dens_sum_ib / 10
	outfile_ib = "analysis.track_density_TCsperyear."+str(year1s[0])+str(year2s[-1])+".png"
	pl.map_composite_data(dens_peryear_ib, lats_overall, lons_overall, outfile_ib, "SIO3", 'Greens', label,bounds,"contour")
	
	
	dens_bias = dens_peryear_nwp - dens_peryear_ib
	outfile_bias = "ukmo_nwp.track_density_bias_TCsperyear.vs_analysis.lt"+str(lt)+"_days."+str(year1s[0])+str(year2s[-1])+".png"
	#bounds_bias = np.linspace(-5, 5, 11)
	bounds_bias = np.array([-5,-4,-3,-2,-1,1,2,3,4,5])
	pl.map_composite_data(dens_bias, lats,lons, outfile_bias, "SIO3", 'BrBG', label,bounds_bias,"contour",'bias') 
	
	



