import sys
#sys.path.append('/usr/lib/python2.7/site-packages/mpl_toolkits/')
#sys.path.append('/usr/lib/python2.7/site-packages/')
print '\n'.join(sys.path)
import picsea_library as pl
import fnmatch
import os
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap, cm
#import mpl_toolkits
#mpl_toolkits.__path__.append('/usr/lib/python2.7/site-packages/mpl_toolkits/')
from mpl_toolkits.basemap import Basemap,cm
import matplotlib.cm as cm
import numpy as np
from matplotlib.patches import Polygon
from netCDF4 import Dataset, num2date, netcdftime
import datetime
import matplotlib.path as mpath
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm, Normalize
import matplotlib.colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

class MidpointNormalize(matplotlib.colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

	#def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		#x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		#return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


fcst_type = 'mean'
error_type = 'track'

#every_TC_database[0] = "TC_season"
#every_TC_database[1] = "track_number"
#every_TC_database[2] = "obs_date"
#every_TC_database[3] = "obs_lon"
#every_TC_database[4] = "obs_lat"
#every_TC_database[5] = "obs_wind_kmh"
#every_TC_database[6] = "obs_mslp"
#every_TC_database[7] = "fcst_date"
#every_TC_database[8] = "fcst_lon"
#every_TC_database[9]  = "fcst_lat"
#every_TC_database[10] = "track_error"
#every_TC_database[11] = "fcst_wind_kmh"
#every_TC_database[12] = "fcst_mslp"
#every_TC_database[13] = "MJO_phase"
#every_TC_database[14] = "MJO_amplitude"

for lt in [0,1,2,3,4,5,6,7]:

	data_file = "2010_2020.ukmo_"+fcst_type+".database_of_each_individual_forecast.obs_and_fcst.location_intensity_trackerror_MJOphase.lead_time"+str(lt)+".using_ibtracs.txt"

	data=np.genfromtxt(data_file, dtype=float,skip_header=1)
	
	obs_lon=data[:,3]
	obs_lat=data[:,4]
	
	track_error=data[:,10]
	
	obs_wind=data[:,5]
	obs_mslp=data[:,6]
	
	fcst_wind=data[:,11]
	fcst_mslp=data[:,12]
	
	MJO_phase=data[:,13]
	MJO_amp=data[:,14]

	############################################
	
	wind_bias = np.zeros(len(fcst_wind))
	
	mslp_bias = np.zeros(len(fcst_mslp))
	
	for i in range(len(fcst_mslp)):
	
		if fcst_mslp[i] > 10000 or obs_mslp[i] > 10000:
		
			mslp_bias[i] = np.nan
		else:
			mslp_bias[i] = fcst_mslp[i] - obs_mslp[i]
			
		if fcst_wind[i] > 10000 or obs_wind[i] > 10000:
		
			wind_bias[i] = np.nan
			
		else:
			wind_bias[i] = fcst_wind[i] - obs_wind[i]
		

	lat1=0
	lat2=-50
	lon1=20
	lon2=110
	
	#track errors
	
	if error_type == 'track':

		fig = plt.figure(figsize=(6, 3))
		ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

		m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')


		#m.fillcontinents(color='white')
		m.drawcoastlines(linewidth=0.4, color='k')
		m.drawcountries(linewidth=0.4, color='k')
					
		cmap = 'PuRd'
		cmap = plt.get_cmap(cmap)
		colors = cmap(np.linspace(0.1, 1, cmap.N))

		cmap2 = LinearSegmentedColormap.from_list(cmap, colors)

		bounds=np.linspace(0,2000,21)
		norm = BoundaryNorm(bounds, ncolors=cmap2.N, clip=True)



		x,y = m(obs_lon, obs_lat)

		im = m.contourf(x,y,track_error, cmap=cmap2, levels=bounds, tri=True, extend="max")
		m.scatter(x,y,c='navy',s=0.2)

		cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="max", ticks=bounds[0::2]) #
		cbar.ax.set_yticklabels([int(i) for i in bounds[0::2]])


		cbar.ax.tick_params(labelsize=14)  # colorbar font size
		plt.text(0.01, 0.03, 'Track Errors (km)\n'+str(lt)+' Days Ahead', transform=ax.transAxes, fontsize=12)
		plt.savefig('ukmo_'+fcst_type+'.spatial_track_errors_map.lead_time_'+str(lt)+'.2010-2020.contour_with_dots.png', bbox_inches='tight', dpi=400)
		plt.close()
		print 'ukmo_'+fcst_type+'.spatial_track_errors_map.lead_time_'+str(lt)+'.2010-2020.png'
	
		#########################################################################
	
		#track errors

		fig = plt.figure(figsize=(6, 3))
		ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

		m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')


		#m.fillcontinents(color='white')
		m.drawcoastlines(linewidth=0.4, color='k')
		m.drawcountries(linewidth=0.4, color='k')
					
		cmap = 'PuRd'
		cmap = plt.get_cmap(cmap)
		colors = cmap(np.linspace(0.1, 1, cmap.N))

		cmap2 = LinearSegmentedColormap.from_list(cmap, colors)

		bounds=np.linspace(0,2000,21)
		norm = BoundaryNorm(bounds, ncolors=cmap2.N, clip=True)



		x,y = m(obs_lon, obs_lat)

		im = m.scatter(x, y, c=track_error, s= 3, cmap=cmap2, norm=norm)

		cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="max", ticks=bounds[0::2]) #
		cbar.ax.set_yticklabels([int(i) for i in bounds[0::2]])


		cbar.ax.tick_params(labelsize=14)  # colorbar font size
		plt.text(0.01, 0.03, 'Track Errors (km)\n'+str(lt)+' Days Ahead', transform=ax.transAxes, fontsize=12)
		plt.savefig('ukmo_'+fcst_type+'.spatial_track_errors_map.lead_time_'+str(lt)+'.2010-2020.dots.png', bbox_inches='tight', dpi=400)
		plt.close()
	

		#########################################################################

	
	elif error_type == 'mslp':
	
		print 'dots'
	
		fig = plt.figure(figsize=(6, 3))
		ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

		m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')


		#m.fillcontinents(color='white')
		m.drawcoastlines(linewidth=0.4, color='k')
		m.drawcountries(linewidth=0.4, color='k')
					
		cmap = 'PiYG'
		cmap = plt.get_cmap(cmap)
		colors = cmap(np.linspace(0, 1, cmap.N))

		cmap2 = LinearSegmentedColormap.from_list(cmap, colors)

		bounds=np.linspace(-50,50,11)
		#norm = BoundaryNorm(bounds, ncolors=cmap2.N, clip=True)
		norm = MidpointNormalize(midpoint=0., vmin=bounds[0], vmax=bounds[-1])

		print np.nanmax(wind_bias)
		print np.nanmin(wind_bias)
	

		x,y = m(obs_lon, obs_lat)

		im = m.scatter(x, y, c=mslp_bias, edgecolor='k',linewidth=0.2,s= 8, cmap=cmap2, norm=norm)

		cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="both") #
		#cbar.ax.set_yticklabels([int(i) for i in bounds])


		cbar.ax.tick_params(labelsize=14)  # colorbar font size
		plt.text(0.01, 0.03, 'MSLP Bias (hPa)\n'+str(lt)+' Days Ahead', transform=ax.transAxes, fontsize=12)
		plt.savefig('ukmo_'+fcst_type+'.spatial_mslp_biases_map.lead_time_'+str(lt)+'.2010-2020.dots.png', bbox_inches='tight', dpi=400)
		plt.close()
	
		#########################################################################
		
		print 'contour'

		fig = plt.figure(figsize=(6, 3))
		ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

		m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')


		#m.fillcontinents(color='white')
		m.drawcoastlines(linewidth=0.4, color='k')
		m.drawcountries(linewidth=0.4, color='k')
					
		cmap = 'PiYG' #'PiYG'
		cmap = plt.get_cmap(cmap)
		colors = cmap(np.linspace(0, 1, cmap.N))

		cmap2 = LinearSegmentedColormap.from_list(cmap, colors)

		bounds=np.linspace(-50,50,11)
		#norm = BoundaryNorm(bounds, ncolors=cmap2.N, clip=True)
		norm = MidpointNormalize(midpoint=0., vmin=bounds[0], vmax=bounds[-1])
	
		print 'attempting to remove nans'
		
		new_mslp_bias = []
		new_obs_lons = []
		new_obs_lats = []


		#need to remove nans for contourf to work
		for q in range(len(mslp_bias)):
		
			if not np.isnan(mslp_bias[q]):
			
				new_mslp_bias = np.append(new_mslp_bias, mslp_bias[q])
				new_obs_lons = np.append(new_obs_lons, obs_lon[q])
				new_obs_lats = np.append(new_obs_lats, obs_lat[q])
				
				
			
				
			
		x,y = m(new_obs_lons, new_obs_lats)
		
		print len(x)
		print len(y)


		print new_mslp_bias 
		
		print np.shape(new_mslp_bias)

		
		im = m.contourf(x,y,new_mslp_bias, cmap=cmap, levels=bounds, tri=True, extend="both")
		
		m.scatter(x,y,c='navy',s=0.2)

		cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="both", ticks=bounds) #
		#cbar.ax.set_yticklabels([int(i) for i in bounds])


		cbar.ax.tick_params(labelsize=14)  # colorbar font size
		plt.text(0.01, 0.03, 'MSLP Error (hPa)\n'+str(lt)+' Days Ahead', transform=ax.transAxes, fontsize=12)
		plt.savefig('ukmo_'+fcst_type+'.spatial_mslp_biases_map.lead_time_'+str(lt)+'.2010-2020.contours_with_dots.png', bbox_inches='tight', dpi=400)
		plt.close()
		
	elif error_type == 'wind':
	
		print "wind"
	
		print 'dots'
	
		fig = plt.figure(figsize=(6, 3))
		ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

		m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')


		#m.fillcontinents(color='white')
		m.drawcoastlines(linewidth=0.4, color='k')
		m.drawcountries(linewidth=0.4, color='k')
					
		cmap = 'PuOr'
		cmap = plt.get_cmap(cmap)
		colors = cmap(np.linspace(0, 1, cmap.N))

		cmap2 = LinearSegmentedColormap.from_list(cmap, colors)

		bounds=np.linspace(-80,80,9)
		#norm = BoundaryNorm(bounds, ncolors=cmap2.N, clip=True)
		norm = MidpointNormalize(midpoint=0., vmin=bounds[0], vmax=bounds[-1])

		print np.nanmax(mslp_bias)
		print np.nanmin(mslp_bias)
	

		x,y = m(obs_lon, obs_lat)

		im = m.scatter(x, y, c=wind_bias, edgecolor='k',linewidth=0.2,s= 8, cmap=cmap2, norm=norm)

		cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="both") #
		#cbar.ax.set_yticklabels([int(i) for i in bounds])


		cbar.ax.tick_params(labelsize=14)  # colorbar font size
		plt.text(0.01, 0.03, 'Wind Speed Bias (km/h)\n'+str(lt)+' Days Ahead', transform=ax.transAxes, fontsize=12)
		plt.savefig('ukmo_'+fcst_type+'.spatial_wind_speed_biases_map.lead_time_'+str(lt)+'.2010-2020.dots.png', bbox_inches='tight', dpi=400)
		plt.close()
	
		#########################################################################
		
		print 'contour'

		fig = plt.figure(figsize=(6, 3))
		ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

		m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')


		#m.fillcontinents(color='white')
		m.drawcoastlines(linewidth=0.4, color='k')
		m.drawcountries(linewidth=0.4, color='k')
					
		cmap = 'PuOr' #'PiYG'
		cmap = plt.get_cmap(cmap)
		colors = cmap(np.linspace(0, 1, cmap.N))

		cmap2 = LinearSegmentedColormap.from_list(cmap, colors)

		bounds=np.linspace(-80,80,9)
		#norm = BoundaryNorm(bounds, ncolors=cmap2.N, clip=True)
		norm = MidpointNormalize(midpoint=0., vmin=bounds[0], vmax=bounds[-1])
	
		print 'attempting to remove nans'
		
		new_wind_bias = []
		new_obs_lons = []
		new_obs_lats = []


		#need to remove nans for contourf to work
		for q in range(len(wind_bias)):
		
			if not np.isnan(wind_bias[q]):
			
				new_wind_bias = np.append(new_wind_bias, wind_bias[q])
				new_obs_lons = np.append(new_obs_lons, obs_lon[q])
				new_obs_lats = np.append(new_obs_lats, obs_lat[q])
				
				
			
				
			
		x,y = m(new_obs_lons, new_obs_lats)
		
		print len(x)
		print len(y)


		print new_wind_bias 
		
		print np.shape(new_wind_bias)

		
		im = m.contourf(x,y,new_wind_bias, cmap=cmap, levels=bounds, tri=True, extend="both")
		
		m.scatter(x,y,c='navy',s=0.2)

		cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="both") #
		#cbar.ax.set_yticklabels([int(i) for i in bounds])


		cbar.ax.tick_params(labelsize=14)  # colorbar font size
		plt.text(0.01, 0.03, 'Wind Speed Bias (km/h)\n'+str(lt)+' Days Ahead', transform=ax.transAxes, fontsize=12)
		plt.savefig('ukmo_'+fcst_type+'.spatial_wind_speed_biases_map.lead_time_'+str(lt)+'.2010-2020.contours_with_dots.png', bbox_inches='tight', dpi=400)
		plt.close()

