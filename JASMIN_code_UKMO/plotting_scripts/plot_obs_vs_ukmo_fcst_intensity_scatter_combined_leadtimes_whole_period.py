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
#error_type = 'mslp'

if fcst_type == 'hres':
	fcst_type_label = 'HRES'
	
elif fcst_type == 'mean':
	fcst_type_label = 'Ensemble Mean'

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

for error_type in ['mslp','wind']:

	c='darkslategray'
	
	colours=['#EF476F', '#F78C6B', '#FFD166', '#06D6A0','#34E4EA','#118AB2','#073B4C','#9B7EDE']

	
	alphas = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
	
	p1 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=c, marker='o', alpha=alphas[0], label='0 Days') 
	p2 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=c, marker='o', alpha=alphas[1],label='1 Day') 
	p3 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=c, marker='o', alpha=alphas[2],label='2 Days') 
	p4 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=c, marker='o', alpha=alphas[3],label='3 Days')
	p5 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=c, marker='o', alpha=alphas[4],label='4 Days')
	p6 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=c, marker='o', alpha=alphas[5],label='5 Days')
	p7 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=c, marker='o', alpha=alphas[6],label='6 Days')
	p8 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=c, marker='o', alpha=alphas[7],label='7 Days')
	p9 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor='silver', marker='o', label='None')
	
	p1 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=colours[0], marker='o', label='MJO 1') 
	p2 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=colours[1], marker='o', label='MJO 2') 
	p3 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=colours[2], marker='o', label='MJO 3') 
	p4 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=colours[3], marker='o', label='MJO 4')
	p5 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=colours[4], marker='o', label='MJO 5')
	p6 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=colours[5], marker='o', label='MJO 6')
	p7 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=colours[6], marker='o', label='MJO 7')
	p8 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor=colours[7], marker='o', label='MJO 8')
	p9 = plt.Line2D((0, 1), (0, 0), color='w', markerfacecolor='silver', marker='o', label='None')
	
	if error_type == 'mslp':
	
		fig, ax = plt.subplots()
		fig.set_size_inches(6,6)
		
		x = []
		y = []
		
		for lt in [7,6,5,4,3,2,1,0]:	
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

		
			for i in range(len(fcst_wind)):
				
				if fcst_mslp[i] > 10000 or obs_mslp[i] > 10000:
					continue
				elif np.isnan(fcst_mslp[i]) or np.isnan(obs_mslp[i]):
					continue
						
				else:
				
					if np.isnan(MJO_phase[i]):
						cc='silver'
						
					else:
					
						cc=colours[int(MJO_phase[i]-1)]
					a = alphas[lt]
				
					plt.scatter(obs_mslp[i], fcst_mslp[i],s=10,c=cc,alpha=a)
						
					x = np.append(x, obs_mslp[i])
					y = np.append(y, fcst_mslp[i])
						
						
		corr = np.corrcoef(x,y)[0,1]
		print corr

		plt.xlim(890,1020)
		plt.xlabel('IBTrACS MSLP (hPa)', fontsize=12)

		plt.ylabel('UKMO '+fcst_type_label+' MSLP (hPa)', fontsize=12)
		plt.ylim(890,1020)
		
		plt.text(0.02,0.02, 'r = {0:.3f}'.format(corr), fontsize=12,transform = ax.transAxes)
	
		legend = ax.legend(handles=[p1,p2,p3,p4,p5,p6,p7,p8,p9], fontsize=9, loc='lower right', title = 'UKMO '+fcst_type_label)
		legend._legend_box.align = "center"
	
	
		plt.savefig("scatter.IBTrACS_vs_UKMO_"+fcst_type+".MSLP.coloured_by_leadtime_and_MJO.2010-2020.png", dpi=400,bbox_inches='tight')
			
		print "scatter.IBTrACS_vs_UKMO_"+fcst_type+".MSLP.coloured_by_leadtime_and_MJO.2010-2020.png"
		
	elif error_type == 'wind':
	
		fig, ax = plt.subplots()
		fig.set_size_inches(6,6)
		
		x = []
		y = []
		
		for lt in [7,6,5,4,3,2,1,0]:	
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

		
			for i in range(len(fcst_wind)):
		
				if fcst_wind[i] > 10000 or obs_wind[i] > 10000:
					continue
				elif np.isnan(fcst_wind[i]) or np.isnan(obs_wind[i]):
					continue
						
				else:
				
					if np.isnan(MJO_phase[i]):
						cc='silver'
						
					else:
					
						cc=colours[int(MJO_phase[i]-1)]
					a = alphas[lt]
				
					plt.scatter(obs_wind[i], fcst_wind[i],s=10,c=cc,alpha=a)
						
					x = np.append(x, obs_wind[i])
					y = np.append(y, fcst_wind[i])
						
					
		corr = np.corrcoef(x,y)[0,1]
		print corr

		plt.xlim(0,300)
		plt.xlabel('IBTrACS Wind Speed (km/h)', fontsize=12)

		plt.ylabel('UKMO '+fcst_type_label+' Wind Speed (km/h)', fontsize=12)
		plt.ylim(0,300)
		
		plt.text(0.8,0.96, 'r = {0:.3f}'.format(corr), fontsize=12,transform = ax.transAxes)
	
		legend = ax.legend(handles=[p1,p2,p3,p4,p5,p6,p7,p8,p9], fontsize=10, loc='upper left', title = 'UKMO '+fcst_type_label)
		legend._legend_box.align = "center"
	
	
		plt.savefig("scatter.IBTrACS_vs_UKMO_"+fcst_type+".wind_speed.coloured_by_leadtime_and_MJO.2010-2020.png", dpi=400,bbox_inches='tight')
			
		print "scatter.IBTrACS_vs_UKMO_"+fcst_type+".wind_speed.coloured_by_leadtime_and_MJO.2010-2020.png"

