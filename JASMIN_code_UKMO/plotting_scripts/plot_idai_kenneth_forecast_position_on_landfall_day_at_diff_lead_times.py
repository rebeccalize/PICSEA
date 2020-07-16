#writing function to plot map of SOI / SE Africa with TC tracks (initially from IBTrACS)

import picsea_library as pl
import fnmatch
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.cm as cm
import numpy as np
import matplotlib.path as mpath
import math
from osgeo import ogr, osr
import matplotlib.colors as colors
from matplotlib.patches import Circle



def get_hurricane_symbol():
    u = np.array([  [2.444,7.553],
                    [0.513,7.046],
                    [-1.243,5.433],
                    [-2.353,2.975],
                    [-2.578,0.092],
                    [-2.075,-1.795],
                    [-0.336,-2.870],
                    [2.609,-2.016]  ])
    u[:,0] -= 0.098
    codes = [1] + [2]*(len(u)-2) + [2] 
    u = np.append(u, -u[::-1], axis=0)
    codes += codes

    return mpath.Path(3*u, codes, closed=False)
    
    
   
    
def map_error_radius(obs_lat,obs_lon,lat_list,lon_list,c,label,outfile):
	"""Plots a map of storm tracks in the Southern Hemisphere
	Plots a map of storm tracks for one storm in the Southern Hemisphere
	from a list of several TRACK files (e.g. one per year) that have been reformatted to suit Python
	Plots all forecasts of one storm, and the reanalysis and ibtracs tracks.
	Plots the tracks so that they change colour from lighter to darker with time
	Note: the colormap used is hardcoded in the function (where 'colors' is assigned)
	Function takes a list of filenames as  input argument"""

	# set up map of region
	#if region == "SH":
		#lat1 = 30
		#lat2 = -60
		#lon1 = -25
		#lon2 = 335

	#elif region == "SIO":
		#lat1=30
		#lat2=-55
		#lon1=-10
		#lon2=130

	
	lat1 = 15
	lat2 = -32
	lon1=25
	lon2=75
	
	lat1 = -2
	lat2 = -33
	lon1=21
	lon2=62
	
	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')

	RSMClats = [-40, 0, 0, -40]
	RSMClons = [30, 30, 90, 90]

	#draw_rectangle(RSMClats,RSMClons,m)

	#m.bluemarble(alpha=0.8)
	m.drawcoastlines(linewidth=0.4, color='darkgray')
	m.drawcountries(linewidth=0.4, color='darkgray')
	m.fillcontinents(color='lightgray', alpha=0.4)
	
	hurricane=get_hurricane_symbol()
	
	#xs, ys = m(42.0, -10.0)
	#m.scatter(xs, ys, marker=hurricane, edgecolors='royalblue', facecolors='None', s=100, linewidth=1)
	
	centrelat=obs_lat
	centrelon=obs_lon
	
	
	for lat,lon,a in zip(lat_list,lon_list,[1.0,0.7,0.5,0.15]): 
	
		#radiusindeg = r/110.574  
	
		x,y=m(lon,lat)
		#x2,y2 = m(centrelon,centrelat+radiusindeg)
		
		#circle = plt.Circle((x, y), y2-y, edgecolor=c,facecolor=c,linewidth = 0.5, alpha=a)
		#ax.add_patch(circle)
		
		m.scatter(x,y ,marker=hurricane,edgecolor=c,facecolor='None',s=60,zorder=10,linewidth=0.5,alpha=a)
	

	x,y = m(centrelon,centrelat)
	m.scatter(x,y ,marker=hurricane,edgecolor='k',facecolor='None',s=60,zorder=10,linewidth=0.5,alpha=0.95)
	

	#save and close the plot
	fig.subplots_adjust(wspace=0)
	
	#plt.text(0.02,0.02, label,fontsize=7,transform = ax.transAxes)
	
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.05, dpi=500)
	plt.close()
	


#map_error_radius([101.3,267.88,499.11,762.9],'b','UKMO HRES\n2010 - 2020\nAverage track error 1,3,5 & 7 days ahead',"track_error_radius_map.2010-2020.average_error.lead_times_1_3_5_7_days.png") #2010-2020 lead times 1,3,5,7
#map_error_radius([84.03,226.1,441.19,719.53],'b','UKMO HRES\n2017 - 2020\nAverage track error 1,3,5 & 7 days ahead',"track_error_radius_map.2017-2020.average_error.lead_times_1_3_5_7_days.png") #2017-2020 lead times 1,3,5,7

#map_error_radius([51.3,99.69,246.28,759.26],'r','UKMO HRES\nCyclone Idai\nTrack error 1,3,5 & 7 days ahead of landfall',"track_error_radius_map.Cyclone_Idai.actual_error_lead_times_1_3_5_7_days_before_landfall.png") #2017-2020 lead times 1,3,5,7

map_error_radius(-12.1,40.9,[-12.4422,-12.2214,-10.3219,-10.7281],[40.4373,43.21,49.5193,52.3347],'r','UKMO HRES\nCyclone Kenneth\nTrack error 1,3,5 & 7 days ahead of landfall',"TC_landfall_actual_position_error.Cyclone_Kenneth.lead_times_1_3_5_7_days_before_landfall.png") #2017-2020 lead times 1,3,5,7

map_error_radius(-19.9,36.3,[-19.5098,-19.3741,-19.3143,-17.4662],[36.0218,35.5539,38.4425,29.8816],'r','UKMO HRES\nCyclone Idai\nTrack error 1,3,5 & 7 days ahead of landfall',"TC_landfall_actual_position_error.Cyclone_Idai.lead_times_1_3_5_7_days_before_landfall.png")



