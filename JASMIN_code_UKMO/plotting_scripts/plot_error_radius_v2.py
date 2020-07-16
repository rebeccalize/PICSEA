
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
    
    
   
    
def map_error_radius(radius,outfile, minradius=None,maxradius=None):
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

	
	lat1 = 5
	lat2 = -22
	lon1=35
	lon2=65
	
	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='f')

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
	
	
	centrelat = -8.0
	centrelon = 52.0
	radiusindeg = radius/110.574  #the average location error 3 days ahead in 2015-2016
	
	x,y=m(centrelon,centrelat)
	x2,y2 = m(centrelon,centrelat+radiusindeg) 
	circle1 = plt.Circle((x, y), y2-y, edgecolor='blue',facecolor='blue',alpha=0.5)
	ax.add_patch(circle1)
	
	if minradius:
		minradiusindeg = minradius/110.574
		x,y=m(centrelon,centrelat)
		x2,y2 = m(centrelon, centrelat+minradiusindeg)
		circle2 = plt.Circle((x,y),y2-y,edgecolor='blue',facecolor='None', linestyle='--', linewidth=0.5,alpha=0.5)
		ax.add_patch(circle2)
		
	if maxradius:
		maxradiusindeg = maxradius/110.574
		x,y=m(centrelon,centrelat)
		x2,y2 = m(centrelon, centrelat+maxradiusindeg)
		circle3 = plt.Circle((x,y),y2-y,edgecolor='blue',facecolor='blue', linestyle='--', linewidth=0.5, alpha=0.2)
		ax.add_patch(circle3)		
	

	x,y = m(centrelon,centrelat)
	m.scatter(x,y ,marker=hurricane,edgecolor='k',facecolor='None',s=50,zorder=10,linewidth=0.5)
	



	#legend
	#if cat == "TS":
		#title = cat+" "+str(name)+"\n"
	#else:
		#title="Category "+str(cat)+" "+str(name)+"\n"+str(startdate)+" - "+str(enddate)
	#ib = plt.Line2D((0, 1), (0, 0), color='white', linestyle='--',linewidth=0.5)
	#an = plt.Line2D((0, 1), (0, 0), color='navy', linestyle='-',linewidth=0.5)
	#nwp = plt.Line2D((0, 1), (0, 0), color='orange',linewidth=0.5)
	#legend = ax.legend((an,ib), ['Analysis Track',' '],title=title, fontsize=5, loc='lower left')
	#plt.setp(legend.get_title(), fontsize='5')
	#legend._legend_box.align = "left"

	#save and close the plot
	fig.subplots_adjust(wspace=0)
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.05, dpi=500)
	plt.close()
	

#outfile = "idai_kenneth_analysis_tracks.png"
#map_nwp_tracks_per_storm(idai_an_file,ken_an_file, outfile, "SIO", "3", "Idai","2019022818","2019032400")

#map_error_radius(325,"cyclone_error_radius_n320_nearseychelles.png", minradius=125,maxradius=750)
#map_error_radius(250,"cyclone_error_radius_n768_nearseychelles.png", minradius=100,maxradius=500)

map_error_radius(250,"cyclone_error_radius_n768_nearseychelles.png", minradius=100,maxradius=500)


#293.5 is the average 3-day lead time error in 2015-2016






