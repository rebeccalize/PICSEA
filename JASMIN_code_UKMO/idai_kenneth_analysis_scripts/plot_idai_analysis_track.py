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

fcstdatadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_idai_forecast_tracks/"
ibdatadir = "/gws/nopw/j04/klingaman/emerton/ibtracs_reformatted_data_and_track_maps/"
andatadir = "/gws/nopw/j04/klingaman/emerton/TRACK/ANALYSES/OP_JUL-APR20182019_VOR_VERTAVG_T63/"

#ibtracs_file = ibdatadir+"idai_ibtracs_reformatted.txt"
idai_an_file = andatadir+"idai_analysis_track_reformatted.txt"
ken_an_file = andatadir+"kenneth_analysis_track_reformatted.txt"



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
    
    
def map_nwp_tracks_per_storm(idai_an_file, ken_an_file, outfile, region,cat,name,startdate,enddate):
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

	
	lat1 = 0
	lat2 = -30
	lon1=20
	lon2=60
	
	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')

	RSMClats = [-40, 0, 0, -40]
	RSMClons = [30, 30, 90, 90]

	#draw_rectangle(RSMClats,RSMClons,m)

	#m.bluemarble(alpha=0.8)
	m.drawcoastlines(linewidth=0.6, color='k')
	m.drawcountries(linewidth=0.6, color='k')
	#m.fillcontinents(color='lightgray', alpha=0.4)
	
	hurricane=get_hurricane_symbol()
	
	idai_an_data = np.genfromtxt(idai_an_file, dtype=float, skip_header=1)

	# count how many track points this storm has
	no_points = len(idai_an_data[:,0])
	# set up an array to hold the lat and lon of each point of this storm
	track_coords = np.zeros((no_points, 2))
	# get the lat and lons from the track file for this storm only
	track_coords[:, 0] = idai_an_data[:, 7]  # longitudes
	track_coords[:, 1] = idai_an_data[:, 8]  # latitudes

	# To have the line changing colour from start to finish, split the coordinates into segments
	points = np.array([track_coords[:, 0], track_coords[:, 1]]).T.reshape(-1, 1, 2)
	segments = np.concatenate([points[:-1], points[1:]], axis=1)
	
	#colors = cm.OrRd(np.linspace(0.2, 1, no_points))
	colors = cm.PuBu(np.linspace(0.25, 1, no_points))

		#iterate over each segment of this track
	for p, cl in zip(range(no_points - 1), colors):
		if p < 61:
			z=10
		elif p > 61:
			z=9
		## need to separate the x and y values for each segment, and plot each segment on the map
		xarr = []
		yarr = []
		xarr.append(segments[p][0][0])
		xarr.append(segments[p][1][0])
		yarr.append(segments[p][0][1])
		yarr.append(segments[p][1][1])
		x, y = m(xarr, yarr)
		#m.plot(x, y, linewidth=1.5, color=cl,zorder=z)
		
		
	pp = range(0,no_points,4)
	colors2=cm.PuBu(np.linspace(0,1,24))
	print pp
	for p, cl in zip(pp,colors2):
		if 37 < p < 61:
			xx,yy = m(idai_an_data[p,7], idai_an_data[p,8])
			#m.scatter(xx,yy,marker=hurricane,edgecolors='steelblue',facecolors='None',s=200,linewidth=1.5,zorder=11)
	
		
	#xs,ys=m(idai_an_data[0,7],idai_an_data[0,8])
	#m.scatter(xs,ys,marker=hurricane,edgecolors='orange',facecolors='None',s=50,linewidth=1)

	ken_an_data = np.genfromtxt(ken_an_file,dtype=float,skip_header=1)

	no_points = len(ken_an_data[:, 0])
	# set up an array to hold the lat and lon of each point of this storm
	track_coords = np.zeros((no_points, 2))
	# get the lat and lons from the track file for this storm only
	track_coords[:, 0] = ken_an_data[:, 7]  # longitudes
	track_coords[:, 1] = ken_an_data[:, 8]  # latitudes

	# To have the line changing colour from start to finish, split the coordinates into segments
	points = np.array([track_coords[:, 0], track_coords[:, 1]]).T.reshape(-1, 1, 2)
	segments = np.concatenate([points[:-1], points[1:]], axis=1)

	colors = cm.PuBu(np.linspace(0.2, 1, no_points))

	# iterate over each segment of this track
	for p, cl in zip(range(no_points - 1), colors):
		## need to separate the x and y values for each segment, and plot each segment on the map
		xarr = []
		yarr = []
		xarr.append(segments[p][0][0])
		xarr.append(segments[p][1][0])
		yarr.append(segments[p][0][1])
		yarr.append(segments[p][1][1])
		x, y = m(xarr, yarr)
		m.plot(x, y, linewidth=1.5, color=cl)
		
		if 12 < p < 20:
			xx,yy = m(ken_an_data[p,7], ken_an_data[p,8])
			m.scatter(xx,yy,marker=hurricane,edgecolors='steelblue',facecolors='None',s=200,linewidth=1.5,zorder=11)
			
	
	
	#xs, ys = m(ken_an_data[0, 7], ken_an_data[0, 8])
	#m.scatter(xs, ys, marker=hurricane, edgecolors='skyblue', facecolors='None', s=50, linewidth=1)



	#save and close the plot
	fig.subplots_adjust(wspace=0)
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.05, dpi=500)
	plt.close()
	

outfile = "kenneth_analysis_tracks.png"
map_nwp_tracks_per_storm(idai_an_file,ken_an_file, outfile, "SIO", "3", "Idai","2019022818","2019032400")








