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

ibtracs_file = ibdatadir+"idai_ibtracs_reformatted.txt"
analysis_file = andatadir+"idai_analysis_track_reformatted.txt"

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
    



def map_nwp_tracks_per_storm(nwp_infile_list,ibtracs_file, outfile, region,cat,name,date):
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

	
	lat1 = -12
	lat2 = -24
	lon1=25
	lon2=49
	
	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')

	RSMClats = [-40, 0, 0, -40]
	RSMClons = [30, 30, 90, 90]

	#draw_rectangle(RSMClats,RSMClons,m)

	#m.bluemarble(alpha=0.8)
	m.drawcoastlines(linewidth=0.4, color='gray')
	m.drawcountries(linewidth=0.4, color='gray')
	m.fillcontinents(color='lightgray',alpha=0.3)
	
	hurricane=get_hurricane_symbol()


	no_tracks=len(nwp_infile_list)
	colors = cm.BuPu(np.linspace(0.3, 1, no_tracks)) #Wistia

	#loop over the files containing the forecast tracks, and over the colours to plot each track a different colour
	for file_in_list,c in zip(nwp_infile_list,colors):

		# load in data from a specified TRACK file (reformatted & interpolated)
		data = np.genfromtxt(file_in_list, dtype=float, skip_header=1)

		#plot this forecast track - loops over colours and plots each track a different colour
		x,y = m(data[:,7], data[:,8])
		m.plot(x,y,linewidth=0.5,color=c)
		xs,ys=m(data[0,7],data[0,8])
		m.scatter(xs,ys,marker=hurricane,edgecolors=c,facecolors='None',s=25,linewidth=0.5)

		#OR the following plots each track so it changes colour along its length

		# count how many track points this storm has
		#no_points = len(data[:,0])
		# set up an array to hold the lat and lon of each point of this storm
		#track_coords = np.zeros((no_points, 2))
		# get the lat and lons from the track file for this storm only
		#track_coords[:, 0] = data[:, 7]  # longitudes
		#track_coords[:, 1] = data[:, 8]  # latitudes

		# To have the line changing colour from start to finish, split the coordinates into segments
		#points = np.array([track_coords[:, 0], track_coords[:, 1]]).T.reshape(-1, 1, 2)
		#segments = np.concatenate([points[:-1], points[1:]], axis=1)

		#colors = cm.BuPu(np.linspace(0.2, 1, no_points))

		# iterate over each segment of this track
		#for p, cl in zip(range(no_points - 1), colors):
			## need to separate the x and y values for each segment, and plot each segment on the map
			#xarr = []
			#yarr = []
			#xarr.append(segments[p][0][0])
			#xarr.append(segments[p][1][0])
			#yarr.append(segments[p][0][1])
			#yarr.append(segments[p][1][1])
			#x, y = m(xarr, yarr)
			#m.plot(x, y, linewidth=0.5, color=cl)


	#plot the analysis and ibtracs tracks


	ibtracs_data = np.genfromtxt(ibtracs_file, dtype=float, skip_header=1)
	x,y = m(ibtracs_data[:,7], ibtracs_data[:,8])
	m.plot(x,y, linewidth=0.75, color='k',linestyle='--')
	xs,ys=m(ibtracs_data[0,7],ibtracs_data[0,8])
	m.scatter(xs,ys,marker=hurricane,edgecolors='k',facecolors='None',s=25,linewidth=0.5)
	
	#analysis_data = np.genfromtxt(analysis_file, dtype=float, skip_header=1)
	
	# count how many track points this storm has
	#no_points = len(analysis_data[:,0])
	# set up an array to hold the lat and lon of each point of this storm
	#track_coords = np.zeros((no_points, 2))
	# get the lat and lons from the track file for this storm only
	#track_coords[:, 0] = analysis_data[:, 7]  # longitudes
	#track_coords[:, 1] = analysis_data[:, 8]  # latitudes

	# To have the line changing colour from start to finish, split the coordinates into segments
	#points = np.array([track_coords[:, 0], track_coords[:, 1]]).T.reshape(-1, 1, 2)
	#segments = np.concatenate([points[:-1], points[1:]], axis=1)
	
	#colors = cm.Greys(np.linspace(0.2, 1, no_points))

		#iterate over each segment of this track
	#for p, cl in zip(range(no_points - 1), colors):
		## need to separate the x and y values for each segment, and plot each segment on the map
		#xarr = []
		#yarr = []
		#xarr.append(segments[p][0][0])
		#xarr.append(segments[p][1][0])
		#yarr.append(segments[p][0][1])
		#yarr.append(segments[p][1][1])
		#x, y = m(xarr, yarr)
		#m.plot(x, y, linewidth=0.5, color=cl)
	
	#x,y = m(analysis_data[:,7], analysis_data[:,8])
	#m.plot(x,y,linewidth=0.75, color='k',linestyle='-')
	#xs,ys=m(analysis_data[0,7],analysis_data[0,8])
	#m.scatter(xs,ys,marker=hurricane,edgecolors='k',facecolors='None',s=25,linewidth=0.5)

	#legend
	if cat == "TS":
		title = cat+" "+str(name)+"\n"+date
	else:
		title="Category "+str(cat)+" "+str(name)+"\n"+str(date)
	ib = plt.Line2D((0, 1), (0, 0), color='k', linestyle='--',linewidth=0.5)
	#an = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-',linewidth=0.5)
	nwp = plt.Line2D((0, 1), (0, 0), color='purple',linewidth=0.5)
	legend = ax.legend((ib, nwp), ['Observed Track', 'UKMO Det. Forecasts'],title=title, fontsize=5, loc='lower left')
	plt.setp(legend.get_title(), fontsize='5')
	legend._legend_box.align = "left"

	#save and close the plot
	fig.subplots_adjust(wspace=0)
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.05, dpi=500)
	plt.close()


		

nwp_infile_list = []
dates = [2019030412, 2019030500, 2019030512,2019030600,2019030612,2019030700,2019030800,2019030812,2019030900,2019030912,2019031000,2019031012,2019031100,2019031112,2019031212,2019031300,2019031312,2019031400]
count=0
for date in dates:
	count+=1
	nwp_infile_list.append(fcstdatadir+"idai_UKMO_"+str(date)+"_reformatted.txt")
	outfile="idai_ukmo_nwp_forecast_tracks_"+str(count)+"_v2.png"
	print "mapping no. ", count
	map_nwp_tracks_per_storm(nwp_infile_list,ibtracs_file,outfile,"SIO","3","Idai",date)

outfile = "idai_ukmo_nwp_forecast_tracks_v2.png"	


map_nwp_tracks_per_storm(nwp_infile_list, ibtracs_file, outfile, "SIO", "3", "Idai"," ")
        
        
        
        
        
        
        
        
        
        
    

