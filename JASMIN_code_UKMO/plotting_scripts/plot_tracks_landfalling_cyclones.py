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
import matplotlib.path as mpath


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


def map_all_obs_tracks_SIO(track_file_list, category_list, outfile, region):
	"""Plots a map of storm tracks in the Southern Hemisphere
	from a list of several TRACK files that have been reformatted to suit Python
	Plots all analysis tracks in the range of years specified, in the SIO
	Function takes a list of filenames as  input argument

	Would like to further include different colours for each storm's category"""

	# set up map of region
	if region == "SH":
		lat1 = 30
		lat2 = -60
		lon1 = -25
		lon2 = 335

	elif region == "SIO":
		lat1=30
		lat2=-55
		lon1=-10
		lon2=130
		
	elif region == "SIO2":
		lat1=0
		lat2=-40
		lon1=20
		lon2=90
		
	hurricane=get_hurricane_symbol()

	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='f')

	def draw_rectangle( lats, lons, m):
		x, y = m( lons, lats )
		xy = zip(x,y)
		poly = Polygon(xy, facecolor='None', edgecolor='darkgrey', linewidth=1,alpha=0.75)
		plt.gca().add_patch(poly)

	RSMClats = [-40, 0, 0, -40]
	RSMClons = [30, 30, 90, 90]

	#draw_rectangle(RSMClats,RSMClons,m)

	m.fillcontinents(color='lightgrey')
	m.drawcoastlines(linewidth=0.3, color='k')
	m.drawcountries(linewidth=0.3, color='k')



	no_tracks=len(track_file_list)
	colors = cm.viridis(np.linspace(0, 1, no_tracks))
	for file_in_list, i in zip(track_file_list, range(len(track_file_list))):

		#print category_list[i]

		if category_list[i] == 0:
			c='#bc5090'
		elif category_list[i] == 1:
			c="#003f5c"
		elif category_list[i] == 2:
			c = "#ffa600"
		#elif category_list[i] == str(3):
			#c="dodgerblue"
		#elif category_list[i] == str(4):
			#c = "darkblue"
		#elif category_list[i] == str(5):
			#c="darkblue"
		else:
			c="darkgrey"

		# load in data from a specified TRACK file (reformatted & interpolated)
		data = np.genfromtxt(file_in_list, dtype=float, skip_header=1)

		#plot this forecast track
		x,y = m(data[:,7], data[:,8])
		m.plot(x,y,linewidth=0.5,color=c)
		
		xs,ys=m(data[0,7],data[0,8])
		m.scatter(xs,ys,marker=hurricane,edgecolors=c,facecolors='None',s=50,linewidth=0.6)

	#legend
	#ts = plt.Line2D((0, 1), (0, 0), color='lightskyblue',linewidth=0.5)
	#c1 = plt.Line2D((0, 1), (0, 0), color='deepskyblue',linewidth=0.5)
	#c2 = plt.Line2D((0, 1), (0, 0), color='dodgerblue', linestyle='--',linewidth=0.5)
	#c3 = plt.Line2D((0, 1), (0, 0), color='dodgerblue',linewidth=0.5)
	#c4 = plt.Line2D((0, 1), (0, 0), color='mediumblue',linewidth=0.5)
	#c5 = plt.Line2D((0, 1), (0, 0), color='darkblue',linewidth=0.5)
	mada = plt.Line2D((0, 1), (0, 0), color='#bc5090',linewidth=0.5)
	moz = plt.Line2D((0, 1), (0, 0), color='#003f5c',linewidth=0.5)
	sey = plt.Line2D((0, 1), (0, 0), color='#ffa600',linewidth=0.5)
	#title = "2006-2018\n"+str(no_tracks)+" Cyclones\nCategory:"
	#legend = ax.legend((ts, c3, c5), ['0-1', '2-3', '4-5'],title=title, fontsize=5, loc='lower left')
	#plt.setp(legend.get_title(), fontsize='5')

	#save and close the plot
	fig.subplots_adjust(wspace=0)
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.05, dpi=500)
	plt.close()

datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"
#datadir="./"

track_file_list=[]
category_list=[]
#find all the track files in this TC season (e.g. 2015-2016)
#year1=[2006, 2007, 2008, 2009, 2010,2011,2012,2013,2014,2015,2016,2017]
#year2=[2007, 2008, 2009, 2010, 2011,2012,2013,2014,2015,2016,2017,2018]

countries = ['madagascar','mozambique','seychelles']

for country in countries:

	#cat_file = datadir + str(y1) + "_" + str(y2) + "/SIO_storms/storm_categories.txt"
	#cf = open(cat_file, 'r')
	#lines = cf.readlines()
	#trnos, cats = ([] for i in range(2))
	#for x in lines:
		#print x
		#trnos.append(x.split(' ')[0])
		#cats.append(x.split(' ')[1])

	for root,dirs,files in os.walk(datadir+str(country)+"_landfalling"):
		#print root
		#print dirs
		#loop over the directories (i.e. loop over the storms)
		for dir in dirs:
		


			#if dir in trnos:
				#si = [trnos.index(i) for i in trnos if dir in i]
				#si = si[0]
				#cat = cats[si]
				#category_list.append(cat)
			#else:
				#category_list.append(np.nan)


			#make a list of all the files in this directory
			list_of_all_files = os.listdir(datadir+str(country)+"_landfalling/"+dir)
			pattern="ibtracs_*.txt"
			#print "directory:", dir
			#for each file in the directory, find the ones that match the forecast filename, analysis filename and ibtracs filename
			for entry in list_of_all_files:
				if fnmatch.fnmatch(entry,pattern):
					track_file_list.append(datadir+str(country)+"_landfalling/"+dir+"/"+entry)
					
					if country == 'madagascar':
						category_list.append(0)
					elif country == 'mozambique':
						category_list.append(1)
					elif country == 'seychelles':
						category_list.append(2)
						
#print track_file_list
#print category_list

#print track_file_list
#outdir="/group_workspaces/jasmin2/klingaman/emerton/ukmo_nwp_analysis/"
outdir="./"
outfile = "all_landfalling_TCs_SIO_2006-2018.png"
map_all_obs_tracks_SIO(track_file_list, category_list, outfile, "SIO2")

