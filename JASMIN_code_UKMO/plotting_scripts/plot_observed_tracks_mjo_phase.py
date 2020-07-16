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

MJO=int(sys.argv[1])

print MJO
print type(MJO)

#MJO_file = "/gws/nopw/j04/klingaman/datasets/MJO_INDICES/MJO_phase"+str(MJO)+".jan1979-may2020_dmeans_ts.nc"

#This MJO file just has the phase and amplitude on every day, not split up by phase, probably easiest to use for any analysis with individual phases
MJO_file = "/gws/nopw/j04/klingaman/datasets/MJO_INDICES/MJO_indices_rmm1_rmm2.jan1979-may2020_dmeans_ts.index_values.nc"

#MJO_file = "/gws/nopw/j04/klingaman/datasets/MJO_INDICES/MJO_rmm1_rmm2.jan-dec_dmeans_ts.1979-2019.nc"

print MJO_file

ffMJO = Dataset(MJO_file,'r')
MJOamp = ffMJO.variables['amplitude'][:]
MJOphase = ffMJO.variables['phase'][:]

#MJOdatesnc = ffMJO.variables['time'][:]
#t_unit = ffMJO.variables['time'].units
#t_cal = ffMJO.variables['time'].calendar
#tvalue = num2date(MJOdatesnc,units=t_unit, calendar=t_cal)
#MJOdates = [i.strftime("%Y-%m-%d") for i in tvalue]

start_time = datetime.datetime(1979,1,1)
MJO_tvalue = np.array([start_time + datetime.timedelta(days=i) for i in xrange(len(MJOamp))])
			
MJOdates = [j.strftime("%Y-%m-%d") for j in MJO_tvalue]

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


def map_all_obs_tracks_SIO(ib_track_file_list,an_track_file_list, outfile, region):
	"""Plots a map of storm tracks (analysis track) in the SWIO, that had their genesis in a given MJO phase (pair)
	Colour-codes the track according to the SWIO intensity category scale, using max winds from IBTrACS"""

	# set up map of region
	if region == "SH":
		lat1 = 30
		lat2 = -60
		lon1 = -25
		lon2 = 335

	elif region == "SIO":
		lat1=10
		lat2=-50
		lon1=10
		lon2=110

	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')

	def draw_rectangle( lats, lons, m):
		x, y = m( lons, lats )
		xy = zip(x,y)
		poly = Polygon(xy, facecolor='None', edgecolor='darkgrey', linewidth=0.75,alpha=0.75)
		plt.gca().add_patch(poly)

	RSMClats = [-40, 0, 0, -40]
	RSMClons = [30, 30, 90, 90]

	draw_rectangle(RSMClats,RSMClons,m)

	m.fillcontinents(color='white')
	m.drawcoastlines(linewidth=0.4, color='k')
	m.drawcountries(linewidth=0.4, color='k')

	colours = ['#fcc200', '#f05238', '#a1005c', '#08025c']

	no_tracks=len(an_track_file_list)
	#colors = cm.viridis(np.linspace(0, 1, no_tracks))
	for file_in_list, i in zip(ib_track_file_list, range(len(ib_track_file_list))):

		# load in data from a specified TRACK file (reformatted & interpolated)
		#an_data = np.genfromtxt(file_in_list, dtype=float, skip_header=1)
		ib_data = np.genfromtxt(ib_track_file_list[i], dtype=float, skip_header=1)
		
		obs_wind=ib_data[:,10]
		
		obs_dates=pl.get_dates(ib_data) #get the start date of the ibtracs track
		
		first_date = datetime.datetime.strptime(str(obs_dates[0]), "%Y%m%d%H")
		
		#print "TC start date: ", first_date.strftime("%Y-%m-%d")
		
		
		if first_date.strftime("%m-%d") == '02-29':
			continue
			
		else:
		
			#print obs_dates[0]
			#print first_date
			
			z = MJOdates.index(first_date.strftime("%Y-%m-%d"))
			#print "Date in MJO file at this index: ", MJOdates[z]
			
			amp = MJOamp[z]
			
			#print amp
			#print "MJO phase: ", MJOphase[z]
			#print "Amplitude: ", MJOamp[z]
			
			
			
			if MJOphase[z] == MJO and MJOamp[z] >= 1.0:
			
				print MJOdates[z]
				print MJOphase[z]
				print MJOamp[z]
				
				
		
				for i in range(len(obs_wind)):
					if obs_wind[i] > 10000:
						obs_wind[i] = np.nan
		
				max_wind = np.nanmax(obs_wind)*3.6
			
				if max_wind < 51:
					c = 'khaki'
					ls = 'dotted'
				elif 51<= max_wind <63:
					c = 'gold'
					ls = 'dotted'
				elif 63 <= max_wind < 89:
					c = 'darkorange'
					ls = 'dotted'
				elif 89 <= max_wind < 118:
					c = 'black'
					ls = 'dotted'
				elif 118 <= max_wind < 166:
					c = 'orangered'
					ls = '-'
				elif 166 <= max_wind < 213:
					c = 'firebrick'
					ls ='-'
				elif max_wind >= 213:
					c = 'k'
					ls = '-'
				elif np.isnan(max_wind):
					c = 'grey'
					ls ='--'
				

				#plot this forecast track
				x,y = m(ib_data[:,7], ib_data[:,8])
				m.plot(x,y,linewidth=1,color=c, linestyle=ls)
			
				hurricane=get_hurricane_symbol()
				xs,ys=m(ib_data[0,7],ib_data[0,8])
				m.scatter(xs,ys,marker=hurricane,edgecolors=c,facecolors='None',s=50,linewidth=0.6)

			else:
				continue
	#legend
	tdf = plt.Line2D((0, 1), (0, 0), color='khaki',linestyle='dotted',linewidth=1)
	td = plt.Line2D((0, 1), (0, 0), color='yellow', linestyle='dotted',linewidth=1)
	tsm = plt.Line2D((0, 1), (0, 0), color='darkorange', linestyle='dotted',linewidth=1)
	tsf = plt.Line2D((0, 1), (0, 0), color='black',linestyle='dotted',linewidth=1)
	tc = plt.Line2D((0, 1), (0, 0), color='orangered',linewidth=1)
	tci = plt.Line2D((0, 1), (0, 0), color='firebrick',linewidth=1)
	tcti = plt.Line2D((0, 1), (0, 0), color='black',linewidth=1)
	title = "2010-2018\n"+str(no_tracks)+" Cyclones\nCategory:"
	#legend = ax.legend((tdf, td, tsm, tsf, tc, tci, tcti), ['Weak Depression', 'Tropical Depression', 'Moderate Tropical Storm', 'Strong Tropical Storm', 'Tropical Cyclone', 'Intense Tropical Cyclone', 'Very Intense Tropical Cyclone'], fontsize=5, loc='lower left')
	#plt.setp(legend.get_title(), fontsize='5')

	#save and close the plot
	fig.subplots_adjust(wspace=0)
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.05, dpi=500)
	plt.close()

datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"
#datadir="./"

an_track_file_list=[]
ib_track_file_list=[]
category_list=[]

#find all the track files in this TC season (e.g. 2015-2016)
year1=[2010,2011,2012,2013,2014,2015,2016,2017,2018,2019]
year2=[2011,2012,2013,2014,2015,2016,2017,2018,2019,2020]
for y1,y2 in zip(year1,year2):

	print y1,y2

	#cat_file = datadir + str(y1) + "_" + str(y2) + "/SIO_storms/storm_categories.txt"
	#cf = open(cat_file, 'r')
	#lines = cf.readlines()
	#trnos, cats = ([] for i in range(2))
	#for x in lines:
		#print x
		#trnos.append(x.split(' ')[0])
		#cats.append(x.split(' ')[1])

	for root,dirs,files in os.walk(datadir+str(y1)+"_"+str(y2)+"/SIO_storms"):
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
			list_of_all_files = os.listdir(datadir+str(y1)+"_"+str(y2)+"/SIO_storms/"+dir)
			pattern1="ibtracs_*.txt"
			pattern2="analysis_*.txt"
			#print "directory:", dir
			#for each file in the directory, find the ones that match the forecast filename, analysis filename and ibtracs filename
			for entry in list_of_all_files:
				if fnmatch.fnmatch(entry,pattern1):
					ib_track_file_list.append(datadir+str(y1)+"_"+str(y2)+"/SIO_storms/"+dir+"/"+entry)
				elif fnmatch.fnmatch(entry,pattern2):
					an_track_file_list.append(datadir+str(y1)+"_"+str(y2)+"/SIO_storms/"+dir+"/"+entry)

#print track_file_list
#outdir="/group_workspaces/jasmin2/klingaman/emerton/ukmo_nwp_analysis/"
outdir="./"
outfile = "observed_TCs_SIO_"+str(year1[0])+"_"+str(year2[-1])+"_MJO_phase_"+str(MJO)+"_ibtracs.no_legend.png"
map_all_obs_tracks_SIO(ib_track_file_list, an_track_file_list, outfile, "SIO")

