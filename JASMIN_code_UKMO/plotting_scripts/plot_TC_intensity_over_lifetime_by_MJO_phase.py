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


def plot_timeseries_intensity(ib_track_file_list,an_track_file_list, outfile):

	fig, ax = plt.subplots()
	fig.set_size_inches(10,6)
	
	xlen=[]



	no_tracks=len(an_track_file_list)
	#colors = cm.viridis(np.linspace(0, 1, no_tracks))
	for file_in_list, i in zip(an_track_file_list, range(len(an_track_file_list))):

		# load in data from a specified TRACK file (reformatted & interpolated)
		#an_data = np.genfromtxt(file_in_list, dtype=float, skip_header=1)
		ib_data = np.genfromtxt(ib_track_file_list[i], dtype=float, skip_header=1)
		
		obs_wind=ib_data[:,10]
		
		
		obs_dates=pl.get_dates(ib_data) 
		
		
		first_date = datetime.datetime.strptime(str(obs_dates[0]), "%Y%m%d%H")
		
		if first_date.strftime("%m-%d") == '02-29':
			continue
			
		else:
			z = MJOdates.index(first_date.strftime("%Y-%m-%d"))
			
			amp = MJOamp[z]
		
			#if amp < 1.0:
				#continue
			
			if MJOphase[z] == MJO and MJOamp[z] >= 1.0:
		
				for i in range(len(obs_wind)):
					if obs_wind[i] > 10000:
						obs_wind[i] = np.nan
					
				print obs_wind
				if all(np.isnan(v) for v in obs_wind):
					print "all nan"
			
					continue
				
				else:
		
					max_wind = np.nanmax(obs_wind)*3.6
			
					obs_wind=obs_wind*3.6
			
					RI = np.zeros(len(obs_wind))
					RIwind = np.zeros(len(obs_wind))
			
					for j in range(len(obs_wind)-4):
				
				
				
						if (obs_wind[j+4] - obs_wind[j]) >= 55.56:
						
							RIwind[j:j+5] = obs_wind[j:j+5]
						if (obs_wind[j+3] - obs_wind[j]) >= 55.56:
							#RI[j:j+3] = obs_wind[j+3] - obs_wind[j]
							RIwind[j:j+4] = obs_wind[j:j+4]	
						if (obs_wind[j+2] - obs_wind[j]) >= 55.56:
							#RI[j:j+2] = obs_wind[j+2] - obs_wind[j]
							RIwind[j:j+3] = obs_wind[j:j+3]	
						if (obs_wind[j+1] - obs_wind[j]) >= 55.56:
							#RI[j:j+1]=obs_wind[j+1] - obs_wind[j]
							RIwind[j:j+2] = obs_wind[j:j+2]
				
					
					#RIwind = np.zeros(len(obs_wind))
			
					for i in reversed(range(len(RIwind))):
				
					
						if RIwind[i] < RIwind[i-1]:
							RIwind[i+1] = 0.0
			
			
					for i in range(len(RIwind)):
			
						if RIwind[i] == 0.0:
							#RI[i] = np.nan
					
							RIwind[i] = np.nan
					
						#elif RI[i] > 0.0:
							#RIwind[i] = obs_wind[i]
					
					#print "RI: ", RI
					print "RIwind: ", RIwind
			
			
			
			
					
			
			
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
				
			
					x = np.arange(0,len(obs_wind),1)
			
					xlen.append(len(x))
			
					print "x: ", x
					print "obs_wind: ", obs_wind
			
					plt.plot(x,obs_wind,color=c,linestyle=ls)
			
					#if RI.any() > 0.:
			
		
					#plt.plot(x[np.where(RI>0.)], obs_wind[np.where(RI>0.)]*3.6, color='blue',linestyle=ls)
					plt.plot(x, RIwind, color='blue',linestyle=ls)
			
			else:
				continue
	print "xlen: ", xlen
	plt.xticks(fontsize=12)
	plt.xlim(0,65)
	xticklocs=np.arange(0,65, step=4)
	plt.ylim(0,300)
	
	
	plt.xticks(xticklocs, xticklocs/4) #sets the location of the xticks (one tick per day = every 4 timesteps), and the values (want to give it in days not timesteps, so /4)
	plt.yticks(fontsize=12)
	plt.xlabel('Days', fontsize = 14)
	plt.ylabel('Maximum Sustained Wind Speed (km/h)', fontsize=14)

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
outfile = "observed_TCs_SIO_"+str(year1[0])+"_"+str(year2[-1])+"_MJO_phase_"+str(MJO)+"_maxwind_timeseries_v2.png"
plot_timeseries_intensity(ib_track_file_list, an_track_file_list, outfile)

