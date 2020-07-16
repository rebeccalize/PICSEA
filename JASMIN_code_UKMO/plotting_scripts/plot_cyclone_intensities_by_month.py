import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.cm as cm
import os
import fnmatch
import matplotlib.lines as mlines



year1=[2006,2007,2008,2009,2010,2011,2012,2013,2014,2015]
year2=[2007,2008,2009,2010,2011,2012,2013,2014,2015,2016]

wind_dict = {}
wind_dict["moz_channel"] = {}
wind_dict["east_of_mada"] = {}

for m in [1,2,3,4,5,6,7,8,9,10,11,12]:
	wind_dict["moz_channel"][m] = []
	wind_dict["east_of_mada"][m] = []
	
	
mslp_dict = {}
mslp_dict["moz_channel"] = {}
mslp_dict["east_of_mada"] = {}

for m in [1,2,3,4,5,6,7,8,9,10,11,12]:
	mslp_dict["moz_channel"][m] = []
	mslp_dict["east_of_mada"][m] = []
	


for y1,y2 in zip(year1,year2):

	print y1, y2
	
	#MOZAMBIQUE CHANNEL

	datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/moz_channel/"

	season_dirs=[]
	for root,dirs,files in os.walk(datadir):
		for dir in dirs:
			season_dirs.append(dir)
	NS = len(season_dirs)
	
	for dir, x in zip(season_dirs,range(NS)):

		list_of_all_files = os.listdir(datadir+dir)
		
		pattern="analysis_*.txt"
		
		for entry in list_of_all_files:
			if fnmatch.fnmatch(entry,pattern):
				ff = datadir+dir+"/"+entry
	
	
		data=np.genfromtxt(ff, dtype=float, skip_header=1)
		mslp=data[:,9]
		wind=data[:,10]
		month=data[:,4]
		
		for i in range(len(wind)):
			if wind[i] > 1000:
				wind[i] = np.nan
				
			if mslp[i] > 10000:
				mslp[i] = np.nan
		
		max_wind = np.nanmax(wind)
		min_mslp = np.nanmin(mslp)
	
		if np.isnan(wind).all():
			continue
		else:
			m = month[np.where(wind == max_wind)][0]

			wind_dict["moz_channel"][m].append(max_wind)
			mslp_dict["moz_channel"][m].append(min_mslp)

			
		
	#EAST OF MADAGASCAR
		
		
	datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/east_of_mada_only/"


	season_dirs=[]
	for root,dirs,files in os.walk(datadir):
		for dir in dirs:
			season_dirs.append(dir)
	NS = len(season_dirs)
	
	for dir, x in zip(season_dirs,range(NS)):
	
		list_of_all_files = os.listdir(datadir+dir)
		
		pattern="analysis_*.txt"
		
		for entry in list_of_all_files:
			if fnmatch.fnmatch(entry,pattern):
				ff = datadir+dir+"/"+entry
	
	
		data=np.genfromtxt(ff, dtype=float, skip_header=1)
		mslp=data[:,9]
		wind=data[:,10]
		month=data[:,4]
		
		for i in range(len(wind)):
			if wind[i] > 1000:
				wind[i] = np.nan
			if mslp[i] > 10000:
				mslp[i] = np.nan
		
		max_wind = np.nanmax(wind)
		min_mslp = np.nanmin(mslp)
		
		if np.isnan(wind).all():
			continue
		else:
			m = month[np.where(wind == max_wind)][0]
					
			wind_dict["east_of_mada"][m].append(max_wind)
			mslp_dict["east_of_mada"][m].append(min_mslp)
		
	
print wind_dict
	
	
fig, ax = plt.subplots()
fig.set_size_inches(6,4)

colors = cm.YlGnBu(np.linspace(0.1, 1, 7))

def get_colour(wind):
	
	w = wind*3.6 #convert m/s to km/h
	
	
	if 0. <= w <= 63.:
		c = colors[0]
	elif 63. <= w < 119.:
		c = colors[1]
	elif 119. <= w < 154.:
		c = colors[2]
	elif 154. <= w < 178.:
		c = colors[3]
	elif 178. <= w <  209.:
		c = colors[4]
	elif 209. <= w < 252.:
		c = colors[5]
	elif w >= 252.:
		c = colors[6]
	else:
		print "something wrong with the wind value?", wind*3.6
	
	return c
	


x = [7,8,9,10,11,12,1,2,3,4,5,6]			
			
for m in [1,2,3,4,5,6,7,8,9,10,11,12]:

	xi = x.index(m)+1

	winds = wind_dict["moz_channel"][m]
	print winds
	
	for i in range(len(winds)):
		if np.isnan(winds[i]):
			continue
		else:
			c = get_colour(winds[i])
			plt.plot(xi, (winds[i]*3.6), marker='o', markerfacecolor='w', markeredgecolor=c)
			
	winds = wind_dict["east_of_mada"][m]
	print winds
	
	for i in range(len(winds)):
		if np.isnan(winds[i]):
			continue
		else:
			c = get_colour(winds[i])
			plt.plot(xi, (winds[i]*3.6), marker='o', markerfacecolor=c, markeredgecolor=c)
			
			
#plt.xlabel('Country', fontsize = 9)
plt.xlim(0.5,12.5)
plt.xticks(np.arange(1,13,1), ('J', 'A', 'S', 'O', 'N', 'D', 'J','F','M','A','M','J'), fontsize=8)
plt.xlabel('Month', fontsize=9)


#plt.xticks(np.arange(1,13,1), fontsize = 8)
plt.ylabel('Maximum Wind Speed (km/h)', fontsize=9)
plt.ylim(50,300)
plt.yticks(np.arange(50,301,25), fontsize=8)


ms = 5
TD = mlines.Line2D([], [], color=colors[0], marker='o', linestyle='None',markersize=ms, label='TD')
TS = mlines.Line2D([], [], color=colors[1], marker='o', linestyle='None',markersize=ms, label='TS')
C1 = mlines.Line2D([], [], color=colors[2], marker='o', linestyle='None',markersize=ms, label='1')
C2 = mlines.Line2D([], [], color=colors[3], marker='o', linestyle='None',markersize=ms, label='2')
C3 = mlines.Line2D([], [], color=colors[4], marker='o', linestyle='None',markersize=ms, label='3')
C4 = mlines.Line2D([], [], color=colors[5], marker='o', linestyle='None',markersize=ms, label='4')
C5 = mlines.Line2D([], [], color=colors[6], marker='o', linestyle='None',markersize=ms, label='5')
mc = mlines.Line2D([], [], marker='o', markerfacecolor='w', markeredgecolor='lightgrey', linestyle='None',markersize=ms, label='Moz. Channel')
eom = mlines.Line2D([], [], marker='o', markerfacecolor='lightgrey', markeredgecolor='lightgrey', linestyle='None',markersize=ms, label='E. of Mada.')

plt.legend(handles=[mc,eom,C5,C4,C3,C2,C1,TS,TD], bbox_to_anchor=(1.34,1.0), title = 'SSHS', fontsize = 8)


#plt.legend((b1[0], b2[0], b3[0], b4[0], b5[0], b6[0], b7[0]), ('TD', 'TS', '1', '2', '3', '4', '5'), bbox_to_anchor=(1.01, 1.0))
#plt.legend((b7[0],b6[0], b5[0], b4[0], b3[0], b2[0], b1[0]), ('5','4','3','2','1','TS','TD'), bbox_to_anchor=(1.0, 1.0))

plt.tight_layout()
	
#plt.show()

plt.savefig("max_wind_speed_ANALYSIS_all_SIO_cyclones_per_month_and_region.png", dpi=400)

plt.close()



fig, ax = plt.subplots()
fig.set_size_inches(6,4)

colors = cm.YlGnBu(np.linspace(0.1, 1, 7))

def get_colour(wind):
	
	w = wind*3.6 #convert m/s to km/h
	
	
	if 0. <= w <= 63.:
		c = colors[0]
	elif 63. <= w < 119.:
		c = colors[1]
	elif 119. <= w < 154.:
		c = colors[2]
	elif 154. <= w < 178.:
		c = colors[3]
	elif 178. <= w <  209.:
		c = colors[4]
	elif 209. <= w < 252.:
		c = colors[5]
	elif w >= 252.:
		c = colors[6]
	else:
		print "something wrong with the wind value?", wind*3.6
	
	return c
	


x = [7,8,9,10,11,12,1,2,3,4,5,6]			
			
for m in [1,2,3,4,5,6,7,8,9,10,11,12]:

	xi = x.index(m)+1

	winds = wind_dict["moz_channel"][m]
	mslp = mslp_dict["moz_channel"][m]
	print winds
	print mslp
	
	for i in range(len(winds)):
		if np.isnan(winds[i]):
			continue
		elif np.isnan(mslp[i]):
			continue
		else:
			c = get_colour(winds[i])
			plt.plot(xi, (mslp[i]), marker='o', markerfacecolor='w', markeredgecolor=c)
			
	winds = wind_dict["east_of_mada"][m]
	mslp = mslp_dict["east_of_mada"][m]
	print winds
	print mslp
	
	for i in range(len(winds)):
		if np.isnan(winds[i]):
			continue
		elif np.isnan(mslp[i]):
			continue
		else:
			c = get_colour(winds[i])
			plt.plot(xi, (mslp[i]), marker='o', markerfacecolor=c, markeredgecolor=c)
			
			
plt.xlim(0.5,12.5)
plt.xticks(np.arange(1,13,1), ('J', 'A', 'S', 'O', 'N', 'D', 'J','F','M','A','M','J'), fontsize=8)
plt.xlabel('Month', fontsize=9)

plt.ylabel('Minimum MSLP (hPa)', fontsize=9)
plt.ylim(1020,900)
plt.yticks(fontsize=8)


ms = 5
TD = mlines.Line2D([], [], color=colors[0], marker='o', linestyle='None',markersize=ms, label='TD')
TS = mlines.Line2D([], [], color=colors[1], marker='o', linestyle='None',markersize=ms, label='TS')
C1 = mlines.Line2D([], [], color=colors[2], marker='o', linestyle='None',markersize=ms, label='1')
C2 = mlines.Line2D([], [], color=colors[3], marker='o', linestyle='None',markersize=ms, label='2')
C3 = mlines.Line2D([], [], color=colors[4], marker='o', linestyle='None',markersize=ms, label='3')
C4 = mlines.Line2D([], [], color=colors[5], marker='o', linestyle='None',markersize=ms, label='4')
C5 = mlines.Line2D([], [], color=colors[6], marker='o', linestyle='None',markersize=ms, label='5')
mc = mlines.Line2D([], [], marker='o', markerfacecolor='w', markeredgecolor='lightgrey', linestyle='None',markersize=ms, label='Moz. Channel')
eom = mlines.Line2D([], [], marker='o', markerfacecolor='lightgrey', markeredgecolor='lightgrey', linestyle='None',markersize=ms, label='E. of Mada.')

plt.legend(handles=[mc,eom,C5,C4,C3,C2,C1,TS,TD], bbox_to_anchor=(1.01,1.0), title = 'SSHS', fontsize = 8)


#plt.legend((b1[0], b2[0], b3[0], b4[0], b5[0], b6[0], b7[0]), ('TD', 'TS', '1', '2', '3', '4', '5'), bbox_to_anchor=(1.01, 1.0))
#plt.legend((b7[0],b6[0], b5[0], b4[0], b3[0], b2[0], b1[0]), ('5','4','3','2','1','TS','TD'), bbox_to_anchor=(1.0, 1.0))

plt.tight_layout()
	
#plt.show()

plt.savefig("min_pressure_ANALYSIS_all_SIO_cyclones_per_month_and_region.png", dpi=400)

	

		

	
	

