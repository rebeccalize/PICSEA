import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.cm as cm
import os
import fnmatch
import matplotlib.lines as mlines



year1=[2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017]
year2=[2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018]

wind_dict = {}
wind_dict["each_season"] = [] #this holds 10 arrays (the number of seasons), each of which contains the max wind speed of each storm in that season
mslp_dict = {}
mslp_dict["each_season"] = []


for y1,y2 in zip(year1,year2):

	datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/"


	season_dirs=[]
	for root,dirs,files in os.walk(datadir):
		for dir in dirs:
			season_dirs.append(dir)
	NS = len(season_dirs)
	
	max_wind_speeds = [None] * NS
	min_mslps = [None] * NS
	
	for dir, x in zip(season_dirs,range(NS)):
	
		#print datadir + dir
		list_of_all_files = os.listdir(datadir+dir)
		
		pattern="ibtracs_*.txt"
		
		for entry in list_of_all_files:
			if fnmatch.fnmatch(entry,pattern):
				analysis_file = datadir+dir+"/"+entry
	
	
		data=np.genfromtxt(analysis_file, dtype=float, skip_header=1)
		mslp=data[:,9]
		wind=data[:,10]
		
		for i in range(len(wind)):
			if wind[i] > 1000:
				wind[i] = np.nan
			if mslp[i] > 10000:
				mslp[i] = np.nan
		
		max_wind_speeds[x] = np.nanmax(wind)
		#print mslp
		print np.nanmin(mslp)
		min_mslps[x] = np.nanmin(mslp)
		
	print min_mslps	
	#print max_wind_speeds
	#print np.nanmax(max_wind_speeds)*3.6
	#wind_arr[np.where(year1==y1)] = max_wind_speeds
	wind_dict["each_season"].append(max_wind_speeds)
	mslp_dict["each_season"].append(min_mslps)
	
print wind_dict
print len(wind_dict["each_season"])

	
	
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



for x in range(len(wind_dict["each_season"])):
	winds = wind_dict["each_season"][x]
	
	for i in range(len(winds)):
		if np.isnan(winds[i]):
			continue
		else:
			c = get_colour(winds[i])
			plt.scatter(x,(winds[i]*3.6), color=c)
			
			
#plt.xlabel('Country', fontsize = 9)
plt.xticks(np.arange(0,12,1), ('2006-2007', '2007-2008', '2008-2009', '2009-2010', '2010-2011', '2011-2012', '2012-2013','2013-2014','2014-2015','2015-2016','2016-2017','2017-2018'), fontsize=8, rotation=45)

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

plt.legend(handles=[C5,C4,C3,C2,C1,TS,TD], bbox_to_anchor=(1.01,1.0), title = 'SSHS')

#plt.legend((b1[0], b2[0], b3[0], b4[0], b5[0], b6[0], b7[0]), ('TD', 'TS', '1', '2', '3', '4', '5'), bbox_to_anchor=(1.01, 1.0))
#plt.legend((b7[0],b6[0], b5[0], b4[0], b3[0], b2[0], b1[0]), ('5','4','3','2','1','TS','TD'), bbox_to_anchor=(1.0, 1.0))

plt.tight_layout()
	
#plt.show()

plt.savefig("max_wind_speed_IBTRACS_all_SIO_cyclones_per_year_2006-2018.png", dpi=400)

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



for x in range(len(mslp_dict["each_season"])):
	winds = wind_dict["each_season"][x]
	mslps = mslp_dict["each_season"][x]
	
	for i in range(len(winds)):
		if np.isnan(winds[i]):
			continue
		else:
			c = get_colour(winds[i])
			plt.scatter(x,(mslps[i]), color=c)
			
			
#plt.xlabel('Country', fontsize = 9)
plt.xticks(np.arange(0,12,1), ('2006-2007', '2007-2008', '2008-2009', '2009-2010', '2010-2011', '2011-2012', '2012-2013','2013-2014','2014-2015','2015-2016','2016-2017','2017-2018'), fontsize=8, rotation=45)

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

plt.legend(handles=[C5,C4,C3,C2,C1,TS,TD], bbox_to_anchor=(1.01,1.0), title = 'SSHS')

#plt.legend((b1[0], b2[0], b3[0], b4[0], b5[0], b6[0], b7[0]), ('TD', 'TS', '1', '2', '3', '4', '5'), bbox_to_anchor=(1.01, 1.0))
#plt.legend((b7[0],b6[0], b5[0], b4[0], b3[0], b2[0], b1[0]), ('5','4','3','2','1','TS','TD'), bbox_to_anchor=(1.0, 1.0))

plt.tight_layout()
	
#plt.show()

plt.savefig("min_mslp_IBTRACS_all_SIO_cyclones_per_year_2006-2018.png", dpi=400)	

		

	
	

