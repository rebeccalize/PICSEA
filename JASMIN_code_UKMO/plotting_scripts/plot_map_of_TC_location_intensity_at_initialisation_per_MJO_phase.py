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


MJO=int(sys.argv[1])


year1s=[2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019] #
year2s=[2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020] # 

tc_lat=[]
tc_lon=[]
tc_wind=[]
tc_mslp=[]



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


for y1, y2 in zip(year1s, year2s):
	
	datadir="/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/" #SIO_storms/"
	 
	 #Find out the number of storms in this TC season and create array of the IDs so we can loop over them
	season_dirs=[]
	for root,dirs,files in os.walk(datadir):
		for dir in dirs:
			season_dirs.append(dir)
			
			NS = len(season_dirs) #total number of storms in this season in the SIO
	
	#Calculate the error of each foreacst at each lead time, for each storm in this TC season / directory, and save
	#Then average across all the forecasts of each storm, and save
	def storm_pos_int(obs_file, nwp_files_list):
		"""This function calculates the statistics for all the forecasts of one storm"""
		
		global tc_lat
		global tc_lon
		global tc_wind
		global tc_mslp
			
		NT = len(nwp_files_list)
		#Get the date, lon, lat and vorticity data for the observed track
		obs_data=np.genfromtxt(obs_file, dtype=float, skip_header=1)
		obs_lon=obs_data[:,7]
		obs_lat=obs_data[:,8]
		obs_mslp=obs_data[:,9]
		obs_wind=obs_data[:,10]
		obs_datelist = pl.get_dates(obs_data)

		#compute the statistics for each forecast of this storm, at each lead time
		"""need to make sure we compare the forecast timestamp with the correct obs timestamp!!"""
		#print NT
		for ff,a in zip(nwp_files_list, range(NT)):
			#print ff
			#get the forecast date, lat, lon and vorticity
			fcst_data=np.genfromtxt(ff,dtype=float,skip_header=1)
			fcst_lon=fcst_data[:,7]
			fcst_lat=fcst_data[:,8]
			fcst_mslp=fcst_data[:,9]
			fcst_wind=fcst_data[:,10]
			fcst_datelist = pl.get_dates(fcst_data)
		
			#check whether first date of this forecast file has an MJO amplitude > 1 (do we want to use it for this MJO phase?)
			#first get the forecast data date into datetime format to compare to the MJO dates
			first_fcst_date = datetime.datetime.strptime(str(fcst_datelist[0]), "%Y%m%d%H")	
		
		
			print first_fcst_date.strftime("%m-%d")
		
			#then find the index of MJO dates at which the date is the first date of the forecast
			if first_fcst_date.strftime("%m-%d") == '02-29': #leap years missing in the MJO data?
		
				continue
			else:		
				z = MJOdates.index(first_fcst_date.strftime("%Y-%m-%d"))
	
				#then find out the MJO amplitude on this date (i.e. at the index of this date)
				amp = MJOamp[z]
				#print "MJO amp", amp
		
				#if the MJO amplitude is less than 1, don't include this foreacst
				#we only want to compute the errors for forecasts started during this MJO phase / phase pair, when the MJO amplitude is >1
				#if amp < 1.0:
					#continue
			
				#if the MJO amplitude was >1 on this day, include this forecast in the stats calculations
				if MJOphase[z] == MJO and MJOamp[z] >= 1.0:			
				
					print fcst_datelist
					print obs_datelist
					
					if fcst_datelist[0] in obs_datelist:
						i = obs_datelist.index(fcst_datelist[0])
					
						print i
						print fcst_datelist[0]
						print obs_datelist[i]
					
						tc_lat.append(obs_lat[i])
						tc_lon.append(obs_lon[i])
						tc_wind.append(obs_wind[i])
						tc_mslp.append(obs_mslp[i])
						
				else:
					continue
						
						
						
	for dir, x in zip(season_dirs,range(len(season_dirs))):
		print datadir+dir
		nwp_files_list = []
		analysis_file=0
		ibtracs_file=0
		#make a list of all the files in this directory
		list_of_all_files = os.listdir(datadir+dir)
		#print list_of_all_files
		pattern="ukmo_nwp*.txt"
		pattern2="analysis_*.txt"
		pattern3="ibtracs_*.txt"
		#for each file in the directory, find the ones that match the forecast filename, analysis filename and ibtracs filename
		for entry in list_of_all_files:
			if fnmatch.fnmatch(entry,pattern):
				nwp_files_list.append(datadir+dir+"/"+entry)
			elif fnmatch.fnmatch(entry,pattern2):
				analysis_file = datadir+dir+"/"+entry
				#print analysis_file
			elif fnmatch.fnmatch(entry,pattern3):
				ibtracs_file = datadir+dir+"/"+entry
		#run the statistics for this storm
		storm_pos_int(ibtracs_file, nwp_files_list)
		
		


for a in range(len(tc_wind)):
	if tc_wind[a] > 10000:
		tc_wind[a] = np.nan
	else:
		tc_wind[a] = tc_wind[a]*3.6
		
	if tc_mslp[a] > 10000:
		tc_mslp[a] = np.nan
	else:
		continue
		


print tc_wind
print tc_mslp

print len(tc_wind)
print len(tc_mslp)
print len(tc_lon)
print len(tc_lat)

print "min tc wind", np.nanmin(tc_wind)
print "max tc wind", np.nanmax(tc_wind)

print "max tc mslp", np.nanmax(tc_mslp)
print "min tc mslp", np.nanmin(tc_mslp)

lat1=10
lat2=-50
lon1=10
lon2=110

fig = plt.figure(figsize=(6, 3))
ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')


#m.fillcontinents(color='white')
m.drawcoastlines(linewidth=0.4, color='k')
m.drawcountries(linewidth=0.4, color='k')
					
cmap = 'YlOrRd'
cmap = plt.get_cmap(cmap)
colors = cmap(np.linspace(0.1, 1, cmap.N))

cmap2 = LinearSegmentedColormap.from_list(cmap, colors)

bounds=np.linspace(50,250,9)
norm = BoundaryNorm(bounds, ncolors=cmap2.N, clip=True)



x,y = m(tc_lon, tc_lat)

im = m.scatter(x, y, c=tc_wind, s= 15, cmap=cmap2, norm=norm)

cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="both", ticks=bounds) #
cbar.ax.set_yticklabels([int(i) for i in bounds])


cbar.ax.tick_params(labelsize=14)  # colorbar font size
plt.text(0.01, 0.02, 'MJO '+str(MJO)+' ('+str(len(tc_lat))+')', transform=ax.transAxes, fontsize=12)
plt.savefig('TC_location_intensity_on_day_forecast_initialised.ibtracs_wind_speed.MJO_'+str(MJO)+'.2010-2019.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
plt.close()

#########################################################################

fig = plt.figure(figsize=(6, 3))
ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')


#m.fillcontinents(color='white')
m.drawcoastlines(linewidth=0.4, color='k')
m.drawcountries(linewidth=0.4, color='k')
					
cmap = 'YlGnBu_r'
cmap = plt.get_cmap(cmap)
colors = cmap(np.linspace(0, 0.9, cmap.N))

cmap2 = LinearSegmentedColormap.from_list(cmap, colors)

bounds=np.linspace(920,1010,10)
norm = BoundaryNorm(bounds, ncolors=cmap2.N, clip=True)



x,y = m(tc_lon, tc_lat)

im = m.scatter(x, y, c=tc_mslp, s= 15, cmap=cmap2, norm=norm)

cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="both", ticks=bounds) #
cbar.ax.set_yticklabels([int(i) for i in bounds])


cbar.ax.tick_params(labelsize=14)  # colorbar font size
plt.text(0.01, 0.02, 'MJO '+str(MJO)+' ('+str(len(tc_lat))+')', transform=ax.transAxes, fontsize=12)
plt.savefig('TC_location_intensity_on_day_forecast_initialised.ibtracs_MSLP.MJO_'+str(MJO)+'.2010-2019.png', bbox_inches='tight', pad_inches=0.5, dpi=400)
plt.close()

