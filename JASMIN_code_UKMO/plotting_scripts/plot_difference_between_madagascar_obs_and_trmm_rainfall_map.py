import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
from calendar import monthrange
from datetime import date
import operator
import math
import matplotlib.cm as cm
import matplotlib as mpl

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


	
#read in observation data
obs_data = np.genfromtxt("/home/users/emerton/MADAGASCAR_RR_daily_2006_2018_CSV.csv", delimiter=',',dtype=str)
	
station_names = obs_data[0,1:]
station_IDs = obs_data[1,1:]
station_lons = obs_data[2,1:]
station_lats = obs_data[3,1:]


obs_precip = obs_data[4:,:] #first value of each row [:,0] is the date! each column [x,:] is a different station

#set up a dictionary to hold the obs and trmm precip at each station location / grid point:
pcp_dict = {}
pcp_dict["obs"] = {}
pcp_dict["trmm"] = {}
pcp_dict["error"] = {}
pcp_dict["perc_error"] = {}

for ID in station_IDs:
	pcp_dict["obs"][ID] = []
	pcp_dict["trmm"][ID] = []
	pcp_dict["error"][ID] = []
	pcp_dict["perc_error"][ID] = []

#use these dictionaries to get the error each day, and then at the end, average over all the days to get the average error at each station, and then map this. 



no_stations = len(station_IDs)

#need to include something to exclude any negative rainfall obs (think missing value uses -99)


trmm_dir="/gws/nopw/j04/klingaman/datasets/TRMM_3B42/V7_NC_daily/"


for year in [2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018]: #,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018
	
	#read in the daily TRMM data for this year
	
	trmm_file = trmm_dir+"3B42_daily."+str(year)+".nc"
	
	fftrmm = Dataset(trmm_file,'r')
	#trmm_pcp = fftrmm.variables['precipitation'][:]
	
	if year < 2016:		
		trmm_pcp_array = fftrmm.variables['r'][:]
		
	elif year >= 2016:
		trmm_pcp_array = fftrmm.variables['precipitation'][:]
		
	trmm_lons = fftrmm.variables['longitude'][:]
	trmm_lats = fftrmm.variables['latitude'][:]
	
	
	#for s in range(1,no_stations+1):
	
		#print obs_precip[:,s]
		#print obs_precip[:,0]
		
		
	d0 = date(year,1,1)
	
	
	for month in [1,2,3,4,5,6,7,8,9,10,11,12]: #,2,3,4,5,6,7,8,9,10,11,12
		
		days=monthrange(year,month)[1]
		
		for day in range(1,days+1): #1,days+1
		
			print year,month,day
			#d = day+1
			
			#get the index of this date in the observation data array
			i = np.where(obs_precip[:,0] == str(year)+str(month).zfill(2)+str(day).zfill(2))
			
			#testing:	
			#print 'looking for ',str(year)+str(month).zfill(2)+str(day).zfill(2)	
			#print 'index: ', i
			#print 'date found: ',obs_precip[i,0]
			
			
			#get index of this date in the TRMM files (no date variable, file starts at 1/1/year, so find how many days this day is from 1st Jan)		
			d1 = date(year,month,day)
			delta = d1-d0
			j=delta.days
			#print j
			
			#for each station, get the observed precip at this index (date):
			for s,ID in zip(range(1,no_stations+1), station_IDs):
			
				#print "station ID: ", ID
			
				station_pcp = obs_precip[i,s][0]
				
				#not sure why this loop is necessary to convert to a float, but it's the only thing I've found that works...				
				for item in station_pcp: 
				
					a = float(item)
					#print "station pcp: ", a
				
					if a >= 0.0:
						pcp_dict["obs"][ID].append(a)
					
					else:
						#print "missing value?: ", a
						pcp_dict["obs"][ID].append(np.ma.masked) #if it doesn't append something, the dates won't match up
						
						
						
				slat = float(station_lats[s-1])
				slon = float(station_lons[s-1])
				
				#get nearest lat and lon values in TRMM to the station lat and lon, and the index of these in the trmm lat and lon arrays
				tlat = min(range(len(trmm_lats)), key=lambda y: abs(trmm_lats[y]-slat))
				tlon = min(range(len(trmm_lons)), key=lambda x: abs(trmm_lons[x]-slon))
				
				#print slat, trmm_lats[tlat]
				#print slon, trmm_lons[tlon]
				
				#get the TRMM precip at this grid point, to match up with the observed precip at the station at that location
				trmm_pcp = trmm_pcp_array[j,tlat,tlon] #[time,lat,lon] is the trmm structure
				
				#print "trmm_pcp: ",trmm_pcp
				pcp_dict["trmm"][ID].append(trmm_pcp)
				
				
				
				
				#########################################################################
				# up to this point, the code has all been about reading in the data, and getting the trmm data at the grid point closest
				# to each rainfall station, on the same dates
				#########################################################################


for ID in station_IDs:
	for i in range(len(pcp_dict["obs"][ID][:])):
	
	
	
		error = pcp_dict["trmm"][ID][i] - pcp_dict["obs"][ID][i]
		#print error
		pcp_dict["error"][ID].append(error)
		
		if pcp_dict["obs"][ID][i] == 0.0: #can't divide by 0, so ignore any days when the observed rainfall is 0. shame, so also look at total bias rather than %
			continue
			
		else:
		
			perc_error = (error / pcp_dict["obs"][ID][i]) * 100
		#print perc_error, "%"

		
		pcp_dict["perc_error"][ID].append(perc_error)
		
		
means_perc=[]
means=[]
maxs=[]
mins=[]
for ID in station_IDs:
	mean_perc = np.nanmean(pcp_dict["perc_error"][ID][:])
	means_perc.append(mean_perc)
	
	mean=np.nanmean(pcp_dict["error"][ID][:])
	means.append(mean)
	
	max=np.nanmax(pcp_dict["error"][ID][:])
	maxs.append(max)
	min=np.nanmin(pcp_dict["error"][ID][:])
	mins.append(min)
	
print means
imax = np.max(means)
imin = abs(np.min(means))	

print "max min:", np.max(mins)

print imax
print imin
					
#just madagascar
lat1=-9 
lat2=-28
lon1=42
lon2=52.5

fig = plt.figure(figsize=(6, 3))
ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])
#ax1 = fig.add_axes([0.8, 0.1, 0.03, 0.8]) 

m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='h')


#m.bluemarble(alpha=0.85)
#m.shadedrelief()
m.drawcoastlines(linewidth=0.4, color='k')
m.drawcountries(linewidth=0.4, color='k')
#m.fillcontinents(color='white')




#normbl = mpl.colors.Normalize(vmin=0, vmax=500)
#cmapbl = cm.Blues
#mbl=cm.ScalarMappable(norm=normbl,cmap=cmapbl)


#normrd = mpl.colors.Normalize(vmin=0, vmax=200)
#cmaprd = cm.Reds
#mrd=cm.ScalarMappable(norm=normrd,cmap=cmaprd)

#plt.colorbar()

axins = inset_axes(ax,
                   width="5%",  # width = 5% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.05, 0., 1, 1),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )


#norm = mpl.colors.Normalize(vmin=0, vmax=200)
norm = mpl.colors.Normalize(vmin=-200, vmax=0)
#cmap = cm.RdBu
#cmap = cm.Blues
cmap = cm.Reds_r
mm = cm.ScalarMappable(norm=norm,cmap=cmap)

#plt.colorbar(cmap=cm.RdBu,vmin=-500,vmax=500)
cb1 = mpl.colorbar.ColorbarBase(axins,cmap=cmap,norm=norm,orientation='vertical', extend='min') #, extend='both'
cb1.set_label('minimum daily rainfall error (mm)', fontsize=6)
cb1.ax.tick_params(labelsize=6)



#xs,ys = m(station_lons,station_lats)
#plt.scatter(xs,ys,10,marker='o',edgecolor='k',facecolor=cmap(norm(means)),linewidth=0.5,zorder=10)


for ID,i in zip(station_IDs, range(no_stations)):

	x, y = m(float(station_lons[i]), float(station_lats[i]))
	z=mins[i]
	
	#print "pcp_dict[perc_error][ID][:] ", pcp_dict["perc_error"][ID][:]
	
	#mean_perc = np.nanmean(pcp_dict["perc_error"][ID][:])
	
	#print "mean_perc: ", mean_perc
	
	if z == 0.:
		
		ax.scatter(x,y,10,marker='o',edgecolor='k',facecolor='white',linewidth=0.5,zorder=10)
		
	else: 
	
		cl=mm.to_rgba(z)
		ax.scatter(x,y,13,marker='o',edgecolor='k',facecolor=cl,linewidth=0.5,zorder=10)
	
	#if z > 0.0:
	
		#cl=mbl.to_rgba(z)
	
		#plt.scatter(x,y,10,marker='o',edgecolor='k',facecolor=cl,linewidth=0.5,zorder=10)
	
	#elif z <= 0.0:
	
		#cl=mrd.to_rgba(z)
		
		#plt.scatter(x,y,10,marker='o',edgecolor='k',facecolor=cl,linewidth=0.5,zorder=10)
		#continue
		

	
	#if station_names[i] == "ANTANANARIVO" or station_names[i] == "MAHANORO" or station_names[i] == "MAHAJANGA" or station_names[i] == "FIANARANTSOA":
		#plt.text(x+0.01,y+0.01,station_names[i],fontsize=3,fontweight='bold', ha='left',va='top')
	#else:
		#plt.text(x+0.01,y+0.01,station_names[i],fontsize=3,fontweight='bold')
	

plt.savefig('madagascar_obs_vs_trmm_min_rainfall_error.png', bbox_inches='tight', dpi=400)




				
					
#for ID in station_IDs:	
	#print ID	
	#print pcp_dict["obs"][ID]
	#print pcp_dict["trmm"][ID]


				
			
	
	



