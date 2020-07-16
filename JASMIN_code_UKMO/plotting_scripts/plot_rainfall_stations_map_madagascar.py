import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
import csv

#Script to plot the location of Madagascar's rainfall observation stations


obs_data = np.genfromtxt("/home/users/emerton/MADAGASCAR_RR_daily_2006_2018_CSV.csv", delimiter=',',dtype=str)
	
station_names = obs_data[0,1:]
station_IDs = obs_data[1,1:]
station_lons = obs_data[2,1:]
station_lats = obs_data[3,1:]

precip = obs_data[4:,:] #first value of each row [:,0] is the date! each column [x,:] is a different station


print station_names
print station_IDs
print station_lons
print station_lats


no_stations = len(station_IDs)

#just madagascar
lat1=-9 
lat2=-28
lon1=42
lon2=52.5

fig = plt.figure(figsize=(6, 3))
ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='h')


#m.bluemarble(alpha=0.85)
m.shadedrelief()
m.drawcoastlines(linewidth=0.4, color='k')
m.drawcountries(linewidth=0.4, color='k')
#m.fillcontinents(color='white')




for i in range(no_stations):
	x, y = m(float(station_lons[i]), float(station_lats[i]))
	plt.scatter(x,y,4,marker='o',color='b')
	if station_names[i] == "ANTANANARIVO" or station_names[i] == "MAHANORO" or station_names[i] == "MAHAJANGA" or station_names[i] == "FIANARANTSOA":
		plt.text(x+0.01,y+0.01,station_names[i],fontsize=3,fontweight='bold', ha='left',va='top')
	else:
		plt.text(x+0.01,y+0.01,station_names[i],fontsize=3,fontweight='bold')
	

plt.savefig('madagascar_rainfall_stations.png', bbox_inches='tight', dpi=400)
