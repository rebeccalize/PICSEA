import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

#Script to plot the location of Mozambique's rainfall observation stations provided by Lelo Tayob (INAM)


station_info=np.zeros((16,3))

print np.shape(station_info)

#0 = station names, 1 = latitude, 2 = longitude, 3 = altitude (m)
station_names = ["Maputo_Mavalane", "Xai-Xai", "Inharrime", "Vilanculos", "Panda", "Inharrime", "Sussundenga", "Beira_Aeroporto", "Chimoio", "Quelimane", "Angoche", "Lumbo", "Nampula", "Montepuez", "Mocimboa da Praia", "Pemba"]
station_info[:,0] = [-25.92, -25.05, -24.38, -22, -24.05, -23.86, -19, -19.8, -19.12, -17.58, -16.22, -15.03, -15.1, -13.13, -11.35, -12.98]
station_info[:,1] = [32.57, 33.63, 35.02, 35.32, 34.72, 35.38, 33.23, 34.54, 33.47, 36.58, 39.9, 40.67, 39.28, 39.03, 40.37, 40.53]
station_info[:,2] = [39,4,43,20,150,14,620,8,731,6,61,10,438,534,27,np.nan]

station_lons = station_info[:,1]
station_lats = station_info[:,0]
station_alts = station_info[:,2]


print station_names
print station_lons
print station_lats


no_stations = len(station_names)

#just madagascar
lat1=-7
lat2=-30
lon1=30
lon2=47

datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"
track_file_list=[]

for root,dirs,files in os.walk(datadir+"mozambique_landfalling"):

	#loop over the directories (i.e. loop over the storms)
	for dir in dirs:
		
		#make a list of all the files in this directory
		list_of_all_files = os.listdir(datadir+"mozambique_landfalling/"+dir)
		pattern="analysis_*.txt"
		#print "directory:", dir
		#for each file in the directory, find the ones that match the forecast filename, analysis filename and ibtracs filename
		for entry in list_of_all_files:
			if fnmatch.fnmatch(entry,pattern):
				track_file_list.append(datadir+"mozambique_landfalling/"+dir+"/"+entry)



for storm in number_of_storms:

	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='h')


	#m.bluemarble(alpha=0.85)
	#m.shadedrelief()
	m.drawcoastlines(linewidth=0.25, color='k')
	m.drawcountries(linewidth=0.25, color='k')
	#m.fillcontinents(color='white')

	


	for i in range(no_stations):
		x, y = m(float(station_lons[i]), float(station_lats[i]))
	
		if station_names[i] == "Pemba":
			plt.text(x+0.01,y+0.01,station_names[i],fontsize=5,fontweight='bold', ha='left',va='top')
			plt.scatter(x,y,6,marker='o',c='k')
		
		
		elif station_names[i] == "Lumbo" or station_names[i] == "Chimoio" or station_names[i] == "Inharrime" or station_names[i] == "Beira_Aeroporto" or station_names[i] == "Maputo_Mavalane":
			plt.text(x+0.01,y+0.01,station_names[i],fontsize=5,fontweight='bold', ha='left',va='top')
		
		elif station_names[i] == "Panda" or station_names[i] == "Xai-Xai" or station_names[i] == "Nampula" or station_names[i] == "Montepuez":
			plt.text(x+0.01,y+0.01,station_names[i],fontsize=5,fontweight='bold', ha='right',va='top')
	
		else:
			plt.text(x+0.01,y+0.01,station_names[i],fontsize=5,fontweight='bold',va='bottom')

	xs,ys = m(station_lons, station_lats)
	im = plt.scatter(xs,ys,6,marker='o',c=station_alts,cmap='viridis_r')


	cbar = fig.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
	cbar.ax.tick_params(labelsize=5)
	cbar.set_label('Altitude [m]', rotation=270, fontsize=5)

	plt.savefig('mozambique_rainfall_stations.png', bbox_inches='tight', dpi=400)
	
	
	
	
	
	
	
	
	
	
	
