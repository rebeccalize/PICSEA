import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches


year1s=[2017, 2018]
year2s=[2018, 2019]

#datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/SIO_storms/"
savedir = "./"#"/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/plots/"

if not os.path.exists(savedir):
    os.makedirs(savedir)


def plot_statistics(stat_type):
	"""stat_type is either "location_error" or "mslp_bias" or "wind_bias """
	fig, ax = plt.subplots()
	fig.set_size_inches(6,6)


	x = np.arange(0,7.25,0.25)

	weak_depressions=[]
	tropical_depressions=[]
	moderate_tropical_storms=[]
	strong_tropical_storms=[]
	tropical_cyclones=[]
	intense_tropical_cyclones=[]
	very_intense_tropical_cyclones=[]
	
	
	
	for y1,y2 in zip(year1s,year2s):
		
		datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/SIO_storms/"
		ibtrackdir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/"
		ib_dir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/"
		
		list_of_files = os.listdir(datadir)
		stats_pattern = ib_pattern = "tr*_average_"+stat_type+"*_ibtracs.txt"
		
		for entry in list_of_files:
			
			if fnmatch.fnmatch(entry,stats_pattern):
				
				#for this storm's statistics, find the ibtracs track data that goes with that storm so we can find the max intensity and get the category 
				ibtrack_pattern = entry[0:6]
				ibtrack_file = ib_dir+entry[0:6]+"/ibtracs_"+entry[0:6]+"_reformatted.txt"
				
				ib_data = np.genfromtxt(ibtrack_file,dtype=float,skip_header=1)
				ib_wind = ib_data[:,10]
				
				for i in range(len(ib_wind)):
					if ib_wind[i] > 10000:
						ib_wind[i] = np.nan
						
				max_wind = np.nanmax(ib_wind)*3.6
				
				if max_wind < 51:
					c = 'khaki'
					ls = 'dotted'
					weak_depressions.append(datadir+entry)
				elif 51<= max_wind <63:
					c = 'gold'
					ls = 'dotted'
					tropical_depressions.append(datadir+entry)
				elif 63 <= max_wind < 89:
					c = 'darkorange'
					ls = 'dotted'
					moderate_tropical_storms.append(datadir+entry)
				elif 89 <= max_wind < 118:
					c = 'black'
					ls = 'dotted'
					strong_tropical_storms.append(datadir+entry)
				elif 118 <= max_wind < 166:
					c = 'orangered'
					ls = '-'
					tropical_cyclones.append(datadir+entry)
				elif 166 <= max_wind < 213:
					c = 'firebrick'
					ls ='-'
					intense_tropical_cyclones.append(datadir+entry)
				elif max_wind >= 213:
					c = 'k'
					ls = '-'
					very_intense_tropical_cyclones.append(datadir+entry)
					
					
	if len(weak_depressions) > 0:
		weak_depression_stats = np.genfromtxt(weak_depressions[0])
		for ff, i in zip(weak_depressions, range(1,len(weak_depressions))):
			data = np.genfromtxt(ff)
			weak_depression_stats = np.vstack([weak_depression_stats,data])
		weak_depression_stats[weak_depression_stats==0] = np.nan
		
	if len(tropical_depressions) > 0:
		tropical_depression_stats = np.genfromtxt(tropical_depressions[0])
		for ff, i in zip(tropical_depressions, range(1,len(tropical_depressions))):
			data = np.genfromtxt(ff)
			tropical_depression_stats = np.vstack([tropical_depression_stats,data])
		tropical_depression_stats[tropical_depression_stats==0] = np.nan
	
	if len(moderate_tropical_storms) > 0:
		moderate_tropical_storm_stats = np.genfromtxt(moderate_tropical_storms[0])
		for ff, i in zip(moderate_tropical_storms, range(1,len(moderate_tropical_storms))):
			data = np.genfromtxt(ff)
			moderate_tropical_storm_stats = np.vstack([moderate_tropical_storm_stats,data])
		moderate_tropical_storm_stats[moderate_tropical_storm_stats==0] = np.nan
	else:
		moderate_tropical_storm_stats=np.zeros((2,29))
		moderate_tropical_storm_stats[:,:] = np.nan
		
	if len(strong_tropical_storms) > 0:	
		strong_tropical_storm_stats = np.genfromtxt(strong_tropical_storms[0])
		for ff, i in zip(strong_tropical_storms, range(1,len(strong_tropical_storms))):
			data = np.genfromtxt(ff)
			strong_tropical_storm_stats = np.vstack([strong_tropical_storm_stats,data])
		strong_tropical_storm_stats[strong_tropical_storm_stats==0] = np.nan
	
	if len(tropical_cyclones) > 0:	
		tropical_cyclone_stats = np.genfromtxt(tropical_cyclones[0])
		for ff, i in zip(tropical_cyclones, range(1,len(tropical_cyclones))):
			data = np.genfromtxt(ff)
			tropical_cyclone_stats = np.vstack([tropical_cyclone_stats,data])
		tropical_cyclone_stats[tropical_cyclone_stats==0] = np.nan
		
	if len(intense_tropical_cyclones) > 0:
		intense_tropical_cyclone_stats = np.genfromtxt(intense_tropical_cyclones[0])
		for ff, i in zip(intense_tropical_cyclones, range(1,len(intense_tropical_cyclones))):
			data = np.genfromtxt(ff)
			intense_tropical_cyclone_stats = np.vstack([intense_tropical_cyclone_stats,data])
		intense_tropical_cyclone_stats[intense_tropical_cyclone_stats==0] = np.nan
		
	if len(very_intense_tropical_cyclones) > 0:
		very_intense_tropical_cyclone_stats = np.genfromtxt(very_intense_tropical_cyclones[0])
		for ff, i in zip(very_intense_tropical_cyclones, range(1,len(very_intense_tropical_cyclones))):
			data = np.genfromtxt(ff)
			very_intense_tropical_cyclone_stats = np.vstack([very_intense_tropical_cyclone_stats,data])
		very_intense_tropical_cyclone_stats[very_intense_tropical_cyclone_stats==0] = np.nan
		
	
	
	weak_depressions_mean, tropical_depressions_mean, moderate_tropical_storms_mean, strong_tropical_storms_mean, tropical_cyclones_mean, intense_tropical_cyclones_mean, very_intense_tropical_cyclones_mean = (np.zeros(29) for i in range(7))
	for lt in range(29):
		#weak_depressions_mean[lt] = np.nanmean(weak_depression_stats[:,lt])
		#tropical_depressions_mean[lt] = np.nanmean(tropical_depression_stats[:,lt])
		moderate_tropical_storms_mean[lt] = np.nanmean(moderate_tropical_storm_stats[:,lt])
		strong_tropical_storms_mean[lt] = np.nanmean(strong_tropical_storm_stats[:,lt])
		tropical_cyclones_mean[lt] = np.nanmean(tropical_cyclone_stats[:,lt])
		intense_tropical_cyclones_mean[lt] = np.nanmean(intense_tropical_cyclone_stats[:,lt])
		very_intense_tropical_cyclones_mean[lt] = np.nanmean(very_intense_tropical_cyclone_stats[:,lt])
		
	print intense_tropical_cyclones_mean
	
	
	

	

	#ib_maxline, ib_minline, an_maxline, an_minline = (np.zeros(29) for i in range(4))
	#for lt in range(29):
		#ib_maxline[lt] = np.nanmax(ib_storms_stats[:,lt])
		#ib_minline[lt] = np.nanmin(ib_storms_stats[:,lt])
		#an_maxline[lt] = np.nanmax(an_storms_stats[:, lt])
		#an_minline[lt] = np.nanmin(an_storms_stats[:, lt])

	#ax.fill_between(x, ib_minline, ib_maxline, facecolor='b', edgecolor = 'b', alpha=0.2,linestyle='--')
	#ax.fill_between(x, an_minline, an_maxline, facecolor='b', edgecolor='b', alpha=0.2)
	
	#plt.plot(x, weak_depressions_mean, color='khaki', linestyle='dotted')
	#plt.plot(x, tropical_depressions_mean, color='gold', linestyle='dotted')
	plt.plot(x, moderate_tropical_storms_mean, color='darkorange', linestyle='dotted')
	plt.plot(x, strong_tropical_storms_mean, color='black', linestyle='dotted')
	plt.plot(x, tropical_cyclones_mean, color='orangered')
	plt.plot(x, intense_tropical_cyclones_mean, color='firebrick')
	plt.plot(x, very_intense_tropical_cyclones_mean, color='black')

	plt.xlim(0,7)
	plt.xticks(fontsize=8)
	plt.yticks(fontsize=8)

	plt.xlabel('Lead Time (Days)', fontsize = 10)
	if stat_type == "location_error":
		plt.ylabel('Track Error (km)', fontsize=10)
		plt.ylim(bottom=0)
	elif stat_type == "mslp_bias":
		plt.ylabel('MSLP Bias (hPa)', fontsize=10)
		plt.hlines(0, x[0], x[-1], color='silver',linewidth=0.5)
	elif stat_type == "wind_bias":
		plt.ylabel(r'10m Wind Bias (ms$^{-1}$)', fontsize=10)
		plt.hlines(0, x[0], x[-1], color='silver',linewidth=0.5)

	#legend
	tdf = plt.Line2D((0, 1), (0, 0), color='khaki',linestyle='dotted',linewidth=1)
	td = plt.Line2D((0, 1), (0, 0), color='yellow', linestyle='dotted',linewidth=1)
	tsm = plt.Line2D((0, 1), (0, 0), color='darkorange', linestyle='dotted',linewidth=1)
	tsf = plt.Line2D((0, 1), (0, 0), color='black',linestyle='dotted',linewidth=1)
	tc = plt.Line2D((0, 1), (0, 0), color='orangered',linewidth=1)
	tci = plt.Line2D((0, 1), (0, 0), color='firebrick',linewidth=1)
	tcti = plt.Line2D((0, 1), (0, 0), color='black',linewidth=1)
	#title = "2010-2018\n"+str(no_tracks)+" Cyclones\nCategory:"
	legend = ax.legend((tsm, tsf, tc, tci, tcti), ['Moderate Tropical Storm ('+str(len(moderate_tropical_storms))+')', 'Strong Tropical Storm ('+str(len(strong_tropical_storms))+')', 'Tropical Cyclone ('+str(len(tropical_cyclones))+')', 'Intense Tropical Cyclone ('+str(len(intense_tropical_cyclones))+')', 'Very Intense Tropical Cyclone ('+str(len(very_intense_tropical_cyclones))+')'], fontsize=10)
	#'Weak Depression ('+str(len(weak_depressions))+')', 'Tropical Depression ('+str(len(tropical_depressions))+')', 
	
	#plt.setp(legend.get_title(), fontsize='5')

	plt.savefig("ukmo_hres_"+str(year1s[0])+"-"+str(year2s[-1])+"_"+stat_type+"_by_intensity_category.png", dpi=400, bbox_inches='tight')

plot_statistics("location_error")
plot_statistics("mslp_bias")
plot_statistics("wind_bias")



