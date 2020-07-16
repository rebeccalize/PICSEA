import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches



datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"
savedir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"


def plot_multiregion_stats(stat_type):
	"""stat_type is either "location_error" or "mslp_bias" or "wind_bias" """
	fig, ax = plt.subplots()
	fig.set_size_inches(6,6)


	moz_ib_error_stats = np.genfromtxt(datadir+"mozambique_landfalling/average_" + str(stat_type) + "_per_lead_time_vs_ibtracs.txt")
	moz_an_error_stats = np.genfromtxt(datadir+"mozambique_landfalling/average_" + str(stat_type) + "_per_lead_time_vs_analysis.txt")

	mada_ib_error_stats = np.genfromtxt(datadir+"madagascar_landfalling/average_" + str(stat_type) + "_per_lead_time_vs_ibtracs.txt")
	mada_an_error_stats = np.genfromtxt(datadir+"madagascar_landfalling/average_" + str(stat_type) + "_per_lead_time_vs_analysis.txt")
	
	sey_ib_error_stats = np.genfromtxt(datadir+"seychelles_landfalling/average_" + str(stat_type) + "_per_lead_time_vs_ibtracs.txt")
	sey_an_error_stats = np.genfromtxt(datadir+"seychelles_landfalling/average_" + str(stat_type) + "_per_lead_time_vs_analysis.txt")

	

	#x = np.arange(0,7)
	x = np.linspace(0,7,29)

	plt.plot(x, moz_ib_error_stats, color = '#003f5c', linestyle='--')
	#plt.plot(x, moz_an_error_stats, color='#003f5c')

	plt.plot(x, mada_ib_error_stats, color = '#bc5090', linestyle='--')
	#plt.plot(x, mada_an_error_stats, color='#bc5090')
	
	plt.plot(x, sey_ib_error_stats, color = '#ffa600', linestyle='--')
	#plt.plot(x, sey_an_error_stats, color = '#ffa600')





	#the following plots the spread of the errors across all the individual storms
	#ib_pattern = "tr*_average_" + stat_type + "*_ibtracs.txt"
	#an_pattern = "tr*_average_" + stat_type + "*_analysis.txt"
	#ib_files=[]
	#an_files=[]
	#for y1, y2, in zip(year1array, year2array):
		#list_of_files = os.listdir(datadir+str(y1)+"_"+str(y2))
		#for entry in list_of_files:
			#if fnmatch.fnmatch(entry, ib_pattern):
				#ib_files.append(datadir+str(y1)+"_"+str(y2)+"/"+entry)
			#elif fnmatch.fnmatch(entry, an_pattern):
				#an_files.append(datadir+str(y1)+"_"+str(y2)+"/"+entry)


	#ib_storms_stats = np.genfromtxt(ib_files[0])
	#for ff, i in zip(ib_files, range(1,len(ib_files))):
		#ib_data = np.genfromtxt(ff)
		#ib_storms_stats = np.vstack([ib_storms_stats, ib_data])

	#an_storms_stats = np.genfromtxt(an_files[0])
	#for ff, i in zip(an_files, range(1,len(an_files))):
		#an_data = np.genfromtxt(ff)
		#an_storms_stats = np.vstack([an_storms_stats, an_data])

	#ib_storms_stats[ib_storms_stats==0] = np.nan
	#an_storms_stats[an_storms_stats == 0] = np.nan

	#ib_maxline, ib_minline, an_maxline, an_minline = (np.zeros(29) for i in range(4))
	#for lt in range(29):
		#ib_maxline[lt] = np.nanmax(ib_storms_stats[:,lt])
		#ib_minline[lt] = np.nanmin(ib_storms_stats[:,lt])
		#an_maxline[lt] = np.nanmax(an_storms_stats[:, lt])
		#an_minline[lt] = np.nanmin(an_storms_stats[:, lt])

	#ax.fill_between(x, ib_minline, ib_maxline, facecolor='b', edgecolor = 'b', alpha=0.2,linestyle='--')
	#ax.fill_between(x, an_minline, an_maxline, facecolor='b', edgecolor='b', alpha=0.2)

	plt.xlim(0,7)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)

	an = plt.Line2D((0, 1), (0, 0), color='k', linewidth=0.5, label = 'Mean error vs. analysis')
	ib = plt.Line2D((0, 1), (0, 0), color='k', linestyle='--', linewidth=0.5, label = 'Mean error vs. IBTrACS')
	gap = plt.Line2D((0, 1), (0, 0), color='white', linewidth=0.5, label = ' ')
	moz = plt.Line2D((0, 1), (0, 0), color='#003f5c', linewidth=0.5, label='Mozambique, 11')
	mada = plt.Line2D((0, 1), (0, 0), color='#bc5090', linewidth=0.5, label='Madagascar, 35')
	sey = plt.Line2D((0, 1), (0, 0), color='#ffa600', linewidth=0.5, label='Seychelles, 4')
	#N320 = plt.Line2D((0, 1), (0, 0), color='#BDBDBD', linewidth=0.5, label='N320 (~40km) 01/07/2006 - 09/03/2010')
	#N512 = plt.Line2D((0, 1), (0, 0), color='#757575', linewidth=0.5, label='N512 (~25km) 09/03/2010 - 15/07/2014')
	#N768 = plt.Line2D((0, 1), (0, 0), color='#424242', linewidth=0.5, label='N768 (~17km) 15/07/2014 - 11/07/2017')
	#N1280 = plt.Line2D((0, 1), (0, 0), color='r', linewidth=0.5, label='N1280 (~10km) 11/07/2017 - present')
	#spread = mpatches.Patch(color='b', alpha=0.2, label='Spread across storms')
	title="South-West Indian Ocean\nLandfalling Cyclones 2006-2016\n"


	plt.xlabel('Lead Time (Days)', fontsize = 15)
	if stat_type == "location_error":
		plt.ylabel('Track Error (km)', fontsize=15)
		plt.ylim(0,1000)
		#legend = ax.legend(handles=[an, ib, gap, N320, N512, N768,N1280], fontsize=8, loc='upper left', title=title)
		#legend = ax.legend(handles=[an, ib, gap, mada,moz,sey], fontsize=8, loc='upper left', title=title)
		#legend = ax.legend(handles=[mada,moz,sey], fontsize=15, loc='upper left', title=title)
		#legend._legend_box.align = "left"
		#plt.setp(legend.get_title(), fontsize='15')
	elif stat_type == "mslp_bias":
		plt.ylabel('MSLP Bias (hPa)', fontsize=15)
		plt.ylim(-20,20)
		#legend = ax.legend(handles=[an, ib, gap, N320, N512, N768,N1280], fontsize=8, loc='lower right', title=title)
		#legend = ax.legend(handles=[an, ib, gap, mada,moz,sey], fontsize=8, loc='lower left', title=title)
		#legend = ax.legend(handles=[mada,moz,sey], fontsize=8, loc='upper left', title=title)
		#legend._legend_box.align = "left"
		#plt.setp(legend.get_title(), fontsize='8')
	elif stat_type == "wind_bias":
		plt.ylabel(r'10m Wind Bias (ms$^{-1}$)', fontsize=15)
		plt.ylim(-15,15)
		#legend = ax.legend(handles=[an, ib, gap, N320, N512, N768,N1280], fontsize=8, loc='upper left', title=title)
		#legend = ax.legend(handles=[an, ib, gap, mada,moz,sey], fontsize=8, loc='upper left', title=title)
		legend = ax.legend(handles=[mada,moz,sey], fontsize=15, loc='upper left', title=title)
		legend._legend_box.align = "left"
		plt.setp(legend.get_title(), fontsize='15')



	plt.savefig("ukmo_nwp_landfalling_cyclones_average_"+stat_type+"_wlegend_ibtracsonly_poster.png", dpi=400, bbox_inches='tight')


#storm_dirs=[]
#for y1,y2 in zip(res1year1,res1year2):
	#storms_dir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/"
	#for root,dirs,files in os.walk(storms_dir):
		#for dir in dirs:
			#storm_dirs.append(dir)
#NS = len(storm_dirs) #number of storms

plot_multiregion_stats("location_error")
plot_multiregion_stats("mslp_bias")
plot_multiregion_stats("wind_bias")




