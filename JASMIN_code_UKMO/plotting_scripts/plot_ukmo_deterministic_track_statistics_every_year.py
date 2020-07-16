import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches


y1s=[2010,2011,2012,2013,2014,2015,2016,2017,2018]
y2s=[2011,2012,2013,2014,2015,2016,2017,2018,2019]


savedir = "./"#"/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/plots/"



def plot_statistics(stat_type):

	fig, ax = plt.subplots()
	fig.set_size_inches(6,6)
	
	
	for y1,y2 in zip(y1s,y2s):

		datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/SIO_storms/"


		storms_dir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/"
		storm_dirs=[]
		for root,dirs,files in os.walk(storms_dir):
			for dir in dirs:
				storm_dirs.append(dir)
		NS = len(storm_dirs) #number of storms
	

		year_data_analysis = np.genfromtxt(datadir+"average_"+str(stat_type)+"_per_lead_time_vs_analysis.txt")
		year_data_ibtracs = np.genfromtxt(datadir+"average_" + str(stat_type) + "_per_lead_time_vs_ibtracs.txt")

		x = np.arange(0,7.25,0.25)
		#plt.plot(x, year_data_analysis, color='k')
		plt.plot(x, year_data_ibtracs[0:29], linestyle='--')

		list_of_files = os.listdir(datadir)
		ib_pattern = "tr*_average_"+stat_type+"*_ibtracs.txt"
		an_pattern="tr*_average_"+stat_type+"*_analysis.txt"
		ib_files_to_plot=[]
		an_files_to_plot=[]
		for entry in list_of_files:
			if fnmatch.fnmatch(entry, ib_pattern):
				ib_files_to_plot.append(datadir+entry)
			elif fnmatch.fnmatch(entry, an_pattern):
				an_files_to_plot.append(datadir + entry)

		ib_storms_stats = np.genfromtxt(ib_files_to_plot[0])
		for ff, i in zip(ib_files_to_plot, range(1,len(ib_files_to_plot))):
			ib_data = np.genfromtxt(ff)
			ib_storms_stats = np.vstack([ib_storms_stats, ib_data])

		an_storms_stats = np.genfromtxt(an_files_to_plot[0])
		for ff, i in zip(an_files_to_plot, range(1,len(an_files_to_plot))):
			an_data = np.genfromtxt(ff)
			an_storms_stats = np.vstack([an_storms_stats, an_data])

		ib_storms_stats[ib_storms_stats==0] = np.nan
		an_storms_stats[an_storms_stats == 0] = np.nan

		ib_maxline, ib_minline, an_maxline, an_minline = (np.zeros(29) for i in range(4))
		for lt in range(29):
			ib_maxline[lt] = np.nanmax(ib_storms_stats[:,lt])
			ib_minline[lt] = np.nanmin(ib_storms_stats[:,lt])
			an_maxline[lt] = np.nanmax(an_storms_stats[:, lt])
			an_minline[lt] = np.nanmin(an_storms_stats[:, lt])

		ax.fill_between(x, ib_minline, ib_maxline, alpha=0.2,linestyle='--') #, facecolor='b', edgecolor = 'b',
		#ax.fill_between(x, an_minline, an_maxline, alpha=0.2) #, facecolor='b', edgecolor='b'

	plt.xlim(0,7)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)

	plt.xlabel('Lead Time (Days)', fontsize = 12)
	if stat_type == "location_error":
		plt.ylabel('Track Error (km)', fontsize=12)
		plt.ylim(bottom=0)
	elif stat_type == "mslp_bias":
		plt.ylabel('MSLP Bias (hPa)', fontsize=12)
	elif stat_type == "wind_bias":
		plt.ylabel(r'10m Wind Bias (ms$^{-1}$)', fontsize=12)

	an = plt.Line2D((0, 1), (0, 0), color='k', linewidth=0.5, label = 'Mean error vs. analysis')
	ib = plt.Line2D((0, 1), (0, 0), color='k', linestyle='--', linewidth=0.5, label = 'Mean error vs. IBTrACS')
	spread = mpatches.Patch(color='k', alpha=0.2, label='Spread across storms')
	legend = ax.legend(handles=[ib, spread], fontsize=8, loc='upper left')
	legend._legend_box.align = "left"

	plt.savefig(savedir+"ukmo_det_2010-2019_each_year_"+stat_type+".png", dpi=400)

plot_statistics("location_error")
#plot_statistics("mslp_bias")
#plot_statistics("wind_bias")



