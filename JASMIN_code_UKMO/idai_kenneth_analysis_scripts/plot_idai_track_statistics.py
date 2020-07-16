import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches


#y1=sys.argv[1]
#y2=sys.argv[2]

ukmodatadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/IDAI/track_statistics/"
ecmwfdatadir = "/gws/nopw/j04/klingaman/emerton/ecmwf/idai_track_statistics/"
#savedir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/IDAI/track_statistics/"

#if not os.path.exists(savedir):
    #os.makedirs(savedir)

#storms_dir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/IDAI/reformatted_idai_forecast_tracks/"


def plot_statistics(stat_type):
	"""stat_type is either "location_error" or "vorticity_bias" """
	fig, ax = plt.subplots()
	fig.set_size_inches(6,6)

	ukmo_avg_data_analysis = np.genfromtxt(ukmodatadir+"idai_average_"+str(stat_type)+"_per_lead_time_vs_analysis.txt")
	ukmo_avg_data_ibtracs = np.genfromtxt(ukmodatadir+"idai_average_" + str(stat_type) + "_per_lead_time_vs_ibtracs.txt")
	
	ecmwf_avg_data_analysis = np.genfromtxt(ecmwfdatadir+"idai_average_"+str(stat_type)+"_per_lead_time_vs_analysis.txt")
	ecmwf_avg_data_ibtracs = np.genfromtxt(ecmwfdatadir+"idai_average_" + str(stat_type) + "_per_lead_time_vs_ibtracs.txt")

	x = np.arange(0,7.25,0.25)
	plt.plot(x, ukmo_avg_data_analysis, color='r')
	plt.plot(x, ukmo_avg_data_ibtracs, color='r', linestyle='--')
	
	plt.plot(x, ecmwf_avg_data_analysis, color='b')
	plt.plot(x, ecmwf_avg_data_ibtracs, color='b', linestyle='--')

	#list_of_files = os.listdir(datadir)
	#ib_pattern = "tr*_average_"+stat_type+"*_ibtracs.txt"
	#an_pattern="tr*_average_"+stat_type+"*_analysis.txt"
	#ib_files_to_plot=[]
	#an_files_to_plot=[]
	#for entry in list_of_files:
		#if fnmatch.fnmatch(entry, ib_pattern):
			#ib_files_to_plot.append(datadir+entry)
		#elif fnmatch.fnmatch(entry, an_pattern):
			#an_files_to_plot.append(datadir + entry)

	#ib_storms_stats = np.genfromtxt(ib_files_to_plot[0])
	#for ff, i in zip(ib_files_to_plot, range(1,len(ib_files_to_plot))):
		#ib_data = np.genfromtxt(ff)
		#ib_storms_stats = np.vstack([ib_storms_stats, ib_data])

	#an_storms_stats = np.genfromtxt(an_files_to_plot[0])
	#for ff, i in zip(an_files_to_plot, range(1,len(an_files_to_plot))):
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
	plt.xticks(fontsize=8)
	plt.yticks(fontsize=8)

	plt.xlabel('Lead Time (Days)', fontsize = 10)
	if stat_type == "location_error":
		plt.ylabel('Track Error (km)', fontsize=10)
		plt.ylim(bottom=0)
	elif stat_type == "mslp_bias":
		plt.ylabel('MSLP Bias (hPa)', fontsize=10)
	elif stat_type == "wind_bias":
		plt.ylabel(r'10m Wind Bias (ms$^{-1}$)', fontsize=10)

	an = plt.Line2D((0, 1), (0, 0), color='k', linewidth=0.5, label = 'Mean error vs. analysis')
	ib = plt.Line2D((0, 1), (0, 0), color='k', linestyle='--', linewidth=0.5, label = 'Mean error vs. IBTrACS')
	ukmo = plt.Line2D((0, 1), (0, 0), color='r', linewidth=0.5, label = 'UKMO')
	ecmwf = plt.Line2D((0, 1), (0, 0), color='b', linewidth=0.5, label = 'ECMWF')
	#spread = mpatches.Patch(color='b', alpha=0.2, label='Spread across storms')
	legend = ax.legend(handles=[an, ib, ukmo, ecmwf], fontsize=8, loc='upper left', title = "Cyclone Idai")
	legend._legend_box.align = "left"
	plt.setp(legend.get_title(), fontsize='8')

	plt.savefig("idai_"+stat_type+".png", dpi=400)

plot_statistics("location_error")
plot_statistics("mslp_bias")
plot_statistics("wind_bias")



