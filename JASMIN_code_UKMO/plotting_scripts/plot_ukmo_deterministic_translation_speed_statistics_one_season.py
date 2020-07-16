import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches


y1=sys.argv[1]
y2=sys.argv[2]

datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/translation_speed_errors/"+str(y1)+"_"+str(y2)+"/SIO_storms/"
savedir = "./"#"/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/plots/"

if not os.path.exists(savedir):
    os.makedirs(savedir)

storms_dir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/"
storm_dirs=[]
for root,dirs,files in os.walk(storms_dir):
	for dir in dirs:
		storm_dirs.append(dir)
NS = len(storm_dirs) #number of storms

def plot_statistics(stat_type, y1, y2):
	"""stat_type is either "error" or "bias" """
	fig, ax = plt.subplots()
	fig.set_size_inches(6,6)

	year_data_analysis = np.genfromtxt(datadir+"average_translation_speed_"+str(stat_type)+"_per_lead_time_vs_analysis.txt")
	year_data_ibtracs = np.genfromtxt(datadir+"average_translation_speed_" + str(stat_type) + "_per_lead_time_vs_ibtracs.txt")

	x = np.arange(0,7.25,0.25)
	plt.plot(x, year_data_analysis, color='k')
	plt.plot(x, year_data_ibtracs, color='k', linestyle='--')

	list_of_files = os.listdir(datadir)
	ib_pattern = "tr*_average_translation_speed_"+stat_type+"*_ibtracs.txt"
	an_pattern="tr*_average_translation_speed_"+stat_type+"*_analysis.txt"
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

	ax.fill_between(x, ib_minline, ib_maxline, facecolor='b', edgecolor = 'b', alpha=0.2,linestyle='--')
	ax.fill_between(x, an_minline, an_maxline, facecolor='b', edgecolor='b', alpha=0.2)

	plt.xlim(0,7)
	plt.xticks(fontsize=8)
	plt.yticks(fontsize=8)

	plt.xlabel('Lead Time (Days)', fontsize = 10)
	if stat_type == "error":
		plt.ylabel(r'Translation Speed Error (ms$^{-1}$)', fontsize=10)
		plt.ylim(bottom=0)
	elif stat_type == "bias":
		plt.ylabel(r'Translation Speed Bias (ms$^{-1}$)', fontsize=10)
s

	an = plt.Line2D((0, 1), (0, 0), color='k', linewidth=0.5, label = 'Mean error vs. analysis')
	ib = plt.Line2D((0, 1), (0, 0), color='k', linestyle='--', linewidth=0.5, label = 'Mean error vs. IBTrACS')
	spread = mpatches.Patch(color='b', alpha=0.2, label='Spread across storms')
	legend = ax.legend(handles=[an, ib, spread], fontsize=8, loc='upper left', title = y1+"-"+y2+", "+str(NS)+" Cyclones")
	legend._legend_box.align = "left"
	plt.setp(legend.get_title(), fontsize='8')

	plt.savefig(savedir+y1+"_"+y2+"_"+stat_type+".png", dpi=400)

plot_statistics("error", y1,y2)
plot_statistics("bias", y1, y2)




