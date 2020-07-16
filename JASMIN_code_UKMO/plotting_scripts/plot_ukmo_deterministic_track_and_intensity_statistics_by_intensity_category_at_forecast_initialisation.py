import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches




datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/intensity_category/"
savedir = "./"#"/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/plots/"

if not os.path.exists(savedir):
    os.makedirs(savedir)


def plot_statistics(stat_type):
	"""stat_type is either "location_error" or "vorticity_bias" """
	fig, ax = plt.subplots()
	fig.set_size_inches(6,6)

	x = np.arange(0,7.25,0.25)
	
	if stat_type == "sample_size":
	
		pre_tc_data = np.genfromtxt(datadir+"2010_2019.pre-TC_number_of_forecasts_included_vs_ibtracs.txt")
		moderate_tropical_storm_data = np.genfromtxt(datadir+"2010_2019.moderate_tropical_storm_number_of_forecasts_included_vs_ibtracs.txt")
		strong_tropical_storm_data = np.genfromtxt(datadir+"2010_2019.strong_tropical_storm_number_of_forecasts_included_vs_ibtracs.txt")
		tropical_cyclone_data = np.genfromtxt(datadir+"2010_2019.tropical_cyclone_number_of_forecasts_included_vs_ibtracs.txt")
		intense_tropical_cyclone_data = np.genfromtxt(datadir+"2010_2019.intense_tropical_cyclone_number_of_forecasts_included_vs_ibtracs.txt")
		very_intense_tropical_cyclone_data = np.genfromtxt(datadir+"2010_2019.very_intense_tropical_cyclone_number_of_forecasts_included_vs_ibtracs.txt")
		
		all_tropical_storms_data = np.genfromtxt(datadir+"2010_2019.all_tropical_storms_number_of_forecasts_included_vs_ibtracs.txt")
		intense_and_vintense_tropical_cyclone_data = np.genfromtxt(datadir+"2010_2019.intense_and_very_intense_tropical_cyclones_number_of_forecasts_included_vs_ibtracs.txt")
		
		plt.plot(x, moderate_tropical_storm_data[0:29], color='darkorange', linestyle='dotted')
		plt.plot(x, strong_tropical_storm_data[0:29], color='black', linestyle='dotted')
		plt.plot(x, tropical_cyclone_data[0:29], color = 'orangered')
		plt.plot(x, intense_tropical_cyclone_data[0:29], color = 'firebrick')
		plt.plot(x, very_intense_tropical_cyclone_data[0:29], color = 'black')
		
		#plt.plot(x, all_tropical_storms_data[0:29], color='black', linestyle = 'dotted')
		#plt.plot(x, tropical_cyclone_data[0:29], color = 'orangered')
		plt.plot(x, pre_tc_data[0:29], color = 'skyblue', linestyle = 'dotted')
		#plt.plot(x, intense_and_vintense_tropical_cyclone_data[0:29], color='darkred')
		
	else:
	
		pre_tc_data = np.genfromtxt(datadir+"2010_2019.pre-TC.average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")
		moderate_tropical_storm_data = np.genfromtxt(datadir+"2010_2019.moderate_tropical_storm.average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")
		strong_tropical_storm_data = np.genfromtxt(datadir+"2010_2019.strong_tropical_storm.average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")
		tropical_cyclone_data = np.genfromtxt(datadir+"2010_2019.tropical_cyclone.average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")
		intense_tropical_cyclone_data = np.genfromtxt(datadir+"2010_2019.intense_tropical_cyclone.average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")
		very_intense_tropical_cyclone_data = np.genfromtxt(datadir+"2010_2019.very_intense_tropical_cyclone.average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")
		
		all_tropical_storms_data = np.genfromtxt(datadir+"2010_2019.all_tropical_storms.average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")
		intense_and_vintense_tropical_cyclone_data = np.genfromtxt(datadir+"2010_2019.intense_and_very_intense_tropical_cyclones.average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")
	
		plt.plot(x, moderate_tropical_storm_data[0:29,0], color='darkorange', linestyle='dotted')
		plt.plot(x, strong_tropical_storm_data[0:29,0], color='black', linestyle='dotted')
		plt.plot(x, tropical_cyclone_data[0:29,0], color = 'orangered')
		plt.plot(x, intense_tropical_cyclone_data[0:29,0], color = 'firebrick')
		plt.plot(x, very_intense_tropical_cyclone_data[0:29,0], color = 'black')
		
		#plt.plot(x, tropical_cyclone_data[0:29,0], color = 'orangered')
		#plt.plot(x, all_tropical_storms_data[0:29,0], color='black', linestyle = 'dotted')
		plt.plot(x, pre_tc_data[0:29,0], color = 'skyblue', linestyle='dotted')
		#plt.plot(x, intense_and_vintense_tropical_cyclone_data[0:29,0], color='darkred')
	
		ax.fill_between(x, moderate_tropical_storm_data[0:29,1], moderate_tropical_storm_data[0:29,2], facecolor='darkorange', edgecolor = None, alpha=0.2)
		ax.fill_between(x, strong_tropical_storm_data[0:29,1], strong_tropical_storm_data[0:29,2], facecolor='black', edgecolor = None, alpha=0.2)
		ax.fill_between(x, tropical_cyclone_data[0:29,1], tropical_cyclone_data[0:29,2], facecolor='orangered', edgecolor = None, alpha=0.2)
		ax.fill_between(x, intense_tropical_cyclone_data[0:29,1], intense_tropical_cyclone_data[0:29,2], facecolor='firebrick', edgecolor = None, alpha=0.2)
		ax.fill_between(x, very_intense_tropical_cyclone_data[0:29,1], very_intense_tropical_cyclone_data[0:29,2], facecolor='black', edgecolor = None, alpha=0.2)
		
		#ax.fill_between(x, all_tropical_storms_data[0:29,1], all_tropical_storms_data[0:29,2], facecolor='black', edgecolor = None, alpha=0.2)
		ax.fill_between(x, pre_tc_data[0:29,1], pre_tc_data[0:29,2], facecolor='skyblue', edgecolor = None, alpha=0.2)
		#ax.fill_between(x, tropical_cyclone_data[0:29,1], tropical_cyclone_data[0:29,2], facecolor='orangered', edgecolor = None, alpha=0.2)
		#ax.fill_between(x, intense_and_vintense_tropical_cyclone_data[0:29,1], intense_and_vintense_tropical_cyclone_data[0:29,2], facecolor='darkred', edgecolor = None, alpha=0.2)
		
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
	elif stat_type == "sample_size":
		plt.ylabel('Number of Forecasts',fontsize=10)

	#legend
	ptc = plt.Line2D((0, 1), (0, 0), color='skyblue', linestyle='dotted',linewidth=1)
	tsm = plt.Line2D((0, 1), (0, 0), color='darkorange', linestyle='dotted',linewidth=1)
	tsf = plt.Line2D((0, 1), (0, 0), color='black',linestyle='dotted',linewidth=1)
	ats = plt.Line2D((0, 1), (0, 0), color='black', linestyle='dotted',linewidth=1)
	tc = plt.Line2D((0, 1), (0, 0), color='orangered',linewidth=1)
	tci = plt.Line2D((0, 1), (0, 0), color='firebrick',linewidth=1)
	tcti = plt.Line2D((0, 1), (0, 0), color='black',linewidth=1)
	ivtc = plt.Line2D((0, 1), (0, 0), color='darkred',linewidth=1)
	#title = "2010-2018\n"+str(no_tracks)+" Cyclones\nCategory:"
	legend = ax.legend((ptc,tsm, tsf, tc, tci, tcti), ['Tropical Depression','Moderate Tropical Storm', 'Strong Tropical Storm', 'Tropical Cyclone', 'Intense Tropical Cyclone', 'Very Intense Tropical Cyclone'], fontsize=10)
	#legend = ax.legend((ptc, ats, tc, ivtc), ['Tropical Depressions', 'Tropical Storms', 'Tropical Cyclones','Intense Tropical Cyclones'], fontsize=10)
	#'Weak Depression ('+str(len(weak_depressions))+')', 'Tropical Depression ('+str(len(tropical_depressions))+')',

	plt.savefig("ukmo_hres.2010-2019."+stat_type+"_by_intensity_category_on_day_forecast_initialised.png", dpi=400, bbox_inches='tight')

plot_statistics("location_error")
plot_statistics("mslp_bias")
plot_statistics("wind_bias")
plot_statistics("sample_size")



