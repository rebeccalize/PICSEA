import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches




datadir = "/perm/mo/more/picsea/ecmwf_track_and_intensity_errors/"



def plot_statistics(stat_type, obs_track):

	
	"""stat_type is either "location_error" or "mslp_bias" or wind_speed_bias or sample_size"""
	fig, ax = plt.subplots()
	fig.set_size_inches(6,6)

	x = np.arange(0,7.25,0.25)
	
	#c="#FCA311"
	
	for fcst_centre in ["ecmwf", "ukmo"]:
	
		if fcst_centre == "ecmwf":
			datadir = "/perm/mo/more/picsea/ecmwf_track_and_intensity_errors/"
			c='blue'
			timeperiod='2012-2020'
			
		elif fcst_centre == "ukmo":
			datadir = "/home/mo/more/PICSEA/UKMO_comparison_stats/"
			c = '#ffa600'
			timeperiod='2010-2020'
			
		for fcst_type in ["det", "ctrl", "mean"]:
	
			if fcst_type == 'det':
				ls = '-'
			elif fcst_type == 'ctrl':
				ls = '-.'
			elif fcst_type == 'mean':
				ls = '--'
			elif fcst_type == 'EPS':
				ls = ':'
				
					
			if stat_type == "sample_size":
		
			
				data = np.genfromtxt(datadir+timeperiod+"number_of_forecasts_included_vs_"+obs_track+"_"+fcst_type+".txt")
			
				plt.plot(x, data[0:29], color=c, linestyle=ls)
			
		
			else:

			
				data = np.genfromtxt(datadir+timeperiod+"average_"+stat_type+"_per_lead_time_vs_"+obs_track+"_"+fcst_type+"_with_confidence_intervals_95.txt")
	
			
				if stat_type == "wind_bias":
			
					plt.plot(x, data[0:29,0]*3.6, color=c, linestyle=ls)

	
					ax.fill_between(x, data[0:29,1]*3.6, data[0:29,2]*3.6, facecolor=c, edgecolor = None, alpha=0.4)
			
			
				else:
					plt.plot(x, data[0:29,0], color=c, linestyle=ls)

	
					ax.fill_between(x, data[0:29,1], data[0:29,2], facecolor=c, edgecolor = None, alpha=0.4)
		
	plt.xlim(0,7)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)

	plt.xlabel('Lead Time (Days)', fontsize = 14)
	if stat_type == "location_error":
		plt.ylabel('Track Error (km)', fontsize=14)
		plt.ylim(bottom=0)
	elif stat_type == "mslp_bias":
		plt.ylabel('MSLP Bias (hPa)', fontsize=14)
		plt.hlines(0, x[0], x[-1], color='silver',linewidth=0.6)
	elif stat_type == "wind_bias":
		plt.ylabel('10m Wind Bias (km/h)', fontsize=14)
		plt.hlines(0, x[0], x[-1], color='silver',linewidth=0.6)
	elif stat_type == "sample_size":
		plt.ylabel('Number of Forecasts',fontsize=14)

	#legend	
	
	ecmwf = plt.Line2D((0, 1), (0, 0), color="blue",linewidth=1, linestyle='-')
	ukmo = plt.Line2D((0, 1), (0, 0), color="#ffa600",linewidth=1, linestyle='-')
	gap = plt.Line2D((0, 1), (0, 0), color='white',linewidth=1, linestyle=ls)
	hres = plt.Line2D((0, 1), (0, 0), color="k",linewidth=1, linestyle='-')
	ctrl = plt.Line2D((0, 1), (0, 0), color="k",linewidth=1, linestyle='-.')
	mean = plt.Line2D((0, 1), (0, 0), color="k",linewidth=1, linestyle='--')
	eps = plt.Line2D((0, 1), (0, 0), color="k",linewidth=1, linestyle=':')
	

	legend = ax.legend((ecmwf, ukmo, gap,hres, ctrl, mean, eps), ['ECMWF', 'UKMO',' ','HRES','Control', 'ENS Mean', 'ENS Members'], fontsize=12)
	
		
	plt.savefig("UKMO_2010-2020.vs.ECMWF_2012-2020_2014missing."+stat_type+".png", dpi=400, bbox_inches='tight')
	
	plt.close()


for obs_track in ["ibtracs"]:
	plot_statistics("location_error", obs_track)
	plot_statistics("mslp_bias", obs_track)
	plot_statistics("wind_bias", obs_track)
	#plot_statistics("sample_size", obs_track)



