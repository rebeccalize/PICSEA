import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches




datadirdet = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"
datadirens = "/gws/nopw/j04/klingaman/emerton/ukmo_TIGGE_analysis/track_and_intensity_errors/"
savedir = "./"#"/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/plots/"

if not os.path.exists(savedir):
    os.makedirs(savedir)

def plot_statistics(stat_type, obs_track):

	
	"""stat_type is either "location_error" or "mslp_bias" or wind_speed_bias """
	fig, ax = plt.subplots()
	fig.set_size_inches(6,6)

	x = np.arange(0,7.25,0.25)
	
	#c="#FCA311"
			
	for fcst_type in ["hres", "ctrl", "mean", "EPS"]:
	
		if fcst_type == 'hres':
			ls = '-'
			c = '#FFA62B'
		elif fcst_type == 'ctrl':
			ls = '-.'
			c='#7E52A0'
		elif fcst_type == 'mean':
			ls = '--'
			c='#0B2027'
		elif fcst_type == 'EPS':
			ls = ':'	
			c='#82C0CC'
		if stat_type == "sample_size":
		
			if fcst_type == "hres":
	
				data = np.genfromtxt(datadirdet+"2010_2020.ukmo_hres.number_of_forecasts_included_vs_"+obs_track+".txt")
				
			elif fcst_type == "EPS":
				data = np.genfromtxt(datadirens+"2010-2020.number_of_forecasts_included_vs_"+obs_track+"_"+fcst_type+".txt")
			
		
			else:
				data = np.genfromtxt(datadirens+"2010-2020number_of_forecasts_included_vs_"+obs_track+"_"+fcst_type+".txt")
			
	
		
			plt.plot(x, data[0:29], color=c, linestyle=ls)
		
		
		else:

			if fcst_type == "hres":
	
				data = np.genfromtxt(datadirdet+"2010_2020.ukmo_hres.average_"+stat_type+"_per_lead_time_vs_"+obs_track+"_with_confidence_intervals_95.txt")
				
			elif fcst_type == "EPS":
			
				data = np.genfromtxt(datadirens+"2010-2020.average_"+stat_type+"_per_lead_time_vs_"+obs_track+"_"+fcst_type+"_with_confidence_intervals_95.txt")
		
			else:
				data = np.genfromtxt(datadirens+"2010-2020average_"+stat_type+"_per_lead_time_vs_"+obs_track+"_"+fcst_type+"_with_confidence_intervals_95.txt")
	
			
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
	
	obs = plt.Line2D((0, 1), (0, 0), color='white',linewidth=1, linestyle=ls)
	hres = plt.Line2D((0, 1), (0, 0), color="#FFA62B",linewidth=1, linestyle='-')
	ctrl = plt.Line2D((0, 1), (0, 0), color="#7E52A0",linewidth=1, linestyle='-.')
	mean = plt.Line2D((0, 1), (0, 0), color="#0B2027",linewidth=1, linestyle='--')
	eps = plt.Line2D((0, 1), (0, 0), color="#82C0CC",linewidth=1, linestyle=':')
	
	if stat_type == "location_error":
		if obs_track == 'ibtracs':
			legend = ax.legend((obs,hres, ctrl, mean, eps), ['Verified against IBTrACS','HRES','Control', 'ENS Mean', 'ENS Members'], fontsize=12)
		elif obs_track == 'analysis':
			legend = ax.legend((obs,hres, ctrl, mean,eps), ['Verified against ECMWF analysis','HRES','Control', 'ENS Mean', 'ENS Members'], fontsize=12)
		
	plt.savefig("ukmo_all_models.2010-2020."+stat_type+".vs_"+obs_track+".with_conf_ints.png", dpi=400, bbox_inches='tight')
	
	plt.close()


for obs_track in ["ibtracs", "analysis"]:
	plot_statistics("location_error", obs_track)
	plot_statistics("mslp_bias", obs_track)
	plot_statistics("wind_bias", obs_track)
	plot_statistics("sample_size", obs_track)



