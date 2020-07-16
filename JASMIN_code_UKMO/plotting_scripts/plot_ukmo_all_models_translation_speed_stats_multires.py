import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches




datadirdet = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/translation_speed_errors/"
datadirens = "/gws/nopw/j04/klingaman/emerton/ukmo_TIGGE_analysis/translation_speed_errors/"
savedir = "./"#"/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/plots/"

if not os.path.exists(savedir):
    os.makedirs(savedir)

def plot_statistics(stat_type, fcst_type, obs_track):

	
	"""converts speed errors from m/s to km/h during plotting by *3.6"""
	
	fig, ax = plt.subplots()
	fig.set_size_inches(6,6)

	x = np.arange(0,7.25,0.25)
	
	res_start_years = [2010, 2014, 2017]
	res_end_years = [2014, 2017, 2020]
	
	if fcst_type == 'hres':
		ls = '-'
	elif fcst_type == 'ctrl':
		ls = '-.'
	elif fcst_type == 'mean':
		ls = '--'
	elif fcst_type == 'EPS':
		ls = ':'
	
	for y1, y2 in zip(res_start_years, res_end_years):
	
		if y1 == 2010:
			c = '#65CBE9' #'#6996b3'
			
		elif y1 == 2014:
			c = '#004c6d'
			
		elif y1 == 2017:
			c = 'red'
			
		
	
		if fcst_type == "hres":
	
			data = np.genfromtxt(datadirdet+str(y1)+"-"+str(y2)+".average_"+stat_type+"_per_lead_time.vs_"+obs_track+".with_confidence_intervals.txt")
		
		elif fcst_type == "EPS":
			data = np.genfromtxt(datadirens+str(y1)+"-"+str(y2)+".average_"+stat_type+"_per_lead_time_vs_"+obs_track+"_"+fcst_type+"_with_confidence_intervals_95.txt")
		else:
			data = np.genfromtxt(datadirens+str(y1)+"-"+str(y2)+".average_"+stat_type+"_per_lead_time.vs_"+obs_track+"."+fcst_type+".with_confidence_intervals.txt")
	
		
		
		plt.plot(x, data[0:29,0]*3.6, color=c, linestyle=ls)

		ax.fill_between(x, data[0:29,1]*3.6, data[0:29,2]*3.6, facecolor=c, edgecolor = None, alpha=0.4)
		
	plt.xlim(0,7)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)

	plt.xlabel('Lead Time (Days)', fontsize = 14)
	if stat_type == "translation_speed_error":
		plt.ylabel('Translation Speed Error (km/h)', fontsize=14)
		plt.ylim(bottom=0)
	elif stat_type == "translation_speed_bias":
		plt.ylabel('Translation Speed Bias (km/h)', fontsize=14)
		plt.hlines(0, x[0], x[-1], color='silver',linewidth=0.6)


	#legend	
	
	obs = plt.Line2D((0, 1), (0, 0), color='white',linewidth=1, linestyle=ls)
	res1 = plt.Line2D((0, 1), (0, 0), color='#65CBE9',linewidth=1, linestyle=ls)
	res2 = plt.Line2D((0, 1), (0, 0), color='#004c6d',linewidth=1, linestyle=ls)
	res3 = plt.Line2D((0, 1), (0, 0), color='red',linewidth=1, linestyle=ls)
	
	if stat_type == "location_error":
		if obs_track == 'ibtracs':
		
			if fcst_type == "hres":
				legend = ax.legend((obs,res1, res2, res3), [fcst_type+' verified against IBTrACS','N512 / ~40km (July 2010 - June 2014)','N768 / ~26km (July 2014 - June 2017)', 'N1280 / ~16km (July 2017 - June 2020)'], fontsize=12)
			else:
				legend = ax.legend((obs,res1, res2, res3), [fcst_type+' verified against IBTrACS','N216 / ~93km (July 2010 - June 2014)','N400 / ~50km (July 2014 - June 2017)', 'N640 / ~31km (July 2017 - June 2020)'], fontsize=12)
		
		
		elif obs_track == 'analysis':
			if fcst_type == "hres":
				legend = ax.legend((obs,res1, res2, res3), [fcst_type+' verified against ECMWF analysis','N512 / ~40km (July 2010 - June 2014)','N768 / ~26km (July 2014 - June 2017)', 'N1280 / ~16km (July 2017 - June 2020)'], fontsize=12)
			else:	
				legend = ax.legend((obs,res1, res2, res3), [fcst_type+' verified against ECMWF analysis','N216 / ~93km (July 2010 - June 2014)','N400 / ~50km (July 2014 - June 2017)', 'N640 / ~31km (July 2017 - June 2020)'], fontsize=12)
				
				
	plt.savefig("ukmo_"+fcst_type+".multi_res."+stat_type+".vs_"+obs_track+".with_conf_ints.2010-2020.png", dpi=400, bbox_inches='tight')
	plt.close()

for fcst_type in ["hres"]: # "hres", "ctrl", "mean", "EPS"
	for obs_track in ["ibtracs", "analysis"]:
		plot_statistics("translation_speed_bias", fcst_type, obs_track)
		#plot_statistics("translation_speed_error", fcst_type, obs_track)
		



