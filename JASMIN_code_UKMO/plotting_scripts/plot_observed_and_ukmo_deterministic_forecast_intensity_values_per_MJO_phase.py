import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches


standard_MJO_pairs = [23, 45,  67, 81]
alt_MJO_pairs = [12, 34, 56, 78]

datadirdet = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/MJO_phase/"
datadirens = "/gws/nopw/j04/klingaman/emerton/ukmo_TIGGE_analysis/track_and_intensity_errors/MJO_phase/"
savedir = "./"#"/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/plots/"

if not os.path.exists(savedir):
    os.makedirs(savedir)

def plot_statistics(stat_type, fcst_type, obs_track):

	standard_MJO_pairs = [23, 45,  67, 81]
	alt_MJO_pairs = [12, 34, 56, 78]

	for MJO_pairs in [standard_MJO_pairs, alt_MJO_pairs]: #standard_MJO_pairs, alt_MJO_pairs

	
		fig, ax = plt.subplots()
		fig.set_size_inches(6,6)

		x = np.arange(0,7.25,0.25)
		
		colours = ['#fcc200', '#f05238', '#a1005c', '#08025c']
	
		#if fcst_type == 'hres':
			#ls = '-'
		#elif fcst_type == 'ctrl':
			#ls = '-.'
		#elif fcst_type == 'mean':
			#ls = '--'
		#elif fcst_type == 'EPS':
			#ls = ':'
	
	
		for MJO, c in zip(MJO_pairs, colours): 
		
			print MJO_pairs
			print MJO
			print c
		
	
			if fcst_type == "hres":
	
				fcstdata = np.genfromtxt(datadirdet+"2010_2020_MJO_phase"+str(MJO)+"_average_UKMO_deterministic_"+stat_type+"_values_per_lead_time_vs_"+obs_track+".with_confidence_intervals.txt")
		
				obsdata = np.genfromtxt(datadirdet+"2010_2020_MJO_phase"+str(MJO)+"_average_ibtracs_"+stat_type+"_values_per_lead_time.with_confidence_intervals.txt")
				
			elif fcst_type == "EPS":
			
				data = np.genfromtxt(datadirens+"2010_2020_MJO_phase"+str(MJO)+"_average_"+stat_type+"_per_lead_time_vs_"+obs_track+"_"+fcst_type+".with_confidence_intervals.txt")
			else:
				data = np.genfromtxt(datadirens+"2010_2020_MJO_phase"+str(MJO)+"_average_"+stat_type+"_per_lead_time_vs_"+obs_track+"_"+fcst_type+".with_confidence_intervals_95.txt")
		
			if stat_type == "wind_speed":
				
				plt.plot(x, fcstdata[0:29,0]*3.6, color=c, linestyle='--')
					
				ax.fill_between(x, fcstdata[0:29,1]*3.6, fcstdata[0:29,2]*3.6, facecolor=c, edgecolor = None, alpha=0.2)
				
				plt.plot(x, obsdata[0:29,0]*3.6, color=c)
				
				ax.fill_between(x, obsdata[0:29,1]*3.6, obsdata[0:29,2]*3.6, facecolor=c, edgecolor=None, alpha=0.2)
				
			else:
				plt.plot(x, fcstdata[0:29,0], color=c, linestyle='--')

				ax.fill_between(x, fcstdata[0:29,1], fcstdata[0:29,2], facecolor=c, edgecolor = None, alpha=0.2)
				
				plt.plot(x, obsdata[0:29,0], color=c)
				
				ax.fill_between(x, obsdata[0:29,1], obsdata[0:29,2], facecolor=c, edgecolor=None, alpha=0.2)
		
		plt.xlim(0,7)
		plt.xticks(fontsize=12)
		plt.yticks(fontsize=12)

		plt.xlabel('Lead Time (Days)', fontsize = 14)
		
		if stat_type == "mslp":
			plt.ylabel('MSLP (hPa)', fontsize=14)
			#plt.hlines(0, x[0], x[-1], color='silver',linewidth=0.6)
			
		elif stat_type == "wind_speed":
			plt.ylabel('10m Wind (km/h)', fontsize=14)
			#plt.hlines(0, x[0], x[-1], color='silver',linewidth=0.6)
		
		#legend	
	
	
		#if obs_track == "ibtracs":
			#obs = plt.Line2D((0, 1), (0, 0), color='white',linewidth=1, linestyle=ls, label=fcst_type+" verified against IBTrACS")
			#title = fcst_type+" verified against IBTrACS"
			
		#elif obs_track == "analysis":
			#obs = plt.Line2D((0, 1), (0, 0), color='white',linewidth=1, linestyle=ls, label=fcst_type+" verified against ECMWF analysis")
			#title = fcst_type+" verified against ECMWF analysis"
			
		#title = fcst_type+" "+stat_type
		
		
		
		if MJO_pairs == standard_MJO_pairs:	
			p1 = plt.Line2D((0, 1), (0, 0), color=colours[0], linewidth=1.5, label='MJO 2-3') 
			p2 = plt.Line2D((0, 1), (0, 0), color=colours[1], linewidth=1.5, label='MJO 4-5') 
			p3 = plt.Line2D((0, 1), (0, 0), color=colours[2], linewidth=1.5, label='MJO 6-7') 
			p4 = plt.Line2D((0, 1), (0, 0), color=colours[3], linewidth=1.5, label='MJO 8-1')
			
			sa = 'standard'
		
		elif MJO_pairs == alt_MJO_pairs:	
			p1 = plt.Line2D((0, 1), (0, 0), color=colours[0], linewidth=1.5, label='MJO 1-2') 
			p2 = plt.Line2D((0, 1), (0, 0), color=colours[1], linewidth=1.5, label='MJO 3-4') 
			p3 = plt.Line2D((0, 1), (0, 0), color=colours[2], linewidth=1.5, label='MJO 5-6') 
			p4 = plt.Line2D((0, 1), (0, 0), color=colours[3], linewidth=1.5, label='MJO 7-8')
			sa = 'alternative'
		
		ob = plt.Line2D((0, 1), (0, 0), color='k', linewidth=1.5, label='Observed (IBTrACS)')
		fc = plt.Line2D((0, 1), (0, 0), color='k', linewidth=1.5, linestyle='--',label='Forecast ('+fcst_type+')')
	
		#if stat_type == "location_error":
			
		
		
		if stat_type == 'wind_speed':
			legend = ax.legend(handles=[p1,p2,p3,p4,fc,ob], fontsize=10, loc='upper left')
			legend._legend_box.align = "left"
			
		elif stat_type == 'mslp':
			legend = ax.legend(handles=[p1,p2,p3,p4,fc,ob], fontsize=10, loc='lower left')
			legend._legend_box.align = "left"
			
			
			#if obs_track == 'ibtracs':
		
				#if fcst_type == "hres":
					#legend = ax.legend((obs,res1, res2, res3), [fcst_type+' verified against IBTrACS','N512 / ~40km (July 2010 - June 2014)','N768 / ~26km (July 2014 - June 2017)', 'N1280 / ~16km (July 2017 - June 2019)'], fontsize=12)
				#else:	
					#legend = ax.legend((obs,res1, res2, res3), [fcst_type+' verified against IBTrACS','N216 / ~93km (July 2010 - June 2014)','N400 / ~50km (July 2014 - June 2017)', 'N640 / ~31km (July 2017 - June 2019)'], fontsize=12)
		
		
			#elif obs_track == 'analysis':
				#if fcst_type == "hres":
					#legend = ax.legend((obs,res1, res2, res3), [fcst_type+' verified against ECMWF analysis','N512 / ~40km (July 2010 - June 2014)','N768 / ~26km (July 2014 - June 2017)', 'N1280 / ~16km (July 2017 - June 2019)'], fontsize=12)
				#else:	
					#legend = ax.legend((obs,res1, res2, res3), [fcst_type+' verified against ECMWF analysis','N216 / ~93km (July 2010 - June 2014)','N400 / ~50km (July 2014 - June 2017)', 'N640 / ~31km (July 2017 - June 2019)'], fontsize=12)
				
				
		
		if MJO_pairs == standard_MJO_pairs:
			plt.savefig("ukmo_"+fcst_type+".2010-2020.standard_MJO_pairs."+stat_type+"_values_and_ibtracs_values.with_conf_ints.png", dpi=400, bbox_inches='tight')
		elif MJO_pairs == alt_MJO_pairs:
			plt.savefig("ukmo_"+fcst_type+".2010-2020.alt_MJO_pairs."+stat_type+"_values_and_ibtracs_values.with_conf_ints.png", dpi=400, bbox_inches='tight')
	
		plt.close()

for fcst_type in ["hres"]: #"hres", "ctrl", "mean", "EPS"
	for obs_track in ["ibtracs"]:
		
		plot_statistics("mslp", fcst_type, obs_track)
		plot_statistics("wind_speed", fcst_type, obs_track)
	



