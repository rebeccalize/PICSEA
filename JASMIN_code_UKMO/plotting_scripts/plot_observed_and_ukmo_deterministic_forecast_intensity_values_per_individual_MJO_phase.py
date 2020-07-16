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

	ind_MJO_phases=[1,2,3,4,5,6,7,8]

	for MJO_pairs in [ind_MJO_phases]: #standard_MJO_pairs, alt_MJO_pairs

	
		
		
		colours = ['#EF476F', '#F78C6B', '#FFD166', '#06D6A0','#34E4EA','#118AB2','#073B4C','#9B7EDE']
	
		#if fcst_type == 'hres':
			#ls = '-'
		#elif fcst_type == 'ctrl':
			#ls = '-.'
		#elif fcst_type == 'mean':
			#ls = '--'
		#elif fcst_type == 'EPS':
			#ls = ':'
	
		p1 = plt.Line2D((0, 1), (0, 0), color=colours[0], linewidth=1.5, label='MJO 1') 
		p2 = plt.Line2D((0, 1), (0, 0), color=colours[1], linewidth=1.5, label='MJO 2') 
		p3 = plt.Line2D((0, 1), (0, 0), color=colours[2], linewidth=1.5, label='MJO 3') 
		p4 = plt.Line2D((0, 1), (0, 0), color=colours[3], linewidth=1.5, label='MJO 4')
		p5 = plt.Line2D((0, 1), (0, 0), color=colours[4], linewidth=1.5, label='MJO 5')
		p6 = plt.Line2D((0, 1), (0, 0), color=colours[5], linewidth=1.5, label='MJO 6')
		p7 = plt.Line2D((0, 1), (0, 0), color=colours[6], linewidth=1.5, label='MJO 7')
		p8 = plt.Line2D((0, 1), (0, 0), color=colours[7], linewidth=1.5, label='MJO 8')
		
		legend_lines=[p1,p2,p3,p4,p5,p6,p7,p8]
	
		for MJO, c, p in zip(MJO_pairs, colours, legend_lines): 
		
			print MJO_pairs
			print MJO
			print c
			
			fig, ax = plt.subplots()
			fig.set_size_inches(6,6)

			x = np.arange(0,7.25,0.25)
		
	
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
				plt.ylim(970,1005)
				plt.ylabel('MSLP (hPa)', fontsize=14)
				#plt.hlines(0, x[0], x[-1], color='silver',linewidth=0.6)
			
			elif stat_type == "wind_speed":
				plt.ylim(50,150)
				plt.ylabel('10m Wind (km/h)', fontsize=14)
				#plt.hlines(0, x[0], x[-1], color='silver',linewidth=0.6)
		
	
		
		
			
			ob = plt.Line2D((0, 1), (0, 0), color='k', linewidth=1.5, label='Observed (IBTrACS)')
			fc = plt.Line2D((0, 1), (0, 0), color='k', linewidth=1.5, linestyle='--',label='Forecast ('+fcst_type+')')
	
			
		
		
			if stat_type == 'wind_speed':
				legend = ax.legend(handles=[p,fc,ob], fontsize=10, loc='upper left')
				legend._legend_box.align = "left"
			
			elif stat_type == 'mslp':
				legend = ax.legend(handles=[p,fc,ob], fontsize=10, loc='lower left')
				legend._legend_box.align = "left"
			
		
		
		
			plt.savefig("ukmo_"+fcst_type+".2010-2020.MJO_phase_"+str(MJO)+"."+stat_type+"_values_and_ibtracs_values.with_conf_ints.png", dpi=400, bbox_inches='tight')
		
	
			plt.close()

for fcst_type in ["hres"]: #"hres", "ctrl", "mean", "EPS"
	for obs_track in ["ibtracs"]:
		
		#plot_statistics("mslp", fcst_type, obs_track)
		plot_statistics("wind_speed", fcst_type, obs_track)
	



