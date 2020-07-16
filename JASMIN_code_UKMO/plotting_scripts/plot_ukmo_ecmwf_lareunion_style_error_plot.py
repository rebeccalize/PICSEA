import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches


y1s=[2010,2011,2012,2013,2014,2015,2016,2017,2018,2019]
y2s=[2011,2012,2013,2014,2015,2016,2017,2018,2019,2020]

lead_time_indices = [0, 2, 4, 8, 12] #, 20    #0 hours, 12 hours, 1 day, 2 days, 3 days, 5 days ahead
colours = ['k', '#390099', '#FF0054', '#FF5400', '#FFBD00'] #, 'skyblue'

LR = np.zeros((6,10))

LR[0,:] = [33, 29, 28, 25, 25, 27, 29, 30, np.nan, np.nan]
LR[1,:] = [79, 68, 58, 55, 66, 65, 70, 62, np.nan, np.nan]
LR[2,:] = [107, 106, 82, 83, 105, 97, 98, 84, np.nan, np.nan]
LR[3, :] = [174, 198, 143, 149, 184, 159, 153, 121, np.nan, np.nan]
LR[4, :] = [254, 290, 212, 244, 272, 231, 210, 143, np.nan, np.nan]
LR[5, :] = np.nan

savedir = "./"#"/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/plots/"


def plot_statistics(stat_type, fcst_type):


	fig, ax = plt.subplots()
	fig.set_size_inches(10,6)
	
	for lt,c, z in zip(lead_time_indices, colours, range(6)):
		for centre in ['ukmo','ecmwf', 'reunion']:
		
			if centre == 'ukmo':
				if fcst_type == 'det':	
					datadir = '/gws/nopw/j04/klingaman/emerton/UKMO_HRES_RESULTS/track_and_intensity_errors/'
					fcst_label = 'hres'
				else:
					datadir = '/gws/nopw/j04/klingaman/emerton/UKMO_ENS_RESULTS/track_and_intensity_errors/'
					fcst_label = fcst_type
				ls = '--'
				
			elif centre == 'ecmwf':
				datadir = '/gws/nopw/j04/klingaman/emerton/ECMWF_RESULTS/copy_fri_3_july_not_final/'
				ls = '-'
				fcst_label = fcst_type
				
		
			mean_data_to_plot = np.zeros(10)
			conf_int_min = np.zeros(10)
			conf_int_max = np.zeros(10)
	
			if centre == 'reunion':
				mean_data_to_plot = LR[z, :]
				conf_int_min[:] = np.nan
				conf_int_max[:] = np.nan
				ls=':'
				
			else:
				for y1,y2,i in zip(y1s,y2s,range(10)):
			
					if centre == 'ecmwf':
						if y1 == 2010 or y1 == 2011:
					
							mean_data_to_plot[i] = np.nan
							conf_int_min[i] = np.nan
							conf_int_max[i] = np.nan
						
						else:
							data = np.genfromtxt(datadir+str(y1)+"-"+str(y2)+"average_location_error_per_lead_time_vs_ibtracs_"+fcst_label+"_with_confidence_intervals_95.txt")
							mean_data_to_plot[i] = data[lt, 0]
							conf_int_min[i] = data[lt, 1]
							conf_int_max[i] = data[lt, 2]
						
						
					elif centre == 'ukmo':
				
						if fcst_type == 'det':
					
							data = np.genfromtxt(datadir+str(y1)+"_"+str(y2)+".ukmo_hres.average_location_error_per_lead_time_vs_ibtracs_with_confidence_intervals_95.txt")
							mean_data_to_plot[i] = data[lt,0]
							conf_int_min[i] = data[lt,1]
							conf_int_max[i] = data[lt,2]
						
						else:
							data = np.genfromtxt(datadir+str(y1)+"-"+str(y2)+"average_location_error_per_lead_time_vs_ibtracs_"+fcst_label+"_with_confidence_intervals_95.txt")
							mean_data_to_plot[i] = data[lt,0]
							conf_int_min[i] = data[lt,1]
							conf_int_max[i] = data[lt,2]
						
				
						
						
			
			x = np.arange(0,10,1)
			plt.plot(x, mean_data_to_plot, color=c,linestyle=ls)
			
			plt.scatter(x, mean_data_to_plot, color=c, s=15)
			ax.fill_between(x, conf_int_min, conf_int_max, facecolor=c, edgecolor=None, alpha=0.1)
			
			for i in range(10):
				if not np.isnan(mean_data_to_plot[i]):
					plt.text(i, mean_data_to_plot[i]+4, str(int(round(mean_data_to_plot[i]))), fontsize=5.5)
						



	

	my_xticks = ['2010-2011 (5)', '2011-2012 (11)', '2012-2013 (10)', '2013-2014 (11)', '2014-2015 (11)', '2015-2016 (8)', '2016-2017 (5)', '2017-2018 (8)', '2018-2019 (14)', '2019-2020 (11)']
	plt.yticks(fontsize=10)
	plt.xticks(x, my_xticks,fontsize=10,rotation=90)
	plt.xlabel('Cyclone Season', fontsize = 12)
	if stat_type == "location_error":
		plt.ylabel('Track Error (km)', fontsize=12)
		plt.ylim(bottom=0)
	elif stat_type == "mslp_bias":
		plt.ylabel('MSLP Bias (hPa)', fontsize=12)
	elif stat_type == "wind_bias":
		plt.ylabel(r'10m Wind Bias (ms$^{-1}$)', fontsize=12)
		
	
		
	p1 = plt.Line2D((0, 1), (0, 0), color=colours[0], linewidth=1.5, label='H+0') 
	p2 = plt.Line2D((0, 1), (0, 0), color=colours[1], linewidth=1.5, label='H+12') 
	p3 = plt.Line2D((0, 1), (0, 0), color=colours[2], linewidth=1.5, label='D+1') 
	p4 = plt.Line2D((0, 1), (0, 0), color=colours[3], linewidth=1.5, label='D+2')
	p5 = plt.Line2D((0, 1), (0, 0), color=colours[4], linewidth=1.5, label='D+3')
	#p6 = plt.Line2D((0, 1), (0, 0), color=colours[5], linewidth=1.5, label='D+5')
	p7 = plt.Line2D((0, 1), (0, 0), color='white', linewidth=1.5, label=' ')
	p8 = plt.Line2D((0, 1), (0, 0), color='grey', linewidth=1.5, label='ECMWF')
	p9 = plt.Line2D((0, 1), (0, 0), color='grey', linewidth=1.5, linestyle='--',label='UKMO')
	p10 = plt.Line2D((0, 1), (0, 0), color='grey', linewidth=1.5, linestyle=':',label='La Reunion')
	
	legend = ax.legend(handles=[p1,p8, p2,p9,p3,p10, p4, p5], fontsize=10, loc='upper center',bbox_to_anchor=(0.5, 1.05), ncol=5, facecolor='white', framealpha=1)
		
	axes = plt.gca()
	axes.yaxis.grid()

	plt.savefig("ukmo_ecmwf_reunion_"+fcst_type+"_2010-2020_"+str(stat_type)+"_annual_evolution_lareunion_style_plot_multi_lead_time.png", bbox_inches='tight', dpi = 400)
	plt.close()

plot_statistics("location_error","det")
#plot_statistics("mslp_bias")
#plot_statistics("wind_bias")



