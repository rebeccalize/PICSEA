import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches


#y1=sys.argv[1]
#y2=sys.argv[2]

datadir = "/gws/nopw/j04/klingaman/emerton/PICSEA/IDAI_KENNETH/IDAI_KENNETH_TRACKS_RESULTS/idai_kenneth_ecwmf_ens_det_ctrl_track_statistics/"


#if not os.path.exists(savedir):
    #os.makedirs(savedir)

#storms_dir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/IDAI/reformatted_idai_forecast_tracks/"


def plot_statistics(stat_type, cyclone):
	"""stat_type is either "location_error" or "mlsp_bias" or "wind_bias"""
	fig, ax = plt.subplots()
	fig.set_size_inches(6,6)

	det_data_ib = np.genfromtxt(datadir+cyclone+"_DET_average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")
	ctrl_data_ib = np.genfromtxt(datadir+cyclone+"_CTRL_average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")
	mean_data_ib = np.genfromtxt(datadir+cyclone+"_MEAN_average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")
	eps_data_ib = np.genfromtxt(datadir+cyclone+"_EPS_average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")

	x = np.arange(0,7.25,0.25)
	
	#bluegreen pallette used for poster
	#plt.plot(x, det_data_ib, color='#49ff0d')
	#plt.plot(x, ctrl_data_ib, color='#00b9b6')
	#plt.plot(x, mean_data_ib, color='#0066cb')
	#plt.plot(x, eps_data_ib, color='#12025c')
	
	#matches maps
	plt.plot(x, det_data_ib, color='#F71735')
	#plt.plot(x, ctrl_data_ib, color='#FF9F1C')
	plt.plot(x, mean_data_ib, color='mediumblue')
	plt.plot(x, eps_data_ib, color='lightblue')

	plt.xlim(0,7)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)

	plt.xlabel('Lead Time (Days)', fontsize = 15)
	if stat_type == "location_error":
		plt.ylabel('Track Error (km)', fontsize=15)
		plt.ylim(0,1000)
	elif stat_type == "mslp_bias":
		plt.ylabel('MSLP Bias (hPa)', fontsize=15)
	elif stat_type == "wind_bias":
		plt.ylabel(r'10m Wind Bias (ms$^{-1}$)', fontsize=15)
		
		
	#det = plt.Line2D((0, 1), (0, 0), color='#49ff0d', linewidth=1, label = 'Deterministic')
	#ctrl = plt.Line2D((0, 1), (0, 0), color='#00b9b6', linewidth=1, label = 'Control')
	#mean = plt.Line2D((0, 1), (0, 0), color='#0066cb', linewidth=1, label = 'Ensemble Mean')
	#eps = plt.Line2D((0, 1), (0, 0), color='#12025c', linewidth=1, label = 'Ensemble Members')
	
	det = plt.Line2D((0, 1), (0, 0), color='#F71735', linewidth=1, label = 'Deterministic')
	ctrl = plt.Line2D((0, 1), (0, 0), color='#FF9F1C', linewidth=1, label = 'Control')
	mean = plt.Line2D((0, 1), (0, 0), color='mediumblue', linewidth=1, label = 'Ensemble Mean')
	eps = plt.Line2D((0, 1), (0, 0), color='lightblue', linewidth=1, label = 'Ensemble Members (Mean Error)')

	if cyclone == "idai":
		cyclonecaps = "Idai"
	elif cyclone == "kenneth":
		cyclonecaps = "Kenneth"
		
	legend = ax.legend(handles=[det,mean,eps], fontsize=14, loc='upper left', title = "Cyclone "+cyclonecaps)
	legend._legend_box.align = "left"
	plt.setp(legend.get_title(), fontsize='14')

	plt.savefig(cyclone+"_"+stat_type+"_ECMWF_eps_det_ctrl_mean_mapcolours_noctrl.png", dpi=400, bbox_inches='tight')

plot_statistics("location_error","idai")
#plot_statistics("mslp_bias","idai")
#plot_statistics("wind_bias","idai")

plot_statistics("location_error","kenneth")
#plot_statistics("mslp_bias","kenneth")
#plot_statistics("wind_bias","kenneth")



