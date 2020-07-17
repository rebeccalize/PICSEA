import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches


y1=sys.argv[1]
y2=sys.argv[2]

datadir = "/perm/mo/more/picsea/track_and_intensity_errors/Y"+str(y1)+str(y2)+"/SIO_storms/"

storms_dir = "/perm/mo/more/TIGGE/Y"+str(y1)+str(y2)+"/reformatted_track_files_per_storm_without_mslp_uv10/SIO_storms/"
storm_dirs=[]
for root,dirs,files in os.walk(storms_dir):
	for dir in dirs:
		storm_dirs.append(dir)
NS = len(storm_dirs) #number of storms


#if not os.path.exists(savedir):
    #os.makedirs(savedir)

#storms_dir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/IDAI/reformatted_idai_forecast_tracks/"


def plot_statistics(stat_type, y1, y2):
	"""stat_type is either "location_error" or "mlsp_bias" or "wind_bias"""
	fig, ax = plt.subplots()
	fig.set_size_inches(6,6)

	det_data_ib = np.genfromtxt(datadir+"average_"+str(stat_type)+"_per_lead_time_vs_ibtracs_det.txt")
	ctrl_data_ib = np.genfromtxt(datadir+"average_"+str(stat_type)+"_per_lead_time_vs_ibtracs_ctrl.txt")
	mean_data_ib = np.genfromtxt(datadir+"average_"+str(stat_type)+"_per_lead_time_vs_ibtracs_mean.txt")
	eps_data_ib = np.genfromtxt(datadir+"average_"+str(stat_type)+"_per_lead_time_vs_ibtracs_EPS.txt")
	
	det_data_an = np.genfromtxt(datadir+"average_"+str(stat_type)+"_per_lead_time_vs_analysis_det.txt")
	ctrl_data_an = np.genfromtxt(datadir+"average_"+str(stat_type)+"_per_lead_time_vs_analysis_ctrl.txt")
	mean_data_an = np.genfromtxt(datadir+"average_"+str(stat_type)+"_per_lead_time_vs_analysis_mean.txt")
	eps_data_an = np.genfromtxt(datadir+"average_"+str(stat_type)+"_per_lead_time_vs_analysis_EPS.txt")

	x = np.arange(0,7.25,0.25)
	
	#bluegreen pallette used for poster
	#plt.plot(x, det_data_ib, color='#49ff0d')
	#plt.plot(x, ctrl_data_ib, color='#00b9b6')
	#plt.plot(x, mean_data_ib, color='#0066cb')
	#plt.plot(x, eps_data_ib, color='#12025c')
	
	#matches maps
	plt.plot(x, det_data_ib[:29], color='#F71735', linestyle='--')
	plt.plot(x, ctrl_data_ib[:29], color='#FF9F1C', linestyle='--')
	plt.plot(x, mean_data_ib[:29], color='darkturquoise', linestyle='--')
	plt.plot(x, eps_data_ib[:29], color='lightblue', linestyle='--')
	
	plt.plot(x, det_data_an[:29], color='#F71735')
	plt.plot(x, ctrl_data_an[:29], color='#FF9F1C')
	plt.plot(x, mean_data_an[:29], color='darkturquoise')
	plt.plot(x, eps_data_an[:29], color='lightblue')
	
	#equidistant colour palette:
	#plt.plot(x, det_data_ib[:29], color='#f2a200', linestyle='--')
	#plt.plot(x, ctrl_data_ib[:29], color='#e83c3a', linestyle='--')
	#plt.plot(x, mean_data_ib[:29], color='#a60062', linestyle='--')
	#plt.plot(x, eps_data_ib[:29], color='#180173', linestyle='--')
	
	#plt.plot(x, det_data_an[:29], color='#f2a200')
	#plt.plot(x, ctrl_data_an[:29], color='#e83c3a')
	#plt.plot(x, mean_data_an[:29], color='#a60062')
	#plt.plot(x, eps_data_an[:29], color='#180173')

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
	mean = plt.Line2D((0, 1), (0, 0), color='darkturquoise', linewidth=1, label = 'Ensemble Mean')
	eps = plt.Line2D((0, 1), (0, 0), color='lightblue', linewidth=1, label = 'Ensemble Members')


	#det = plt.Line2D((0, 1), (0, 0), color='#f2a200', linewidth=1, label = 'Deterministic')
	#ctrl = plt.Line2D((0, 1), (0, 0), color='#e83c3a', linewidth=1, label = 'Control')
	#mean = plt.Line2D((0, 1), (0, 0), color='#a60062', linewidth=1, label = 'Ensemble Mean')
	#eps = plt.Line2D((0, 1), (0, 0), color='#180173', linewidth=1, label = 'Ensemble Members')
	
		
	legend = ax.legend(handles=[det,ctrl,mean,eps], fontsize=12, loc='upper left', title = "SWIO, "+str(y1)+"-"+str(y2)+", "+str(NS)+" Cyclones")
	legend._legend_box.align = "left"
	plt.setp(legend.get_title(), fontsize='12')

	plt.savefig(str(y1)+"_"+str(y2)+"_"+stat_type+"_ECMWF_eps_det_ctrl_mean.png", dpi=400, bbox_inches='tight')

#plot_statistics("location_error","idai")
#plot_statistics("mslp_bias","idai")
#plot_statistics("wind_bias","idai")

#plot_statistics("location_error","kenneth")
#plot_statistics("mslp_bias","kenneth")
#plot_statistics("wind_bias","kenneth")#

plot_statistics("location_error",y1, y2)
#plot_statistics("mslp_bias",y1, y2)
#plot_statistics("wind_bias",y1, y2)



