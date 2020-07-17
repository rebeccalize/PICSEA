import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches



#year1=[2006,2007, 2008, 2009,2010, 2011, 2012, 2013,2014, 2015, 2016, 2017]
#year2=[2007,2008, 2009, 2010,2011, 2012, 2013, 2014,2015, 2016, 2017, 2018]

year1=[2014,2015,2016,2017]
year2=[2015,2016,2017,2018]



datadirloc = "/perm/mo/more/picsea/track_location_errors/"
datadirint = "/perm/mo/more/picsea/intensity_errors/"
#savedir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"

#if not os.path.exists(savedir):
    #os.makedirs(savedir)

#FOR THE 2009-2010 SEASON:
#The resolution was updated in March 2010 (check exact date),
#Will need to add something here to check, for each storm in the 2009-2010 season, whether forecasts were started
#before or after the resolution change (within each tr directory, the filenames have the date in them, use this?)


storm_dirs=[]
for y1,y2 in zip(year1,year2):
	storms_dir = "/perm/mo/more/TIGGE/Y"+str(y1)+str(y2)+"/reformatted_track_files_per_storm_without_mslp_uv10/SIO_storms/"
	for root,dirs,files in os.walk(storms_dir):
		for dir in dirs:
			storm_dirs.append(dir)
print "storm_dirs_loc: ", storm_dirs
print len(storm_dirs)
NSloc = len(storm_dirs) #number of storms


storm_dirs_int=[]
for y1,y2 in zip(year1,year2):
	storms_dir_INT = "/perm/mo/more/TIGGE/Y"+str(y1)+str(y2)+"/reformatted_track_files_per_storm_WITH_INTENSITY/SIO_storms/"
	for root,dirs,files in os.walk(storms_dir_INT):
		for dir in dirs:
			storm_dirs_int.append(dir)
print "storm_dirs_int: ", storm_dirs_int
print len(storm_dirs_int)
NSint = len(storm_dirs_int) #number of storms



def plot_multiyear_stats(stat_type, year1array, year2array, NS):

	if stat_type == "location_error":
		datadir = datadirloc
	else:
		datadir = datadirint
	
	fig, ax = plt.subplots()
	fig.set_size_inches(6,6)

	all_years_ib_error_stats_ec_det=np.genfromtxt(datadir+"Y"+str(year1array[0])+str(year2array[0])+"/SIO_storms/average_" + str(stat_type) + "_per_lead_time_vs_ibtracs_det.txt")
	all_years_ib_error_stats_ec_ctrl=np.genfromtxt(datadir+"Y"+str(year1array[0])+str(year2array[0])+"/SIO_storms/average_" + str(stat_type) + "_per_lead_time_vs_ibtracs_ctrl.txt")
	all_years_ib_error_stats_ec_mean=np.genfromtxt(datadir+"Y"+str(year1array[0])+str(year2array[0])+"/SIO_storms/average_" + str(stat_type) + "_per_lead_time_vs_ibtracs_mean.txt")
	all_years_ib_error_stats_ec_eps=np.genfromtxt(datadir+"Y"+str(year1array[0])+str(year2array[0])+"/SIO_storms/average_" + str(stat_type) + "_per_lead_time_vs_ibtracs_EPS.txt")
	all_years_ib_error_stats_ec_eps = all_years_ib_error_stats_ec_eps[:29]
	ukmo_ib_error_stats_2014_2017=np.genfromtxt("UKMO_DET_2014_2018_average_"+str(stat_type)+"_per_lead_time_vs_ibtracs.txt")
	
	
	#all_years_an_error_stats=np.genfromtxt(datadir+str(year1array[0])+"_"+str(year2array[0])+"/SIO_storms/v2/average_" + str(stat_type) + "_per_lead_time_vs_analysis.txt")
	
	
	
	for y1, y2, in zip(year1array[1:], year2array[1:]):
		print y1,y2
		#year_data_analysis = np.genfromtxt(datadir+str(y1)+"_"+str(y2)+"/SIO_storms/v2/average_"+str(stat_type)+"_per_lead_time_vs_analysis.txt")
		#all_years_an_error_stats = np.vstack([all_years_an_error_stats, year_data_analysis])
		
		
		year_data_ib_ec_det = np.genfromtxt(datadir+"Y"+str(y1)+str(y2)+"/SIO_storms/average_" + str(stat_type) + "_per_lead_time_vs_ibtracs_det.txt")
		all_years_ib_error_stats_ec_det = np.vstack([all_years_ib_error_stats_ec_det, year_data_ib_ec_det[:42]])
		
		year_data_ib_ec_ctrl = np.genfromtxt(datadir+"Y"+str(y1)+str(y2)+"/SIO_storms/average_" + str(stat_type) + "_per_lead_time_vs_ibtracs_ctrl.txt")
		all_years_ib_error_stats_ec_ctrl = np.vstack([all_years_ib_error_stats_ec_ctrl, year_data_ib_ec_ctrl[:42]])
		
		year_data_ib_ec_mean = np.genfromtxt(datadir+"Y"+str(y1)+str(y2)+"/SIO_storms/average_" + str(stat_type) + "_per_lead_time_vs_ibtracs_mean.txt")
		all_years_ib_error_stats_ec_mean = np.vstack([all_years_ib_error_stats_ec_mean, year_data_ib_ec_mean[:42]])
		
		year_data_ib_ec_eps = np.genfromtxt(datadir+"Y"+str(y1)+str(y2)+"/SIO_storms/average_" + str(stat_type) + "_per_lead_time_vs_ibtracs_EPS.txt")
		print all_years_ib_error_stats_ec_eps
		print year_data_ib_ec_eps[:29]
		print '###################################'
		all_years_ib_error_stats_ec_eps = np.vstack([all_years_ib_error_stats_ec_eps, year_data_ib_ec_eps[:29]])


	mean_all_years_ib_error_stats_ec_det, mean_all_years_ib_error_stats_ec_ctrl, mean_all_years_ib_error_stats_ec_mean, mean_all_years_ib_error_stats_ec_eps = (np.zeros(29) for i in range(4))
	for lt in range(29):
		mean_all_years_ib_error_stats_ec_det[lt] = np.nanmean(all_years_ib_error_stats_ec_det[:,lt])
		mean_all_years_ib_error_stats_ec_ctrl[lt] = np.nanmean(all_years_ib_error_stats_ec_ctrl[:,lt])
		mean_all_years_ib_error_stats_ec_mean[lt] = np.nanmean(all_years_ib_error_stats_ec_mean[:,lt])
		mean_all_years_ib_error_stats_ec_eps[lt] = np.nanmean(all_years_ib_error_stats_ec_eps[:,lt])
		
		
	np.savetxt(str(year1array[0])+"_"+str(year2array[-1])+"_average_"+str(stat_type)+"_per_lead_time_vs_ibtracs_ECMWF_det.txt", mean_all_years_ib_error_stats_ec_det[:], '%.4f')
	np.savetxt(str(year1array[0])+"_"+str(year2array[-1])+"_average_"+str(stat_type)+"_per_lead_time_vs_ibtracs_ECMWF_ctrl.txt", mean_all_years_ib_error_stats_ec_ctrl[:], '%.4f')
	np.savetxt(str(year1array[0])+"_"+str(year2array[-1])+"_average_"+str(stat_type)+"_per_lead_time_vs_ibtracs_ECMWF_mean.txt", mean_all_years_ib_error_stats_ec_mean[:], '%.4f')
	np.savetxt(str(year1array[0])+"_"+str(year2array[-1])+"_average_"+str(stat_type)+"_per_lead_time_vs_ibtracs_ECMWF_EPS.txt", mean_all_years_ib_error_stats_ec_eps[:], '%.4f')
	
	
	
	x = np.arange(0,7.25,0.25)
	
	#ECMWF stats
	#plt.plot(x, mean_all_years_ib_error_stats_ec_det[:29], color='#F71735', linestyle='--', linewidth=1.5)
	#plt.plot(x, mean_all_years_ib_error_stats_ec_ctrl[:29], color='#FF9F1C', linestyle='--', linewidth=1.5)
	#plt.plot(x, mean_all_years_ib_error_stats_ec_mean[:29], color='darkturquoise', linestyle='--', linewidth=1.5)
	#plt.plot(x, mean_all_years_ib_error_stats_ec_eps[:29], color='lightblue', linestyle='--', linewidth=1.5)
	
	
	plt.plot(x, mean_all_years_ib_error_stats_ec_det[:29], color='#003f5c', linewidth=1.5)
	plt.plot(x, mean_all_years_ib_error_stats_ec_ctrl[:29], color='#003f5c', linestyle='dashdot', linewidth=1.5)
	plt.plot(x, mean_all_years_ib_error_stats_ec_mean[:29], color='#003f5c', linestyle='--', linewidth=1.5)
	plt.plot(x, mean_all_years_ib_error_stats_ec_eps[:29], color='#003f5c', linestyle='dotted', linewidth=1.5)
	
	
	#UKMO stats	
	#plt.plot(x, ukmo_ib_error_stats_2014_2017[:29], color='#58508d',linestyle='--', linewidth=1.5)
	plt.plot(x, ukmo_ib_error_stats_2014_2017[:29], color='#ffa600', linewidth=1.5)

	plt.xlim(0,7)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)

	plt.xlabel('Lead Time (Days)', fontsize = 14)
	if stat_type == "location_error":
		plt.ylabel('Track Error (km)', fontsize=14)
		plt.ylim(0,1000)
	elif stat_type == "mslp_bias":
		plt.ylabel('MSLP Bias (hPa)', fontsize=14)
		#plt.ylim(-25,70)
		plt.hlines(0,0,7,colors='gainsboro',linewidth=0.75)
	elif stat_type == "wind_bias":
		plt.ylabel(r'10m Wind Bias (ms$^{-1}$)', fontsize=14)
		plt.ylim(-12,2)
		plt.hlines(0,0,7,colors='gainsboro',linewidth=0.75)

	an = plt.Line2D((0, 1), (0, 0), color='k', linewidth=0.5, label = 'Analysis')
	ib = plt.Line2D((0, 1), (0, 0), color='k', linestyle='--', linewidth=0.5, label = 'IBTrACS')
	gap = plt.Line2D((0, 1), (0, 0), color='white', linewidth=0.5, label = ' ')
	
	#ec_det = plt.Line2D((0, 1), (0, 0), color='#F71735', linewidth=1.5, label = 'EC HRES')
	#ec_ctrl = plt.Line2D((0, 1), (0, 0), color='#FF9F1C', linewidth=1.5, label = 'EC Control')
	#ec_mean = plt.Line2D((0, 1), (0, 0), color='darkturquoise', linewidth=1.5, label = 'EC ENS Mean')
	#ec_eps = plt.Line2D((0, 1), (0, 0), color='lightblue', linewidth=1.5, label = 'EC ENS Members')
	
	ecmwf = plt.Line2D((0, 1), (0, 0), color='#003f5c', linewidth=1.5, label = 'ECMWF ('+str(NS)+' TCs)')
	ukmo = plt.Line2D((0, 1), (0, 0), color='#ffa600', linewidth=1.5, label = 'UKMO (32 TCs)')
	
	
	det = plt.Line2D((0, 1), (0, 0), color='darkgrey', linewidth=1.5, label = 'HRES')
	ctrl = plt.Line2D((0, 1), (0, 0), color='darkgrey', linewidth=1.5, linestyle='dashdot', label = 'Control')
	mean = plt.Line2D((0, 1), (0, 0), color='darkgrey', linewidth=1.5, linestyle='--', label = 'ENS Mean')
	eps = plt.Line2D((0, 1), (0, 0), color='darkgrey', linewidth=1.5, linestyle='dotted', label = 'ENS Members')
	
	#uk_det = plt.Line2D((0, 1), (0, 0), color='#58508d', linewidth=1.5, label = 'UKMO HRES (32 TCs)')
	#uk_det = plt.Line2D((0, 1), (0, 0), color='#ffa600', linewidth=1.5, label = 'UKMO HRES (32 TCs)')
	
	
	title=str(year1array[0]) + "-" + str(year2array[-1]) + ", IBTrACS"

	if stat_type == "location_error":
		legend = ax.legend(handles=[ecmwf,ukmo,gap,det,ctrl,mean,eps], fontsize=12, loc='upper left', title = title)
		legend._legend_box.align = "left"
		plt.setp(legend.get_title(), fontsize='12')

	plt.savefig(str(year1array[0])+"_"+str(year2array[-1])+"_"+stat_type+"_ECMWF_and_UKMOdet.png", bbox_inches='tight', dpi=400)




plot_multiyear_stats("location_error", year1,year2, NSloc)
plot_multiyear_stats("mslp_bias", year1, year2, NSint)
plot_multiyear_stats("wind_bias", year1, year2, NSint)



