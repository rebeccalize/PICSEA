import sys
import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import matplotlib.patches as mpatches


y1s=[2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018]
y2s=[2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019]


savedir = "./"#"/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/plots/"

lead_time_days=3
lead_time_index=lead_time_days*4

#c = 'k' #boxplot colour 

def plot_statistics(stat_type):

	if stat_type == "wind_bias":
		c = 'b'
	elif stat_type == "mslp_bias":
		c='g'
	elif stat_type == "location_error":
		c='k'

	boxdata20062007 = []
	boxdata20072008 = []
	boxdata20082009 = []
	boxdata20092010 = []
	boxdata20102011 = []
	boxdata20112012 = []
	boxdata20122013 = []
	boxdata20132014 = []
	boxdata20142015 = []
	boxdata20152016 = []
	boxdata20152016 = []
	boxdata20162017 = []
	boxdata20172018 = []
	boxadta20182019 = []

	fig, ax = plt.subplots()
	fig.set_size_inches(12,6)
	
	
	for y1,y2 in zip(y1s,y2s):
	
		box_data=[]

		datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/"+str(y1)+"_"+str(y2)+"/SIO_storms/"


		storms_dir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/"
		storm_dirs=[]
		for root,dirs,files in os.walk(storms_dir):
			for dir in dirs:
				storm_dirs.append(dir)
		NS = len(storm_dirs) #number of storms
	

		#year_data_analysis = np.genfromtxt(datadir+"average_"+str(stat_type)+"_per_lead_time_vs_analysis.txt")
		#year_data_ibtracs = np.genfromtxt(datadir+"average_" + str(stat_type) + "_per_lead_time_vs_ibtracs.txt")

		#x = np.arange(0,7.25,0.25)
		#plt.plot(x, year_data_analysis, color='k')
		#plt.plot(x, year_data_ibtracs[0:29], linestyle='--')

		list_of_files = os.listdir(datadir)
		ib_pattern = "tr*_average_"+stat_type+"*_ibtracs.txt"
		an_pattern="tr*_average_"+stat_type+"*_analysis.txt"
		ib_files_to_plot=[]
		an_files_to_plot=[]
		for entry in list_of_files:
			if fnmatch.fnmatch(entry, ib_pattern):
				ib_files_to_plot.append(datadir+entry)
			elif fnmatch.fnmatch(entry, an_pattern):
				an_files_to_plot.append(datadir + entry)

		#ib_storms_stats = np.genfromtxt(ib_files_to_plot[0])
		for ff, i in zip(ib_files_to_plot, range(len(ib_files_to_plot))):
			ib_data = np.genfromtxt(ff)
			print ib_data[lead_time_index]
			
			if np.isnan(ib_data[lead_time_index]):
				continue
			else:
				box_data.append(ib_data[lead_time_index])
			
		print box_data
			
		if y1 == 2006:
			boxdata20062007 = box_data
		elif y1 == 2007:
			boxdata20072008 = box_data
		elif y1 == 2008:
			boxdata20082009 = box_data
		elif y1 == 2009:
			boxdata20092010 = box_data
		elif y1 == 2010:
			boxdata20102011 = box_data
		elif y1 == 2011:
			boxdata20112012 = box_data
		elif y1 == 2012:
			boxdata20122013 = box_data
		elif y1 == 2013:
			boxdata20132014 = box_data
		elif y1 == 2014:
			boxdata20142015 = box_data
		elif y1 == 2015:
			boxdata20152016 = box_data
		elif y1 == 2016:
			boxdata20162017 = box_data
		elif y1 == 2017:
			boxdata20172018 = box_data
		elif y1 == 2018:
			boxdata20182019 = box_data

	
	plot_data=[boxdata20062007,boxdata20072008,boxdata20082009,boxdata20092010,boxdata20102011,boxdata20112012,boxdata20122013,boxdata20132014,boxdata20142015,boxdata20152016,boxdata20162017,boxdata20172018,boxdata20182019]
	
	

	#TO USE 5TH AND 95TH PERCENTILES, CHANGE whis=100 to whis=[5,95]
	bp = plt.boxplot(plot_data, sym='+', vert=1, whis=100, widths=0.25, patch_artist=True)
	plt.setp(bp['boxes'], color='black',linewidth=2)
	plt.setp(bp['whiskers'], color='black',linewidth=2,linestyle='solid')
	#plt.setp(bp['fliers'], color='red', marker='o',alpha=0.5)
	plt.setp(bp['medians'], color='black',linewidth=1)
	plt.setp(bp['caps'],color='black',linewidth=2)

	#box = plt.boxplot(data, patch_artist=True)
	#colors = ['#1A5899','#1A5899','#9C1127','#529DC8','#BBDAEA','#529DC8','#529DC8','white','#DB6B55','#DB6B55','#529DC8','#9C1127','#9C1127','#DB6B55' ,]
	colors = [c,c,c,c,c,c,c,c,c,c,c,c,c] 
	for patch, color,flier,median in zip(bp['boxes'], colors,bp['fliers'],bp['medians']):
		patch.set(facecolor=c,alpha=0.5)
		flier.set(marker='o', color=c, alpha=0.5,markersize=3)
		median.set(color='k',linewidth=2)


	#for label positions
	xpos=[1,2,3,4,5,6,7,8,9,10,11,12,13]

	my_xticks = ['2006-2007 ('+str(len(boxdata20062007))+')','2007-2008 ('+str(len(boxdata20072008))+')','2008-2009 ('+str(len(boxdata20082009))+')','2009-2010 ('+str(len(boxdata20092010))+')','2010-2011 ('+str(len(boxdata20102011))+')', '2011-2012 ('+str(len(boxdata20112012))+')','2012-2013 ('+str(len(boxdata20122013))+')','2013-2014 ('+str(len(boxdata20132014))+')','2014-2015 ('+str(len(boxdata20142015))+')','2015-2016 ('+str(len(boxdata20152016))+')','2016-2017 ('+str(len(boxdata20162017))+')','2017-2018 ('+str(len(boxdata20172018))+')','2018-2019 ('+str(len(boxdata20182019))+')']
	#my_yticks = ['0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60',' ','70',' ','80',' ','90',' ','100']
	plt.yticks(fontsize=12)
	plt.xticks(xpos, my_xticks,fontsize=12,rotation=40)
	#plt.xlabel('Cyclone Season', fontsize = 14)
	if stat_type == "location_error":
		plt.ylabel('Track Error (km)', fontsize=14)
		plt.ylim(bottom=0)
	elif stat_type == "mslp_bias":
		plt.ylabel('MSLP Bias (hPa)', fontsize=14)
	elif stat_type == "wind_bias":
		plt.ylabel(r'10m Wind Bias (ms$^{-1}$)', fontsize=14)
		
		
	spread = mpatches.Patch(color=c, alpha=0.5, label=str(lead_time_days)+' days lead time')
	legend = ax.legend(handles=[spread], fontsize=10, loc='upper right')
		


	plt.savefig('ukmo_det_2006-2019_'+str(stat_type)+'_each_year_box_plots_'+str(lead_time_days)+'_days_lead_time.png', bbox_inches='tight', dpi = 400)
	plt.close()

#plot_statistics("location_error")
plot_statistics("mslp_bias")
plot_statistics("wind_bias")



