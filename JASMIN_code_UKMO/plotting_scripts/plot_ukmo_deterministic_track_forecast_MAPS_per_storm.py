#writing function to plot map of SOI / SE Africa with TC tracks (initially from IBTrACS)

import picsea_library as pl
import fnmatch
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.cm as cm
import numpy as np

datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_TIGGE_analysis/reformatted_track_files_per_storm/Y"

#plot the UKMO NWP tracks for one storm at a time, alongside reanalysis and ibtracs tracks
y1=sys.argv[1]
y2=sys.argv[2]

directory = datadir+str(y1)+str(y2)+"/"
print directory

#cat_file = datadir+str(y1)+"_"+str(y2)+"/SIO_storms/storm_categories.txt"
#cf = open(cat_file,'r')
#lines=cf.readlines()
#trnos, cats, names = ([] for i in range(3))
#for x in lines:
	#print x
	#trnos.append(x.split(' ')[0])
	#cats.append(x.split(' ')[1])
	#namen = x.split(' ')[2]
	#names.append(namen.rstrip())


#find all the track files in this TC season (e.g. 2015-2016)
for root,dirs,files in os.walk(directory): #+"/SIO_storms"
	#loop over the directories (i.e. loop over the storms)
	for dir in dirs:
		print dir
		#if dir in trnos:
			#si = [trnos.index(i) for i in trnos if dir in i]
			#si = si[0]
			#cat = cats[si]
			#if cat == "0":
			#cat = "TS"
			#name = names[si]
		#else:
		cat = np.ma.masked
		name = np.ma.masked


		nwp_infile_list = []
		analysis_file=0
		ibtracs_file=0
		#make a list of all the files in this directory
		list_of_all_files = os.listdir(directory+"/"+dir+"/") #
		pattern="ukmo_tigge_mean*.txt"
		pattern2="analysis_*.txt"
		pattern3="ibtracs_*.txt"
		#print "directory:", dir
		#for each file in the directory, find the ones that match the forecast filename, analysis filename and ibtracs filename
		for entry in list_of_all_files:
			if fnmatch.fnmatch(entry,pattern):
				nwp_infile_list.append(directory+dir+"/"+entry) #+"/SIO_storms/"
			if fnmatch.fnmatch(entry,pattern2):
				analysis_file = directory+dir+"/"+entry #+"/SIO_storms/"
				#print analysis_file
			elif fnmatch.fnmatch(entry,pattern3):
				ibtracs_file = directory+dir+"/"+entry #+"/SIO_storms/"
				#print ibtracs_file
		#print nwp_infile_list
		outfile = directory+"plot_of_ensmean_tracks_"+dir+".png" #/SIO_storms
		#print outfile

		#create the plot for this storm
		pl.map_nwp_tracks_per_storm(nwp_infile_list, analysis_file, ibtracs_file, outfile, "SIO",cat,name)
		print "Saved ", dir
        
        
        
        
        
        
        
        
        
        
    

