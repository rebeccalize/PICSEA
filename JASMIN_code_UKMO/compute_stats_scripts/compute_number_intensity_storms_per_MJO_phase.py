import track_statistics_library as ts
import picsea_library as pl
import sys
import os
import fnmatch
import numpy as np
from netCDF4 import Dataset, num2date, netcdftime
import datetime

#This script, for one TC season (e.g. 2015-2016) in the SIO:
#Uses functions in track_statistics_library.py, to:
#Compute the track location error for every consecutive deterministic forecast of each storm (great circle distance in km), at each lead time
#Compute the intensity bias (vorticity) for every consecutive deterministic forecast of each storm, at each lead time
#Compute the average track location error and intensity bias across all the forecast for each storm
#Compute the average track location error and intensity bias across all the forecasts and across all the storms in this TC season

#y1=sys.argv[1]
#y2=sys.argv[2]
MJO=sys.argv[1] #MJO phase / phase pair [12,23,34,45,56,67,78,81]


array_len=31 #number of forecast timesteps (7 days x 4 forecasts)

#CURRENTLY RUNNING OTHER SH STORMS, NOT SIO - CHECK BEFORE RUNNING

#FOR DIFFERENT SEASONS, or across like 2006 - 2016, use this:

#savedir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/MJO_phase/"

#FOR LANDFALLING CYCLONES IN DIFFERENT COUNTRIES, USE THESE:
#datadir="/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/mozambique_landfalling/" #SIO_storms/"
#savedir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_and_intensity_errors/mozambique_landfalling/"

#if not os.path.exists(savedir):
    #os.makedirs(savedir)


	
	
MJO_file = "/gws/nopw/j04/klingaman/datasets/MJO_INDICES/MJO_phase"+str(MJO)+".jan-dec_dmeans_ts.1979-2019.nc"
print MJO_file

ffMJO = Dataset(MJO_file,'r')
MJOamp = ffMJO.variables['phase_ts'][:]
MJOdatesnc = ffMJO.variables['time'][:]

t_unit = ffMJO.variables['time'].units
t_cal = ffMJO.variables['time'].calendar
tvalue = num2date(MJOdatesnc,units=t_unit, calendar=t_cal)

MJOdates = [i.strftime("%Y-%m-%d") for i in tvalue] #converts dates from netcdf date format to python's datetime format - bit easier to use later


number_of_storms=[]
storm_max_winds=[]
storm_min_mslp=[]
	
	
def mjo_tc_freq_int(ibtracs_file, analysis_file, trno):

	ib_file=ibtracs_file
	ib_data=np.genfromtxt(ib_file, dtype=float, skip_header=1)
	ib_lon=ib_data[:,7]
	ib_lat=ib_data[:,8]
	ib_mslp=ib_data[:,9]
	ib_wind=ib_data[:,10]
	ib_datelist = pl.get_dates(ib_data)
	
	an_file=analysis_file
	an_data=np.genfromtxt(an_file, dtype=float, skip_header=1)
	an_lon=an_data[:,7]
	an_lat=an_data[:,8]
	an_mslp=an_data[:,9]
	an_wind=an_data[:,10]
	an_datelist = pl.get_dates(an_data)
	
	#want to get the tracks that started in a certain MJO phase, but based on the whole lifecycle (i.e. the analysis) rather than just TC stages
	first_date = datetime.datetime.strptime(str(an_datelist[0]), "%Y%m%d%H")	
			
			
	#print first_date.strftime("%m-%d")
			
	#then find the index of MJO dates at which the date is the first date of the analysis track
	if not first_date.strftime("%m-%d") == '02-29': #leap years missing in the MJO data?
			
		z = MJOdates.index(first_date.strftime("%Y-%m-%d"))
		
	
			
	#then find out the MJO amplitude on this date (i.e. at the index of this date)
	amp = MJOamp[z]
	#print "MJO amp", amp
			
	#if the MJO amplitude is less than 1, don't include this foreacst
	#we only want to compute the errors for forecasts started during this MJO phase / phase pair, when the MJO amplitude is >1
	#if amp < 1.0:
		#continue
				
	#if the MJO amplitude was >1 on this day, add 1 to the number of storms included for this MJO phase pair
	#and get the intensity data from the ibtracs track, to plot the range of max intensities of the storms in this MJO phase pair
	if amp >= 1.0:
	
		number_of_storms.append(1)
		storm_max_winds.append(np.nanmax(an_wind))
		storm_min_mslp.append(np.nanmin(an_mslp))
		
		print trno
		
		
		
year1s=[2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017]
year2s=[2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018]



for y1,y2 in zip(year1s,year2s):

	datadir="/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/" #SIO_storms/"

	print y1,y2

	#Find out the number of storms in this TC season and create array of the IDs so we can loop over them
	season_dirs=[]
	for root,dirs,files in os.walk(datadir):
		for dir in dirs:
			season_dirs.append(dir)
	NS = len(season_dirs) #total number of storms in this season in the SIO

	#loop over storm directories in this TC season, and in each directory find each ibtracs file
	#then call the function for each storm using the identified files
	for dir, x in zip(season_dirs,range(len(season_dirs))):
		#print datadir+dir

		ibtracs_file=0
		#make a list of all the files in this directory
		list_of_all_files = os.listdir(datadir+dir)

		pattern_ib="ibtracs_*.txt"
		pattern_an="analysis_*.txt"
		
		#for each file in the directory, find the ones that match the forecast filename, analysis filename and ibtracs filename
		for entry in list_of_all_files:
			if fnmatch.fnmatch(entry,pattern_ib):
				ibtracs_file = datadir+dir+"/"+entry
			elif fnmatch.fnmatch(entry,pattern_an):
				analysis_file = datadir+dir+"/"+entry
		trno=dir

		#run the statistics for this storm
		mjo_tc_freq_int(ibtracs_file, analysis_file, trno)
		
print number_of_storms
print storm_max_winds
print storm_min_mslp

for i in range(len(storm_max_winds)):
	if storm_max_winds[i] > 10000:
		storm_max_winds[i] = np.nan
		
	if storm_min_mslp[i] > 10000:
		storm_min_mslp[i] = np.nan

np.savetxt("intensity_of_TCs_ANALYSIS_072006-062018_MJO_phase_"+str(MJO)+".txt", np.c_[storm_max_winds,storm_min_mslp],'%.4f')
		




