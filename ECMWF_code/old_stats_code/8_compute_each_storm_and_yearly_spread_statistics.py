import track_statistics_library as ts
import picsea_library as pl
import sys
import os
import fnmatch
import numpy as np

es=[3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52]

#This script calculates the spread-error relationship for location, intensity and translation speed
#Spread is calculated by first taking the difference between each ensemble member and the ensemble mean
#And then taking the average of this difference across all the ensemble members, at each timestep
#Then plot this alongside the error in the ensemble mean

y1=sys.argv[1]
y2=sys.argv[2]

array_len=42


datadir="/perm/mo/more/TIGGE/Y"+str(y1)+str(y2)+"/reformatted_track_files_per_storm_WITH_INTENSITY/SIO_storms/"

savedir = "/perm/mo/more/picsea/forecast_spread/Y"+str(y1)+str(y2)+"/SIO_storms/"

#Find out the number of storms in this TC season and create array of the IDs so we can loop over them
season_dirs=[]
for root,dirs,files in os.walk(datadir):
	for dir in dirs:
		season_dirs.append(dir)
NS = len(season_dirs) 




all_storms_avg_loc_spread, all_storms_avg_wind_spread,all_storms_avg_mslp_spread,all_storms_avg_speed_spread, all_storms_wgt = (np.zeros(array_len) for i in range(5))

all_storms_ind_loc_spread, all_storms_ind_wind_spread, all_storms_ind_mslp_spread, all_storms_ind_speed_spread = (np.zeros((NS,array_len)) for i in range(4))


def storm_stats(eps_files_list, mean_files_list, array_len, NT_eps, storm_no):
	
	storm_loc_spread, storm_wind_spread, storm_mslp_spread, storm_speed_spread, storm_wgt = (np.zeros(array_len) for i in range(5))
		
	#first put the avg diff ebtween each ens member and the mean, for each forecast track for this storm, in here
	#then afterwards, average across the avg difference for each forecast track, to get the average for each storm
	#always at each lead time...
	each_forecast_loc_spread, each_forecast_wind_spread, each_forecast_mslp_spread, each_forecast_speed_spread, each_forecast_wgt = (np.zeros((NT_eps, array_len)) for i in range(5))
		
		
	#load in the ensemble mean file for this date
		
	#load in the ensemble file for this date
	#then loop over the ensemble forecast tracks and calculate the difference between that and the mean track
		
	#for each pair of forecasts, calculate the difference between the two forecasts
	for ff_mean,ff_eps, tno in zip(mean_files_list, eps_files_list, range(NT_eps)):
	
		print ff_mean
		print ff_eps
		
		all_ens_members_loc_diff = np.zeros((50, array_len))
		all_ens_members_wind_diff = np.zeros((50, array_len))
		all_ens_members_mslp_diff = np.zeros((50, array_len))
		all_ens_members_speed_diff = np.zeros((50, array_len))
		
		mean_data=np.genfromtxt(ff_mean,dtype=float,skip_header=1,usecols=np.arange(0,11))
		mean_lon=mean_data[:,7]
		mean_lat=mean_data[:,8]
		mean_mslp=mean_data[:,9]
		mean_wind=mean_data[:,10]
		mean_datelist = pl.get_dates(mean_data)
		
		mean_speed = ts.prop_speed_vals(mean_lon, mean_lat)
		print "mean_speed: ", mean_speed
		ML = len(mean_speed)
		print ML
		print "mean dates: ", mean_datelist
		
		#here, will need to load in separate file with the mean forecast translation speed
		
			
		eps_data=np.genfromtxt(ff_eps,dtype=float,skip_header=1,usecols=np.arange(0,11))
		
		#here, will need to load in separate file with the eps forecast translation speed
		
		for e, j in zip(es, range(len(es))):
		
			print "ensemble member number ", e,j
		
			eps_lon=eps_data[np.where(eps_data[:,0] == e),7][0]
			eps_lat=eps_data[np.where(eps_data[:,0] == e),8][0]
			eps_mslp=eps_data[np.where(eps_data[:,0] == e),9][0]
			eps_wind=eps_data[np.where(eps_data[:,0] == e),10][0]
			eps_datelist = pl.get_dates(eps_data[np.where(eps_data[:,0] == e),:][0])
			
			#print "eps_lon", eps_lat
			#print "eps_lat", eps_lat
				
			eps_speed = ts.prop_speed_vals(eps_lon, eps_lat)
			print "eps_speed: ", eps_speed	
			EL = len(eps_speed)
			print EL
			print "eps dates: ", eps_datelist
			
			fcst_len = np.min([ML,EL])
			print "fcst_len: ", fcst_len
			
			if fcst_len > 41:
				fcst_len = 41
				
			if not fcst_len == 0:
			
			
				#loc_diff, wind_diff, mslp_diff, speed_diff = (np.zeros(array_len) for i in range(4))
				for lt in range(fcst_len+1):
					#print lt
			
					if lt == fcst_len:
						all_ens_members_loc_diff[j,lt] = ts.trerr([mean_lon[lt],mean_lat[lt]], [eps_lon[lt],eps_lat[lt]])
						all_ens_members_wind_diff[j,lt] = abs(mean_wind[lt] - eps_wind[lt])
						all_ens_members_mslp_diff[j,lt] = abs(mean_mslp[lt] - eps_mslp[lt])
					else:
						all_ens_members_loc_diff[j,lt] = ts.trerr([mean_lon[lt],mean_lat[lt]], [eps_lon[lt],eps_lat[lt]])
						all_ens_members_wind_diff[j,lt] = abs(mean_wind[lt] - eps_wind[lt])
						all_ens_members_mslp_diff[j,lt] = abs(mean_mslp[lt] - eps_mslp[lt])
						all_ens_members_speed_diff[j,lt] = ts.prop_speed_abs_err(mean_speed[lt], eps_speed[lt])
					
				
		for lt in range(array_len):
			for e,j in zip(es, range(len(es))):
				if not np.isnan(all_ens_members_loc_diff[j,lt]):
					storm_wgt[lt] += 1
					all_storms_wgt[lt] += 1

			each_forecast_loc_spread[tno,lt] = np.nanmean(all_ens_members_loc_diff[:,lt])
			each_forecast_wind_spread[tno,lt] = np.nanmean(all_ens_members_wind_diff[:,lt])
			each_forecast_mslp_spread[tno,lt] = np.nanmean(all_ens_members_mslp_diff[:,lt])
			each_forecast_speed_spread[tno,lt] = np.nanmean(all_ens_members_speed_diff[:,lt])
			
	
	for lt in range(array_len):
		
			
		storm_loc_spread[lt] = np.nanmean(each_forecast_loc_spread[:,lt])
		storm_wind_spread[lt] = np.nanmean(each_forecast_wind_spread[:,lt])
		storm_mslp_spread[lt] = np.nanmean(each_forecast_mslp_spread[:,lt])
		storm_speed_spread[lt] = np.nanmean(each_forecast_speed_spread[:,lt])
		
	np.savetxt(savedir + dir + "_loc_spread_per_lead_time.txt", storm_loc_spread[:], '%.4f')
	np.savetxt(savedir + dir + "_wind_spread_per_lead_time.txt", storm_wind_spread[:], '%.4f')
	np.savetxt(savedir + dir + "_mslp_spread_per_lead_time.txt", storm_mslp_spread[:], '%.4f')
	np.savetxt(savedir + dir + "_speed_spread_per_lead_time.txt", storm_speed_spread[:], '%.4f')
	np.savetxt(savedir + dir + "_number_ens_members_included_in_spread_calcs_per_lead_time.txt", storm_wgt[:], '%.4f')
	
		
	all_storms_ind_loc_spread[storm_no,:] = storm_loc_spread[:]
	all_storms_ind_wind_spread[storm_no,:] = storm_wind_spread[:]
	all_storms_ind_mslp_spread[storm_no,:] = storm_mslp_spread[:]
	all_storms_ind_speed_spread[storm_no,:] = storm_speed_spread[:]
		

#this bit should go after the function really, just trying to get my head around it
for dir, x in zip(season_dirs,range(len(season_dirs))):
	
	orig_eps_files_list = []
	orig_mean_files_list = []
	matched_eps_files_list = []
	matched_mean_files_list = []
	
	list_of_all_files = os.listdir(datadir+dir)
	#print list_of_all_files
	pattern_eps="ecmwf_ens*.txt"
	pattern_mean = "ecmwf_mean*.txt"
	
	#for each file in the directory, find the ones that match the forecast filename, analysis filename and ibtracs filename
	for entry in list_of_all_files:
		if fnmatch.fnmatch(entry,pattern_eps):
			orig_eps_files_list.append(datadir+dir+"/"+entry)
		elif fnmatch.fnmatch(entry,pattern_mean):
			orig_mean_files_list.append(datadir+dir+"/"+entry)
		else:
			continue
	
	#sometimes, there are more mean files than eps files, or vice versa perhaps
	#this bit finds the date in the filename, goes through the mean files looking for the same date, and if it finds a match, saves both those filenames for use
	#so if there's a mean file for a date that doesn't have an ensemble forecast, we don't use it 
	#although they don't end up in the same order, this shouldn't really matter
	
	eps_dates=[]
	mean_dates=[]
	
	for eps_filename in orig_eps_files_list:
		#print eps_filename
		eps_date = eps_filename[107:117] #characters of the date in the filename, which includes the dir string
		
		for mean_filename in orig_mean_files_list:
			mean_date = mean_filename[108:118]
			if mean_date == eps_date:
				matched_mean_files_list.append(mean_filename)
				mean_dates.append(mean_date)
				
				matched_eps_files_list.append(eps_filename)
				eps_dates.append(eps_date)
			else:
				#print mean_filename," not found in eps_file_list"
				continue
				
			
		
	print eps_dates
	print mean_dates
	
	NT_eps = len(matched_eps_files_list)
	
	print NT_eps

	
	#run the statistics for this storm
	storm_stats(matched_eps_files_list, matched_mean_files_list, array_len, NT_eps, x)
		
		

		
for lt in range(array_len):

	all_storms_avg_loc_spread[lt] = np.nanmean(all_storms_ind_loc_spread[:,lt])
	all_storms_avg_wind_spread[lt] = np.nanmean(all_storms_ind_wind_spread[:,lt])
	all_storms_avg_mslp_spread[lt] = np.nanmean(all_storms_ind_mslp_spread[:,lt])
	all_storms_avg_speed_spread[lt] = np.nanmean(all_storms_ind_speed_spread[:,lt])
	
	
	
np.savetxt(savedir + "average_loc_spread_per_lead_time.txt", all_storms_avg_loc_spread[:], '%.4f')
np.savetxt(savedir + "average_wind_spread_per_lead_time.txt", all_storms_avg_wind_spread[:], '%.4f')
np.savetxt(savedir + "average_mslp_spread_per_lead_time.txt", all_storms_avg_mslp_spread[:], '%.4f')
np.savetxt(savedir + "average_speed_spread_per_lead_time.txt", all_storms_avg_speed_spread[:], '%.4f')


np.savetxt(savedir + "overall_number_ens_members_included_in_spread_calcs_per_lead_time.txt", all_storms_wgt[:], '%.4f')	
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		 
	
	
		
	
		



