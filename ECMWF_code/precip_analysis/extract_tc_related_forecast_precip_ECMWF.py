import numpy as np
from netCDF4 import Dataset
import picsea_library as pl
import sys
import os
from netCDF4 import Dataset, num2date
import datetime
import fnmatch

precipdir = "/vol/floods/more/rebecca_TC_do_not_delete/PICSEA_DATA/PRECIP/"
savedir = "/vol/floods/more/rebecca_TC_do_not_delete/PICSEA_DATA/TC-RELATED_PRECIP/"

#function to decide whether a given coordinate is within a given degree radius of a given point
def isInside(circle_x,circle_y, radius, x, y):	
	if ((x - circle_x)**2) + ((y - circle_y)**2) <= (radius**2):
		return True
	else:
		return False


#currently only have precip for 2013 - 2018...
year1s = [2012, 2018] #maybe 2010, 2011
year2s = [2013, 2019] #maybe 2011, 2012


for fcst_type in ["det"]: #"det","mean", "ctrl"
	print fcst_type
		
	for y1,y2 in zip(year1s, year2s):
		
		print y1, y2
		
		trackdir="/perm/mo/more/TIGGE/Y"+str(y1)+str(y2)+"/reformatted_track_files_per_storm_correct/SIO_storms/" 
		print trackdir
		
		#Find out the number of storms in this TC season and create array of the IDs so we can loop over them
		season_dirs=[]
		for root,dirs,files in os.walk(trackdir):
			for dir in dirs:
				print dir
				season_dirs.append(dir)
		NS = len(season_dirs) #total number of storms in this season in the SIO
		print NS
		
		precip_file = precipdir+"ECMWF_PRECIP_Jul"+str(y1)+"-Jun"+str(y2)+"_"+fcst_type+".nc"
		print precip_file
		precip_data = Dataset(precip_file, 'r')
		plons = precip_data.variables['longitude'][:]
		#print plons
		plats = precip_data.variables['latitude'][:]
		#print plats
		psteps = precip_data.variables['step'][:]
		#print psteps
						
		#structure of pcp is [date, time, step, lat, lon]
		#where date is YMD, time is 00 or 12, step is every 6 hours from 0 to 240 (41 steps)
		#This should match the forecast tracks, which start at the start date and go out to x steps every 6 hours
		#netcdf file too big to read in all the data, try reading in just the required day later... 
		#pcp = precip_data.variables['tp'][:]
						
		start_time = datetime.datetime(y1, 7, 1)
		
		tvalues = np.array([start_time + datetime.timedelta(days=i) for i in xrange(365)])
		
		pcp_dates = [j.strftime("%Y-%m-%d") for j in tvalues]
			
		lons_, lats_ = np.meshgrid(plons, plats)
		
			
		for dir, x in zip(season_dirs,range(len(season_dirs))):
			print trackdir+dir
			nwp_files_list = []
			analysis_file=0
			ibtracs_file=0
			#make a list of all the files in this directory
			list_of_all_files = os.listdir(trackdir+dir)
			#print list_of_all_files
			pattern="ecmwf_"+fcst_type+"*.txt"
			#pattern2="analysis_*.txt"
			#pattern3="ibtracs_*.txt"
			#for each file in the directory, find the ones that match the forecast filename, analysis filename and ibtracs filename
			for entry in list_of_all_files:
				if fnmatch.fnmatch(entry,pattern):
					nwp_files_list.append(trackdir+dir+"/"+entry)
				#elif fnmatch.fnmatch(entry,pattern2):
					#analysis_file = datadir+dir+"/"+entry
					#print analysis_file
				#elif fnmatch.fnmatch(entry,pattern3):
					#ibtracs_file = datadir+dir+"/"+entry
				else:
					continue
						
						
			print nwp_files_list
			
			for ff in nwp_files_list:
					
				track = np.genfromtxt(ff, dtype=float, skip_header=1)
				
				tracklons = track[:,7]
				tracklats = track[:,8]
				tracktime=track[:,6]
				
				print tracktime
				print type(tracktime[0])
					
				if tracktime[0] == 12:
					pcp_time = 1
				elif tracktime[0] == 0:
					pcp_time = 0
				else:
					print "ERROR: hang on, i thought the start time had to be 00 or 12???, but it was ", tracktime[0]
					
				fcst_datelist=pl.get_dates(track)
				print "fcst_datelist: ", fcst_datelist
				
				first_fcst_date = datetime.datetime.strptime(str(fcst_datelist[0]), "%Y%m%d%H")
			
				z = pcp_dates.index(first_fcst_date.strftime("%Y-%m-%d"))
				
				#Here, just read in the precip forecast data for this date's forecast, as file too big to read in all the data for the year
				pcp = precip_data.variables['tp'][z,:,:,:,:]
				
				this_fcst_tc_pcp = np.zeros((41, len(plats), len(plons)))
				
				for lt in range(1,len(tracklons)): #track might not go out to 41 timesteps, otherwise this would be lead time out to 41
				
					#print lt
					if lt < 41:
					
						#ECMWF archives accumulated precip through the whole forecast period, we just want to save the precip accumulated
						#in each timestep or the precip within 5 degrees will look strange as it will be the precip accumulated all the way up
						#to that lead timethroughout the forecast ,but only around that track points...
						#so get the precip accumulation just for the 6 hours that go with this track point...
						
						#Subtract the previous lead time's accumulation from this lead time's accumulation, 
						#to get the accumulation in the 6 hours since the last track point
						six_hour_pcp_accum = pcp[pcp_time, lt, :, :] - pcp[pcp_time, lt-1, :, :]
							
						track_point_lat = tracklats[lt]
						track_point_lon = tracklons[lt]
			
						#get the precip within a 5' radius of the track point at this lead time
						pcp_arr = np.ma.where(np.hypot(lons_-track_point_lon, lats_-track_point_lat) > 5, 0.,six_hour_pcp_accum)
						this_fcst_tc_pcp[lt,:,:] = pcp_arr
						
					else:
						continue
						#sometimes the track forecast goes out to further ahead than our det. precip forecast. how?
								
				
				#Save a tc-related precip file (nc) for each 'ff', to have the equivalent precip forecast as we have track forecasts
					
				outfile = Dataset(savedir+"ECMWF_det_"+fcst_datelist[0]+str(pcp_time)+"_"+dir+"_TC-related-PRECIP.nc", 'w', format='NETCDF4')
				outfile.description = "TC-related precip (precip within a 5' radius of a TC) for the forecast initiated on given date of given TC, for all lead times"
					
				outfile.createDimension('latitude', len(plats))
				outfile.createDimension('longitude', len(plons))
				outfile.createDimension('step', 41)
					
				latitude = outfile.createVariable('latitude', 'f4', 'latitude')
				longitude = outfile.createVariable('longitude', 'f4', 'longitude')
				step = outfile.createVariable('step', 'int', 'step')
				tc_pcp = outfile.createVariable('tc_pcp', 'float32', ('step', 'latitude', 'longitude'))
				
				latitude[:] = plats
				longitude[:] = plons
				step[:] = psteps
				tc_pcp[:,:,:] = this_fcst_tc_pcp
					
				latitude.units = 'degrees_north'
				longitude.units = 'degrees_east'
				step.units = 'hours'
				tc_pcp.units = 'm'
				tc_pcp.long_name = 'Accumulated precip since start time'
					
				outfile.close()
					
				print "file saved: ", savedir+"ECMWF_det_"+fcst_datelist[0]+str(pcp_time)+"_"+dir+"_TC-related-PRECIP.nc"
							

		precip_data.close()	
			
		

 


