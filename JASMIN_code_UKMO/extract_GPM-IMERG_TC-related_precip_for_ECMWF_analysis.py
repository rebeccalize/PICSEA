import numpy as np
from netCDF4 import Dataset
import picsea_library as pl
import sys
import os
from netCDF4 import Dataset, num2date, date2num
import datetime
import fnmatch

precipdir = "/gws/nopw/j04/klingaman/datasets/GPM_IMERG/six_hourly/"
savedir = "/gws/nopw/j04/klingaman/emerton/GPM_IMERG_TC-related_precip/"

#function to decide whether a given coordinate is within a given degree radius of a given point
def isInside(circle_x,circle_y, radius, x, y):	
	if ((x - circle_x)**2) + ((y - circle_y)**2) <= (radius**2):
		return True
	else:
		return False


#currently only have precip for 2013 - 2018...
year1s = [2016, 2017, 2018] #maybe 2019
year2s = [2017, 2018, 2019] #maybe 2020


for obs_track in ["analysis"]: #, "analysis"

		
	for y1,y2 in zip(year1s, year2s):
		
		print y1, y2
		
		trackdir="/gws/nopw/j04/klingaman/emerton/reformatted_sh_track_output/UKMO/HRES/reformatted_track_files_by_storm/"+str(y1)+"_"+str(y2)+"/SIO_storms/" 
		print trackdir
		
		#Find out the number of storms in this TC season and create array of the IDs so we can loop over them
		season_dirs=[]
		for root,dirs,files in os.walk(trackdir):
			for dir in dirs:
				print dir
				season_dirs.append(dir)
		NS = len(season_dirs) #total number of storms in this season in the SIO
		print NS
		
		#precip_file = precipdir+"ECMWF_PRECIP_Jul"+str(y1)+"-Jun"+str(y2)+"_"+fcst_type+".nc"
		#print precip_file
		#precip_data = Dataset(precip_file, 'r')
		#plons = precip_data.variables['longitude'][:]
		#print plons
		#plats = precip_data.variables['latitude'][:]
		#print plats
		#psteps = precip_data.variables['step'][:]
		#print psteps
						
		#structure of pcp is [date, time, step, lat, lon]
		#where date is YMD, time is 00 or 12, step is every 6 hours from 0 to 240 (41 steps)
		#This should match the forecast tracks, which start at the start date and go out to x steps every 6 hours
		#netcdf file too big to read in all the data, try reading in just the required day later... 
		#pcp = precip_data.variables['tp'][:]
						
		#start_time = datetime.datetime(y1, 7, 1)
		
		#tvalues = np.array([start_time + datetime.timedelta(days=i) for i in xrange(365)])
		
		#pcp_dates = [j.strftime("%Y-%m-%d") for j in tvalues]
			
		#lons_, lats_ = np.meshgrid(plons, plats)
		
			
		for dir, x in zip(season_dirs,range(len(season_dirs))):
			print trackdir+dir
			nwp_files_list = []
			analysis_file=0
			ibtracs_file=0
			#make a list of all the files in this directory
			list_of_all_files = os.listdir(trackdir+dir)
			#print list_of_all_files
			#pattern="ecmwf_"+fcst_type+"*.txt"
			#pattern2="analysis_*.txt"
			pattern = obs_track+"_*.txt"
			#for each file in the directory, find the ones that match the forecast filename, analysis filename and ibtracs filename
			for entry in list_of_all_files:
				if fnmatch.fnmatch(entry,pattern):
					#nwp_files_list.append(trackdir+dir+"/"+entry)
					track_file = trackdir+dir+"/"+entry
				#elif fnmatch.fnmatch(entry,pattern2):
					#analysis_file = datadir+dir+"/"+entry
					#print analysis_file
				#elif fnmatch.fnmatch(entry,pattern3):
					#ibtracs_file = datadir+dir+"/"+entry
				else:
					continue
						
						
			
					
			track = np.genfromtxt(track_file, dtype=float, skip_header=1)
				
			tracklons = track[:,7]
			tracklats = track[:,8]
			trackyears = track[:,3]
			trackmonths = track[:,4]
			trackdays = track[:,5]
			tracktime = track[:,6]
			
			print tracktime
			print type(tracktime[0])
					
				
			track_datelist=pl.get_dates(track)
			print "track_datelist: ", track_datelist
			
			first_track_date = datetime.datetime.strptime(str(track_datelist[0]), "%Y%m%d%H")
			
			this_track_pcp = np.zeros((len(track_datelist), 1800, 3600))
			
			for date, i in zip(track_datelist, range(len(track_datelist))):
			
				if int(trackmonths[i]) <= 6:
					yyyy = y2
				elif int(trackmonths[i]) >= 7:
					yyyy = y1
				
				print trackyears[i]
				print trackmonths[i]
				print trackdays[i]
				
				precip_file = precipdir+str(yyyy)+"/3B-HHR.MS.MRG.3IMERG."+str(int(trackyears[i]))+str(int(trackmonths[i])).zfill(2)+str(int(trackdays[i])).zfill(2)+".6hr_means.V06B.nc"
				
				print precip_file
				precip_data = Dataset(precip_file, 'r')
				plons = precip_data.variables['longitude'][:]
				#print plons
				plats = precip_data.variables['latitude'][:]
				#print plats
						
				#structure of pcp is [time, lat, lon]
				#where time is 0, 6, 12, 18, i.e. the average mm/day in each 6-hour period
				
				pcp = precip_data.variables['precipitationCal'][:]
						
				#start_time = datetime.datetime(y1, 7, 1)
		
				#tvalues = np.array([start_time + datetime.timedelta(days=i) for i in xrange(365)])
		
				#pcp_dates = [j.strftime("%Y-%m-%d") for j in tvalues]
			
				lons_, lats_ = np.meshgrid(plons, plats)
				
				if int(tracktime[i]) == 0:
					z = 0
				elif int(tracktime[i]) == 6:
					z = 1
				elif int(tracktime[i]) == 12:
					z = 2
				elif int(tracktime[i]) == 18:
					z = 3
					
				
				pcp = pcp[z,:,:]
				
				track_point_lat = tracklats[i]
				track_point_lon = tracklons[i]
				
				pcp_arr = np.ma.where(np.hypot(lons_-track_point_lon, lats_-track_point_lat) > 5, 0., pcp)
				this_track_pcp[i,:,:] = pcp_arr
				
				precip_data.close()
				
				
			
			#Save a tc-related precip file (nc) for each track, to have the equivalent GPM IMERG precipas we have ibtracs/analysis trac
				
			outfile = Dataset(savedir+"GPM-IMERG_"+obs_track+"_Y"+str(y1)+"-"+str(y2)+"_"+dir+"_startdate_"+str(int(trackyears[0]))+str(int(trackmonths[0])).zfill(2)+str(int(trackdays[0])).zfill(2)+str(int(tracktime[0])).zfill(2)+"_TC-related-PRECIP.nc", 'w', format='NETCDF4')
			outfile.description = "TC-related precip (precip within a 5' radius of a TC) from GPM IMERG associated with this observed track of given TC,from start to end date of TC every 6 hours."
				
			outfile.createDimension('latitude', len(plats))
			outfile.createDimension('longitude', len(plons))
			outfile.createDimension('time', len(track_datelist))
		
				
			latitude = outfile.createVariable('latitude', 'f4', 'latitude')
			longitude = outfile.createVariable('longitude', 'f4', 'longitude')
			time = outfile.createVariable('time', 'int', 'time')
			tc_pcp = outfile.createVariable('tc_pcp', 'float32', ('time', 'latitude', 'longitude'))
			
			
			#create proper netcdf times
			time.units = "hours since 1970-01-01 00"
			times = []
			for t in range(len(track_datelist)):
				times.append(datetime.datetime(int(trackyears[t]),int(trackmonths[t]),int(trackdays[t]),int(tracktime[t])))
			times = date2num(times, time.units)
			
			#add the data to the file :) 
			time[:] = times
			latitude[:] = plats
			longitude[:] = plons
			tc_pcp[:,:,:] = this_track_pcp
		
			#add units etc to the file		
			latitude.units = 'degrees_north'
			longitude.units = 'degrees_east'
			tc_pcp.units = 'mm/hr'
			tc_pcp.long_name = 'Average 6-hourly precipitation in mm/hr (from calibrated GPM-IMERG 3B-HHR V06B)'
			
			#this num2date(time[:],time.units) returns the correct dates from the netcdf file variable...
			#print num2date(time[:],time.units)
				
			outfile.close()
				
			print "file saved: ", savedir+"GPM-IMERG_"+obs_track+"_"+str(y1)+str(y2)+"_"+dir+"_TC-related-PRECIP.nc"
							

		#precip_data.close()	
			
		

 


