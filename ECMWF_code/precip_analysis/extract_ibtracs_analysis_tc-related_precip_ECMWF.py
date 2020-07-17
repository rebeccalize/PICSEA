import numpy as np
from netCDF4 import Dataset
import picsea_library as pl
import sys
import os
from netCDF4 import Dataset, num2date
import datetime
import fnmatch

precipdir = "/scratch/mo/more/GPM-IMERG-FILES-CORRECT/"
savedir = "/vol/floods/more/rebecca_TC_do_not_delete/PICSEA_DATA/TC-RELATED_PRECIP/"

#function to decide whether a given coordinate is within a given degree radius of a given point
def isInside(circle_x,circle_y, radius, x, y):	
	if ((x - circle_x)**2) + ((y - circle_y)**2) <= (radius**2):
		return True
	else:
		return False


#currently only have precip for 2013 - 2018...
year1s = [2013] #maybe 2010, 2011
year2s = [2014] #maybe 2011, 2012


for track_type in ["ibtracs", "analysis"]: #"det","mean", "ctrl"
	print track_type
		
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
		
		precip_file = precipdir+"GPM_IMERG_DAILY.TRMM_grid.Jul"+str(y1)+"-Jun"+str(y2)+".nc"
		print precip_file
		precip_data = Dataset(precip_file, 'r')
		plons = precip_data.variables['longitude'][:]
		#print plons
		plats = precip_data.variables['latitude'][:]
		#print plats
		
		ptime = precip_data.variables['time'][:]
	
						
		#structure of pcp is [date, time, step, lat, lon]
		#where date is YMD, time is 00 or 12, step is every 6 hours from 0 to 240 (41 steps)
		#This should match the forecast tracks, which start at the start date and go out to x steps every 6 hours
		#netcdf file too big to read in all the data, try reading in just the required day later... 
		#pcp = precip_data.variables['tp'][:]
						
		start_time = datetime.datetime(y1, 7, 1)
		
		tvalues = np.array([start_time + datetime.timedelta(days=i) for i in xrange(len(ptime))])
		
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
			#pattern="ecmwf_"+fcst_type+"*.txt"
			if track_type == "ibtracs":
				pattern = "ibtracs_*.txt"
			elif track_type == "analysis":
				pattern = "analysis_*.txt"
		
			#for each file in the directory, find the ones that match the forecast filename, analysis filename and ibtracs filename
			for entry in list_of_all_files:
				if fnmatch.fnmatch(entry,pattern):
				
					track_file = trackdir+dir+"/"+entry
			
				else:
					continue
						
			
					
			track = np.genfromtxt(track_file, dtype=float, skip_header=1)
				
			tracklons = track[:,7]
			tracklats = track[:,8]
			tracktime = track[:,6]
			
			print tracktime
			print type(tracktime[0])
			
			track_datelist=pl.get_dates(track)
			print "track_datelist: ", track_datelist
			
			first_track_date = datetime.datetime.strptime(str(track_datelist[0]), "%Y%m%d%H")
		
			#z = pcp_dates.index(first_track_date.strftime("%Y-%m-%d"))
			
			#Here, just read in the precip forecast data for this date rather than the whole file
			#pcp = precip_data.variables['tp'][z,:,:]
			
			this_track_tc_pcp = np.zeros((len(tracktime), len(plats), len(plons)))
			
			
			
			print "length: ", len(track_datelist)
			for date, di in zip(track_datelist, range(len(track_datelist))):
			
				print "di = ", di
			
				this_date = datetime.datetime.strptime(str(track_datelist[di]), "%Y%m%d%H")
			
				#The first track point could be at 00 or 06 or 12 or 18 hours - but the daily precip runs 00 - 00
				#for a normal day, we want to loop over the 4 track points and get the rainfall around each
				#But the first track point might not have 4 points in the day, 
				#it might start at 18 so only have 2 points (at 18 and 00)
				#or it might start at 06 and have 4 points (at 06, 12, 18 and 00)
				#it might start at 12 and have 3 points (12, 18 and 00)
				#or it might start at 00 and have 1 point (00)
				
				if di == 0:
					print "first date: ", this_date
					z = pcp_dates.index(this_date.strftime("%Y-%m-%d"))
					pcp = precip_data.variables['precipitationCal'][z,:,:]
					pcp = np.transpose(pcp)
				
					if this_date.strftime("%H") == '00':
						timesteps = 1
						offset = 0
					elif this_date.strftime("%H") == '06':
						timesteps = 4
						offset = 3
					elif this_date.strftime("%H") == '12':
						timesteps = 3
						offset = 2
					elif this_date.strftime("%H") == '18':
						timesteps = 2
						offset = 1
						
					#circles=[]
						
					for ts in range(timesteps):
					
						track_point_lat = tracklats[ts]
						track_point_lon = tracklons[ts]
						
						print track_point_lat, track_point_lon
						
						#circle = np.hypot(lons_-track_point_lon, lats_-track_point_lat)
						#circles.append(circles)
						
						print "shape circle: ",np.shape(np.hypot(lons_-track_point_lon, lats_-track_point_lat))
						print "shape pcp: ",np.shape(pcp)
						
						#return elements from x (0) where we are outside 5 degree circle, otherwise return value from pcp
						this_ts_pcp = np.ma.where(np.hypot(lons_-track_point_lon, lats_-track_point_lat) > 5, 0., pcp)
						
						for x in range(len(plats)):
							for y in range(len(plons)):
							
								if this_ts_pcp[x,y] > 0.:
									print x, y, this_ts_pcp[x,y]
									
									this_track_tc_pcp[di+offset,x,y] = this_ts_pcp[x,y]
									
								
										
						
							
						
					#print circles[0]
					#pcp_arr = np.ma.where(any(circles) > 5, 0., pcp)
					#this_track_tc_pcp[di+offset] = pcp_arr
					#print pcp_arr
					
				
				elif di == len(track_datelist) - 1:
				
					print this_date
				
					z = pcp_dates.index(this_date.strftime("%Y-%m-%d"))
					pcp = precip_data.variables['precipitationCal'][z,:,:]
					pcp = np.transpose(pcp)
				
					if this_date.strftime("%H") == '00':
						continue
					elif this_date.strftime("%H") == '06':
						tsa = [di]
					
					elif this_date.strftime("%H") == '12':
						tsa = [di-1, di]
						
					elif this_date.strftime("%H") == '18':
						tsa = [di-2, di-1, di]
						
						
					#circles=[]
						
					for ts in tsa:
					
						track_point_lat = tracklats[ts]
						track_point_lon = tracklons[ts]
						
						#circle = np.hypot(lons_-track_point_lon, lats_-track_point_lat)
						#circles.append(circles)
						
						this_ts_pcp = np.ma.where(np.hypot(lons_-track_point_lon, lats_-track_point_lat) > 5, 0., pcp)
						
						for x in range(len(plats)):
							for y in range(len(plons)):
							
								if this_ts_pcp[x,y] > 0.:
									
									this_track_tc_pcp[di,x,y] = this_ts_pcp[x,y]
									
								
					#print circles
					#pcp_arr = np.ma.where(any(circles) > 5, 0., pcp)
					#this_track_tc_pcp[di] = pcp_arr
					#print pcp_arr
						

				#If this isn't the first date in list, check if the day is the same as the one before,
				#If the day is the same as the one before, we don't want to run it again
				else:
				
					last_date = datetime.datetime.strptime(str(track_datelist[di-1]), "%Y%m%d%H")
					if this_date.strftime("%Y-%m-%d") == last_date.strftime("%Y-%m-%d"):
						
						print "dates the same: ", this_date, last_date
						
						continue
						
					else:
					
						if di >= len(track_datelist) - 4:
							continue
						else:
							print "getting precip", this_date
							z = pcp_dates.index(this_date.strftime("%Y-%m-%d"))
							pcp = precip_data.variables['precipitationCal'][z,:,:]
							pcp = np.transpose(pcp)
						
							#circles = []
						
							for ts in [di +1, di+2, di+3, di+4]: #in theory, this should give the track points on this date at 06, 12 ,18, 24 (00 of next date as taking the accum of the previous 6 hours)
				
								track_point_lat = tracklats[ts]
								track_point_lon = tracklons[ts]
						
								#circle = np.hypot(lons_-track_point_lon, lats_-track_point_lat)
								#circles.append(circles)
								print "shape circle: ",np.shape(np.hypot(lons_-track_point_lon, lats_-track_point_lat))
								print "shape pcp: ",np.shape(pcp)
								
								this_ts_pcp = np.ma.where(np.hypot(lons_-track_point_lon, lats_-track_point_lat) > 5, 0., pcp)
						
								for x in range(len(plats)):
									for y in range(len(plons)):
							
										if this_ts_pcp[x,y] > 0.:
									
											this_track_tc_pcp[di+4,x,y] = this_ts_pcp[x,y]
											
									
						
						
							#pcp_arr = np.ma.where(any(circles) > 5, 0., pcp)
							#this_track_tc_pcp[di+4] = pcp_arr
			
			
			
			
				
			outfile = Dataset(savedir+"GPM-IMERG_"+track_type+"_"+str(y1)+"_"+str(y2)+"_"+dir+"_TC-related-PRECIP.nc", 'w', format='NETCDF4')
			outfile.description = "GPM-IMERG TC-related precip (precip within a 5' radius of a TC) for "+track_type+" track. 24-hour accumulated precip around each 6-hourly observed track point in the past 24 hours. Precip should be saved at each 00 timestep, and this refers to the precip over the previous 24 hours. Except the start and end dates, which may not start/end at 00Z"
				
			outfile.createDimension('latitude', len(plats))
			outfile.createDimension('longitude', len(plons))
			outfile.createDimension('time', len(track_datelist))
				
			latitude = outfile.createVariable('latitude', 'f4', 'latitude')
			longitude = outfile.createVariable('longitude', 'f4', 'longitude')
			time = outfile.createVariable('time','i4','time')
			tc_pcp = outfile.createVariable('tc_pcp', 'float32', ('time', 'latitude', 'longitude'))
			
			latitude[:] = plats
			longitude[:] = plons
			time[:] = track_datelist
			tc_pcp[:,:,:] = this_track_tc_pcp
				
			latitude.units = 'degrees_north'
			longitude.units = 'degrees_east'
			time.units = '6-hourly starting from '+first_track_date.strftime("%Y-%m-%d-%H")
			tc_pcp.units = 'm'
			tc_pcp.long_name = 'Accumulated 24-hour precip within 5 degree radius around each 6-hourly observed track point in the last 24 hours'
				
			outfile.close()
				
			print "file saved: ", savedir+"GPM-IMERG_"+track_type+"_"+str(y1)+"_"+str(y2)+"_"+dir+"_TC-related-PRECIP.nc"
							

		precip_data.close()	
			
		

 


