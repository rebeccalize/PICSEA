import numpy as np
from netCDF4 import Dataset

trackdir = "/perm/mo/more/picsea/reformatted_idai_kenneth_tracks/"
precipdir = "/vol/floods/more/rebecca_TC_do_not_delete/PICSEA_DATA/PRECIP/nc/"
savedir = "/scratch/mo/more/idai_kenneth_ecmwf_forecast_tc_related_precip/"

#function to decide whether a given coordinate is within a given degree radius of a given point
def isInside(circle_x,circle_y, radius, x, y):	
	if ((x - circle_x)**2) + ((y - circle_y)**2) <= (radius**2):
		return True
	else:
		return False


for d in range(0): #29
	day = str(d).zfill(2)
	print d
	
	for h in [0,12]:
		hr = str(h).zfill(2)
		
		date = "201903"+day+hr
		print date, "deterministic"
		
		
		######################################################################################################
		
		#open track file
		det_track_file=trackdir+date+"_idai_ECMWF_det_3way_matched_reformatted.txt"
		det_track=np.genfromtxt(det_track_file, dtype = float, skip_header=1, usecols=np.arange(0,11))
		
		#get lats and lons of the right track (file contains the analysis and ibtracs as well as the forecast)
		NP = len(det_track[:,0])
		tid = 3
		xyi = []
		for i in range(NP):
			if det_track[i,0] == float(tid):
				xyi.append(i)
				
		print xyi
		print len(xyi)
	
		lats=[]
		lons=[]
		for i in xyi:
			lons.append(det_track[i,7])
			lats.append(det_track[i,8])

		#open precip file (deterministic)
		det_precip_file = precipdir+"ECMWF_PRECIP_"+date+"_detrm.nc"
		det_precip = Dataset(det_precip_file,'r')
		plons = det_precip.variables['lon'][:]
		plats = det_precip.variables['lat'][:]
		det_TP = det_precip.variables['TP'][:]
		flt = det_precip.variables['time'][:]
		
		print np.shape(det_TP)
		
		#for each lead time where we have a track point
		for lt in range(len(xyi)): #the track won't always go out that far ahead, so want to get the precip for the lead times when we have a track...
			
			#print np.max(det_TP[lt,:,:])
			
			track_point_lat = lats[lt]
			track_point_lon = lons[lt]
			
			
			lons_, lats_ = np.meshgrid(plons, plats)
			pcp0_arr = np.ma.where(np.hypot(lons_-track_point_lon, lats_-track_point_lat) > 5, 0.,det_TP[lt,:,:])
			tcp = pcp0_arr
			
			#print "track point lat: ", track_point_lat
			#print "track point lon: ", track_point_lon
			
			#create array for the TC-related precip in same shape as full array
			#tcp = np.zeros((2560, 5136))
			
			#for each lat lon in the original file, see if the coordinate is within a 5 degree radius of the track point at this lead time
			#for x in range(2560): #2560
				#for y in range(5136): #5136
			
					#tested this by printing out every incidence of "true" and it works fine!
					#insideradius = isInside(track_point_lat,track_point_lon,5.,plats[x],plons[y])
						
					#if the point is within the 5 degree radius, store the precip for that location
					#if insideradius == True:
						#tcp[x,y] = det_TP[lt,x,y]
						#print track_point_lat
						#print plats[x]
						#print track_point_lon
						#print plons[y]
						
					#if not, set the new array to nan outside the 5 degree radius 
					#else:
						#tcp[x,y] = np.nan
						
			#print np.nanmax(tc_precip)
			
			#outfile = savedir+"idai_tc_related_precip_ecmwf_det_"+date+"_lt_"+str(int(flt[lt])).zfill(3)+".txt"
			#np.savetxt(outfile, tcp)
			
			outfile = Dataset(savedir+"idai_tc_related_precip_ecmwf_det_"+date+"_lt_"+str(int(flt[lt])).zfill(3)+"_v2.nc",'w')
			src = Dataset(det_precip_file)
			
			for name in src.ncattrs():
				outfile.setncattr(name, src.getncattr(name))
			#excludedim=['time']
			for name, dimension in src.dimensions.iteritems():
				#if name not in excludedim:
				outfile.createDimension(name, (len(dimension) if not dimension.isunlimited else None))
			excludevar=['TP']
			for name, variable in src.variables.iteritems():
				if name not in excludevar:
					x = outfile.createVariable(name, variable.datatype, variable.dimensions)
					outfile.variables[name][:] = src.variables[name][:]
					
			tc_precip = outfile.createVariable('tc_precip','float32',('lat','lon'))
			tc_precip.setncattr('units','m')
			tc_precip[:] = tcp
			
			print "file saved: ", savedir+"idai_tc_related_precip_ecmwf_det_"+date+"_lt_"+str(int(flt[lt])).zfill(3)+"_v2.nc"	
			
			outfile.close()	
			src.close()
			


for d in range(1,20): #20
	day = str(d).zfill(2)
	print d
	
	for h in [0,12]:
		hr = str(h).zfill(2)
		
		date = "201903"+day+hr
		print date, "ensemble mean"
			
		#######################################################################################################	
			
		ensmean_track_file=trackdir+date+"_idai_ECMWF_mean_3way_matched_reformatted.txt"
		ensmean_track=np.genfromtxt(ensmean_track_file, dtype = float, skip_header=1)
		
		#get lats and lons of the right track (file contains the analysis and ibtracs as well as the forecast)
		#NP = len(ensmean_track[:,0])
		#tid = 3
		#xyi = []
		#for i in range(NP):
			#if det_track[i,0] == float(tid):
				#xyi.append(i)
				
		#print xyi
		#print len(xyi)
	
		ensmeanlats=ensmean_track[:,8]
		ensmeanlons=ensmean_track[:,7]
		#for i in xyi:
			#lons.append()
			#lats.append(det_track[i,8])

		#open precip file (deterministic)
		ensmean_precip_file = precipdir+"ECMWF_PRECIP_"+date+"_ensmean.nc"
		ensmean_precip = Dataset(ensmean_precip_file,'r')
		ensmeanplons = ensmean_precip.variables['lon'][:]
		ensmeanplats = ensmean_precip.variables['lat'][:]
		ensmean_TP = ensmean_precip.variables['TP'][:]
		ensmeanflt = ensmean_precip.variables['time'][:]
		
		print np.shape(ensmean_TP)
		
		#for each lead time where we have a track point
		for lt in range(len(ensmean_track)): #the track won't always go out that far ahead, so want to get the precip for the lead times when we have a track...
			
			print len(ensmean_track)
			
			#print np.max(det_TP[lt,:,:])
			
			ensmeantrack_point_lat = ensmeanlats[lt]
			ensmeantrack_point_lon = ensmeanlons[lt]
			
			
			ensmeanlons_, ensmeanlats_ = np.meshgrid(ensmeanplons, ensmeanplats)
			ensmeanpcp0_arr = np.ma.where(np.hypot(ensmeanlons_-ensmeantrack_point_lon, ensmeanlats_-ensmeantrack_point_lat) > 5, 0.,ensmean_TP[lt,:,:])
			ensmeantcp = ensmeanpcp0_arr
			
			#print "track point lat: ", track_point_lat
			#print "track point lon: ", track_point_lon
			
			#create array for the TC-related precip in same shape as full array
			#tcp = np.zeros((2560, 5136))
			
			#for each lat lon in the original file, see if the coordinate is within a 5 degree radius of the track point at this lead time
			#for x in range(2560): #2560
				#for y in range(5136): #5136
			
					#tested this by printing out every incidence of "true" and it works fine!
					#insideradius = isInside(track_point_lat,track_point_lon,5.,plats[x],plons[y])
						
					#if the point is within the 5 degree radius, store the precip for that location
					#if insideradius == True:
						#tcp[x,y] = det_TP[lt,x,y]
						#print track_point_lat
						#print plats[x]
						#print track_point_lon
						#print plons[y]
						
					#if not, set the new array to nan outside the 5 degree radius 
					#else:
						#tcp[x,y] = np.nan
						
			#print np.nanmax(tc_precip)
			
			#outfile = savedir+"idai_tc_related_precip_ecmwf_det_"+date+"_lt_"+str(int(flt[lt])).zfill(3)+".txt"
			#np.savetxt(outfile, tcp)
			
			ensmeanoutfile = Dataset(savedir+"idai_tc_related_precip_ecmwf_ensmean_"+date+"_lt_"+str(int(ensmeanflt[lt])).zfill(3)+".nc",'w')
			ensmeansrc = Dataset(ensmean_precip_file)
			
			for name in ensmeansrc.ncattrs():
				ensmeanoutfile.setncattr(name, ensmeansrc.getncattr(name))
			#excludedim=['time']
			for name, dimension in ensmeansrc.dimensions.iteritems():
				#if name not in excludedim:
				ensmeanoutfile.createDimension(name, (len(dimension) if not dimension.isunlimited else None))
			excludevar=['TP']
			for name, variable in ensmeansrc.variables.iteritems():
				if name not in excludevar:
					x = ensmeanoutfile.createVariable(name, variable.datatype, variable.dimensions)
					ensmeanoutfile.variables[name][:] = ensmeansrc.variables[name][:]
					
			ensmeantc_precip = ensmeanoutfile.createVariable('tc_precip','float32',('lat','lon'))
			ensmeantc_precip.setncattr('units','m')
			ensmeantc_precip[:] = ensmeantcp
			
			print "file saved: ", savedir+"idai_tc_related_precip_ecmwf_ensmean_"+date+"_lt_"+str(int(ensmeanflt[lt])).zfill(3)+".nc"	
			
			ensmeanoutfile.close()	
			ensmeansrc.close()
			



for d in range(0): #29
	day = str(d).zfill(2)
	print d
	
	for h in [0,12]:
		hr = str(h).zfill(2)
		
		date = "201903"+day+hr
		print date, "control"
			
		
		#######################################################################################################
		
		ctrl_track_file=trackdir+date+"_idai_ECMWF_eps_3way_matched_reformatted.txt"
		ctrl_track=np.genfromtxt(ctrl_track_file, dtype = float, skip_header=1, usecols=np.arange(0,11))
		
		#get lats and lons of the right track (file contains the analysis and ibtracs as well as the forecast)
		NP = len(ctrl_track[:,0])
		tid = 3
		xyi = []
		for i in range(NP):
			if ctrl_track[i,0] == float(tid):
				xyi.append(i)
				
		print xyi
		print len(xyi)
	
		ctrllats=[]
		ctrllons=[]
		for i in xyi:
			ctrllons.append(ctrl_track[i,7])
			ctrllats.append(ctrl_track[i,8])

		#open precip file (deterministic)
		ctrl_precip_file = precipdir+"ECMWF_PRECIP_"+date+"_ctrl.nc"
		ctrl_precip = Dataset(ctrl_precip_file,'r')
		ctrlplons = ctrl_precip.variables['lon'][:]
		ctrlplats = ctrl_precip.variables['lat'][:]
		ctrl_TP = ctrl_precip.variables['TP'][:]
		ctrlflt = ctrl_precip.variables['time'][:]
		
		print np.shape(ctrl_TP)
		
		#for each lead time where we have a track point
		for lt in range(len(xyi)): #the track won't always go out that far ahead, so want to get the precip for the lead times when we have a track...
			
			#print np.max(det_TP[lt,:,:])
			
			ctrltrack_point_lat = ctrllats[lt]
			ctrltrack_point_lon = ctrllons[lt]
			
			
			ctrllons_, ctrllats_ = np.meshgrid(ctrlplons, ctrlplats)
			ctrlpcp0_arr = np.ma.where(np.hypot(ctrllons_-ctrltrack_point_lon, ctrllats_-ctrltrack_point_lat) > 5, 0.,ctrl_TP[lt,:,:])
			ctrltcp = ctrlpcp0_arr
			
			
			ctrloutfile = Dataset(savedir+"idai_tc_related_precip_ecmwf_ctrl_"+date+"_lt_"+str(int(ctrlflt[lt])).zfill(3)+".nc",'w')
			ctrlsrc = Dataset(ctrl_precip_file)
			
			for name in ctrlsrc.ncattrs():
				ctrloutfile.setncattr(name, ctrlsrc.getncattr(name))
			#excludedim=['time']
			for name, dimension in ctrlsrc.dimensions.iteritems():
				#if name not in excludedim:
				ctrloutfile.createDimension(name, (len(dimension) if not dimension.isunlimited else None))
			excludevar=['TP']
			for name, variable in ctrlsrc.variables.iteritems():
				if name not in excludevar:
					x = ctrloutfile.createVariable(name, variable.datatype, variable.dimensions)
					ctrloutfile.variables[name][:] = ctrlsrc.variables[name][:]
					
			ctrltc_precip = ctrloutfile.createVariable('tc_precip','float32',('lat','lon'))
			ctrltc_precip.setncattr('units','m')
			ctrltc_precip[:] = ctrltcp
			
			print "file saved: ", savedir+"idai_tc_related_precip_ecmwf_ctrl_"+date+"_lt_"+str(int(ctrlflt[lt])).zfill(3)+".nc"	
			
			ctrloutfile.close()	
			ctrlsrc.close()
			
			
		

 


