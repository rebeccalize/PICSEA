import picsea_library as pl
import datetime
import iris
import itertools
import numpy as np
import os
import warnings
from netCDF4 import Dataset

MODEL_DIR = '/gws/nopw/j04/klingaman/fascinate/ukmo_nwp/precip/%Y/%Y%m%d%H'
#MODEL_DIR_2 = '/gws/nopw/j04/klingaman/emerton/ukmo_nwp_precip_data_notimebounds_variable/'
#MODEL_DIR_2 = '/gws/nopw/j04/klingaman/emerton/ukmo_nwp_data_nobounds_variables/'

MODEL_FILE = 'prods_op_gl-mn_%Y%m%d_%H_???.nc'
MODEL_FILE_2 = 'prods_op_gl-mn_%Y%m%d_%H_???.calc.nc'
MODEL_FILE_3 = 'prods_op_gl-mn_%Y%m%d_%H_006.calc.nc'
MODEL_FILE_4 = 'prods_op_gl-mn_%Y%m%d_%H_006.nc'
#before some date in 2007, the data looks different (stored in one file rather than split into many)
SINGLE_MODEL_FILE = 'prods_op_gl-mn_%Y%m%d_%H.nc'

COMP_PCP_TOT_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_precip_composites/"
COMP_PCP_TOT_FILE = 'ukmo_nwp.comp_pcp_tc.%d_days.%d%02d.nc'

COMP_PCP_TOT_DIR_TRMM = "/gws/nopw/j04/klingaman/emerton/ibtracs_precip_composites/"
COMP_PCP_TOT_FILE_TRMM = 'ibtracs.comp_pcp_tc.%d%02d.nc'

COMP_PCP_TC_DIR_TRMM_ANALYSIS = "/gws/nopw/j04/klingaman/emerton/analysis_trmm_tc_precip_composites/"
COMP_PCP_TC_FILE_TRMM_ANALYSIS = 'analysis.trmm.comp_pcp_tc.%d%02d.nc'

#This script contains functions to compute various different composites using the UKMO NWP forecasts, and TC tracks (that have been tracked and matched already)
#The functions are primarily from Simon's code from the fascinate project; script ukmo_nwp.py, and also use some miscellaneous functions from picsea_library.py

#To use these functions, the track output files should need to be rewritten and interpolated first, into a python-usable format (rewrite_track_files) interpolated to halfway between the existing times
#Also, before using these, a list of dates to exclude at a given lead time is needed, due to missing NWP data (use exclude_days function)

#There are functions to compute the following:

#Compute composites of precip due to:
#   All TCs at a particular lead time, for a particular year and month (composite_pcp_tc_year_month)
#   All TCs in the forecast initialised at a given date/time (single forecast) (composite_pcp_tc_single_forecast)
#   All TCs in forecasts for each MJO phase (composite_pcp_tc_forecasts_mjo)
#   TC-related precip in all forecast windows (??) for the given year and month (composite_pcp_tc_forecast_windows_month)
#   TC-related precip in all forecast windows, using file already generated for each year and month above (composite_pcp_tc_forecast_windows_all)
#   TCs at a particular forecast lead time, for a particular year (composite_pcp_tc_year)
#   As above, but for a year in which the resolution changes midway through (composite_pcp_tc_year_separate_resolutions)
#   TCs at a particular forecast lead time, for a particular month in the given years (composite_pcp_tc_month)
#   As above, but taking into account resolution changes (composite_pcp_tc_month_separate_resolutions)
#   
#   Note: Computing the composites requires the functions that (1) produce a list of excluded dates (2) calculate accumulated precip (3) add cubes (see common.py in simon's code)

def trmm_3b42_pcp_accumulation(year, month, day, hour, lon=None, lat=None,
							   quiet=False):
	"""Returns a map of 6-hour accumulated precipitation in mm from TRMM 3B42,
	over the 6-hour period CENTRED ON the given time in UTC.  Hour must be one
	of [3, 9, 15, 21].

	Values more than 5 degrees from the given lat/lon location set to zero (if
	lat/lon are given).

	3B42 is provided every 3 hours as a mean over the 3-hour period CENTERED on
	the given time.  Hence, if hour=9 then we sum:
		  TRMM 06UTC file x1.5
		+ TRMM 09UTC file x3.0
		+ TRMM 12UTC file x1.5
	"""

	# Check we have a valid time (same as track points)
	valid_hours = range(3, 22, 6)
	if hour not in valid_hours:
		raise ValueError('Invalid hour %s; must be one of %s' % (hour, valid_hours))

	# Iterate for each time
	t = datetime.datetime(year, month, day, hour)
	dt = datetime.timedelta(hours=3)
	pcp_list = iris.cube.CubeList()
	for mult, time in itertools.izip([1.5, 3., 1.5], [t - dt, t, t + dt]):
		if not quiet:
			print '   ', time.strftime('%Y/%m/%d %H:%M'), 'x%.01f' % mult

		# Open TRMM 3B42 file for given time
		pcp = pl.read_trmm_3b42(time.year, time.month, time.day, time.hour)

		##################################################################
		# Set to zero outside 5 degree radius
		# Tracks use longitude values 0 -> 360 but TRMM uses -180 -> 180
		##################################################################
		if lat is not None and lon is not None:
			lats = pcp.coord(axis='Y').points
			lons1 = pcp.coord(axis='X').points
			lons = np.where(lons1 < 0, lons1 + 360., lons1)
			lons_, lats_ = np.meshgrid(lons, lats)
			pcp.data = np.ma.where(np.hypot(lons_ - lon, lats_ - lat) > 5, 0.,
								   pcp.data * mult)
		else:
			pcp.data = pcp.data * mult
		pcp_list.append(pcp)

	# Return sum
	tot_pcp = pl.add_cubes(pcp_list, contributing_days=False)
	tot_pcp.units = 'mm'
	tot_pcp.standard_name = 'lwe_thickness_of_precipitation_amount'
	tot_pcp.long_name = 'accumulated precipitation'
	tot_pcp.var_name = 'pcp'
	return tot_pcp

def nwp_pcp_accumulation(forecast_date, forecast_hours, lon=None, lat=None):
	"""Returns a map of 6-hour accumulated precipitation in mm from the
	forecast initialized at *forecast_date*, over the 6-hour period CENTRED ON
	*forecast_hours* hours after the initialization.  Hour of validity time must
	be one of [3, 9, 15, 21].

	Values more than 5 degrees from the given lat/lon location set to zero (if
	lat/lon are given).

	Forecast accumulations are given a timestamp at the END of the 6-hour
	period, so this function adds 3 hours to get the appropriate time co-
	ordinate value.

	Arguments:

	*forecast_date*
		`datetime.datetime`

	*forecast_hours*
		`integer`

	*lat*, *lon*=None
		`float`s (if None, returns entire field without setting any values to
		zero)
	"""

	# Work out required time
	gt = forecast_date + datetime.timedelta(hours=forecast_hours)

	# Check we have a valid time (same as track points)
	valid_hours = range(3, 22, 6)
	if gt.hour not in valid_hours:
		raise ValueError('Invalid hour %s for forecast_date+forecast_hours; '
						 'must be one of %s' % (gt.hour, valid_hours))

	# Callback function
	def cb(cube, field, filename):
		for att in ['runid', 'history', 'NCO']:
			try:
				del cube.attributes[att]
			except KeyError:
				pass

	# Open model output file for the given forecast initialization time
	if forecast_date >= datetime.datetime(2007, 11, 28, 0) and forecast_date <= datetime.datetime(2017,12,17,0):
		model_path = forecast_date.strftime(os.path.join(MODEL_DIR, MODEL_FILE))
		print model_path
		
	elif forecast_date >= datetime.datetime(2017,12,17,12) and forecast_date <= datetime.datetime(2018,02,13,0):
		model_path = forecast_date.strftime(os.path.join(MODEL_DIR, MODEL_FILE_2))
		
	elif forecast_date >= datetime.datetime(2018,02,13,12) and forecast_date <= datetime.datetime(2018,9,25,0):
		#print "DATE IS PAST 13 FEB 2018!"
		#to get the cubes, as normal:
		model_path = forecast_date.strftime(os.path.join(MODEL_DIR, MODEL_FILE_2))
		
		#to get just one file for this forecast date, to get the correct bounds value and create the time_points array manually		
		infile_weird_data = forecast_date.strftime(os.path.join(MODEL_DIR, MODEL_FILE_3))
		 
		#the netcdf files past this date have a strange time_bounds variable, with 2 values. we want the second, but the cubes pick up the first
		#the second value matches with all the files from previous dates, following the same conventions
		#so get this second value, and then create an array of the 28 values (the bounds), going up by 0.25 monotonically... 
		#later, replace the weird bounds values in the cubes, with the correct values that we actually want to use and fit with the norm		
		ff_weird_data = Dataset(infile_weird_data,'r')
		bounds_weird_data = ff_weird_data.variables['time_bounds'][:]
		time_points_start = bounds_weird_data[1]
		print time_points_start
		
		manual_time_points = np.zeros(28)
		manual_time_points[0] = time_points_start
		for i in range(1,28):
			manual_time_points[i] += manual_time_points[i-1]+0.25
		print manual_time_points
		
	elif forecast_date >= datetime.datetime(2018,9,25,12):
		model_path = forecast_date.strftime(os.path.join(MODEL_DIR, MODEL_FILE))
		
			#to get just one file for this forecast date, to get the correct bounds value and create the time_points array manually		
		infile_weird_data = forecast_date.strftime(os.path.join(MODEL_DIR, MODEL_FILE_4))
		 
		#the netcdf files past this date have a strange time_bounds variable, with 2 values. we want the second, but the cubes pick up the first
		#the second value matches with all the files from previous dates, following the same conventions
		#so get this second value, and then create an array of the 28 values (the bounds), going up by 0.25 monotonically... 
		#later, replace the weird bounds values in the cubes, with the correct values that we actually want to use and fit with the norm		
		ff_weird_data = Dataset(infile_weird_data,'r')
		bounds_weird_data = ff_weird_data.variables['time_bounds'][:]
		time_points_start = bounds_weird_data[1]
		print time_points_start
		
		manual_time_points = np.zeros(28)
		manual_time_points[0] = time_points_start
		for i in range(1,28):
			manual_time_points[i] += manual_time_points[i-1]+0.25
		print manual_time_points
		
		
	else:
		model_path = forecast_date.strftime(os.path.join(MODEL_DIR, SINGLE_MODEL_FILE))
		
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		if forecast_date < datetime.datetime(2018,9,25,12):
			cubes = iris.cube.CubeList([c for c in iris.load(model_path, callback=cb) if c.var_name.startswith('UM_m10s20i013_vn') and not c.var_name.endswith('_1')])
		elif forecast_date >= datetime.datetime(2018,9,25,12):
			cubes = iris.cube.CubeList([c for c in iris.load(model_path, callback=cb) if c.var_name.startswith('UM_m01s05i216_vn') and not c.var_name.endswith('_1')])
			
		#print "ORIGINAL CUBES FROM FILELIST:"
		#print cubes
	if not len(cubes):
		print model_path
		raise IOError('No precip data for %s T+%03d' % (forecast_date.strftime('%Y%m%d %H:00'), forecast_hours))

	
	#From 2018021312, the cubes are loaded in with an anonymous coordinate (--: 28), that should be the time dim coord, and usually this is "time: 28", 
	#with a time dim coord, and no time aux coord. however, we now get the anonymous coordinate, no time dim coord, and a time aux coord
	#how to fix this? ideally, promote the time coord to a dim coord, but this gives the 'the bounds array must be strictly monotonic' error :( 
	# Some Cubes have scalar time coord; promote to dim coord
	#for ii, icube in enumerate(cubes):
		#if len(icube.dim_coords[:]) < 3:
		
			#print icube
			#icube.coord('time').bounds = None
			
			#This gives: ValueError: Unequal lengths. Cube dimension 0 => 1; coord u'time' => 28.
			#cubes[ii] = iris.util.new_axis(icube, 'time')
			
			#####
			#this works, without any of the "new_axis" stuff BUT 
			#then when concatenating later, get "AttributeError: 'NoneType' object has no attribute 'metadata'"
			#cubes[ii] = iris.util.promote_aux_coord_to_dim_coord(icube, 'time')
			#print icube
			#####
			
			#Tried the following from https://groups.google.com/forum/#!topic/scitools-iris/FsT14rFtumg
			#But gives the error "ValueError: The dim_coord may not be an AuxCoord instance."
			#icube.coord('time').bounds = None
			#newCoord = icube.coord('time')			
			#icube.remove_coord('time')
			#icube.add_dim_coord(newCoord,0)
			
			#Tried the following, but the anonymous coord isn't actually 0, so this just renames the latitude coordinate
			#coord = icube.dim_coords[0]
			#coord.rename('time')
			#print icube
	
	#this next bit was in Simon's code, but doesn't seem to ever be used
	
	# Get rid of duplicate times (at the moment, works only if they appeared in
	# a Cube with scalar time dimension)
	#for ii, icube in enumerate(cubes):
		#if len(icube.coord(axis='T').points) == 1:
			#for jcube in (cubes[:ii] + cubes[ii+1:]):
				#if icube.coord(axis='T').points[0] in jcube.coord(axis='T').points:
					#cubes.remove(icube)
					#break
	
	
	################################################################################################################
	#This loop works, AND then the concatenation works, and all the info printed looks the way it normally does when the code runs,
	#BUT the constraint function then doesn't work (I think it's that bit), and the resulting pcp cube is 'None'
	#so also need to fix the time_points values, below
	
	for cube in cubes:
		if len(cube.dim_coords[:]) < 3:
			print cube
			print cube.coord('time').points
			cube.coord('time').bounds = None
			iris.util.promote_aux_coord_to_dim_coord(cube,'time')
			print cube
			print cube.coord('time').points
			
	#the time_bounds are weird in the files past 13th feb 2018
	#the earlier files, have one bounds value, which is always a multiple of .25 (e.g. 44.25, 44.5, 44.75, 45), and increases by .25 for each new forecast date
	#but in the later files, the bounds have 2 values, and the second one is the one we want (matches previos files), but the script picks up the first one
	#so when we read in the files earlier, we also read in the bounds variable, chose the correct value, and created our own time points array
	#so here, we have to replace the time points with the correct ones, and then everything seems to work... 
	
	cube = cubes.concatenate_cube()
	if forecast_date >= datetime.datetime(2018,02,13,12):
		cube.coord('time').points = manual_time_points

	print "CONCATENATED CUBE:"
	print cube
	print cube.coord('time').points
	###################################################################################################################	
	
		
	# Get correct time
	tcoord = cube.coord(axis='T')
	#print tcoord
	con_func = lambda tt: tt.point == (gt + datetime.timedelta(hours=3))
	con = iris.Constraint(**{tcoord.name(): con_func})
	#with iris.FUTURE.context(cell_datetime_objects=True):
	     #pcp = cube.extract(con)
	pcp = cube.extract(con)
	print "PCP CUBE:"
	print pcp
	if pcp is None:
		raise IOError('No precip data for %s T+%03d' % (forecast_date.strftime('%Y%m%d %H:00'), forecast_hours))

	# Set precip to zero outside 5 degree radius of TC track
	if lat is not None and lon is not None:
		pcp = iris.util.squeeze(pcp)
		lats = pcp.coord(axis='Y').points
		lons = pcp.coord(axis='X').points
		lons_, lats_ = np.meshgrid(lons, lats)
		pcp0_arr = np.ma.where(np.hypot(lons_-lon, lats_-lat) > 5, 0.,
							   pcp.data)
		pcp.data = pcp0_arr

	# Set metadata and return
	pcp.units = 'mm'
	pcp.standard_name = 'lwe_thickness_of_precipitation_amount'
	pcp.long_name = 'accumulated precipitation'
	pcp.var_name = 'pcp'
	return pcp

def composite_pcp_tot_month_separate_resolutions(month, lead):
	"""Computes the total precipitation at a particular forecast lead time (in
	days) for a particular month in the given years."""

	# Iterate for each resolution
	#y1 = [2006, 2008][lead==6]
	#y2 = [2016, 2017][month <= 7]
	#y2=2016
	#res_list = {2006: ['n320', 'n512', 'n768'], 2008: ['n320']}[y1]
	res_list = ['n320', 'n512', 'n768','n1280']
	#res_list = ['n320', 'n768']
	for res in res_list:
		print res
		# Iterate for each year
		if res == 'n320':
			if lead == 6:
				start_year = 2009
			elif month <= 6:
				start_year = 2007
			else:
				start_year = 2006
			if month <= 3:
				years = range(start_year, 2011)
			else:
				years = range(start_year, 2010)
		elif res == 'n512':
			start_year = 2006
			if month <= 2:
				years = range(2011, 2015)
			elif month <= 7:
				years = range(2010, 2015)
			else:
				years = range(2010, 2014)
		elif res == 'n768':
			start_year = [2006, 2008][lead==6]

			if month <= 6:
				years = range(2015, 2018)
			elif month == 7:
				years = range(2014, 2018)
			else:
				years = range(2014, 2017)
				
				
		elif res == 'n1280':
			start_year = [2006, 2008][lead==6]
			
			end_year=2018 # THIS WILL BE 2019 ONCE HAVE THE 2018-2019 IBTRACS ETC (FULL DATASET_
			
			if month <= 6:
				years=range(2018,end_year+1)
			elif month == 7:
				years=range(2018,end_year+1)
			else:
				years=range(2017,2018)
				
			


		print years
		infiles = []
		for y in years:
			tot_file = COMP_PCP_TOT_FILE % (lead, y, month)
			if (y, month) in [(2010, 3), (2014, 7), (2017,7)]:
				tot_file = tot_file.replace('.nc', '.%s.nc' % res)
			#print tot_file
			infiles.append(os.path.join(COMP_PCP_TOT_DIR, str(lead), str(y), tot_file))

		print infiles

		outfile = infiles[0].replace('/%d' % years[0], '')
		outfile = outfile.replace('%d%02d' % (years[0], month), '%d_%d.%02d' % (years[0], years[-1], month))
		if res not in outfile:
			outfile = outfile.replace('.nc', '.%s.nc' % res)
		print outfile
		pl.add_files(infiles, outfile, deal_with_masks=True, cube_count=len(years))

def composite_pcp_tot_month_separate_resolutions_trmm(month):
	"""Computes the total precipitation at a particular forecast lead time (in
	days) for a particular month in the given years."""

	# Iterate for each resolution
	#y1 = [2006, 2008][lead==6]
	#y2 = [2016, 2017][month <= 7]
	#y2=2016
	#res_list = {2006: ['n320', 'n512', 'n768'], 2008: ['n320']}[y1]
	res_list = ['n320', 'n512', 'n768','n1280']
	#res_list = ['n320', 'n768']
	for res in res_list:
		print res
		# Iterate for each year
		if res == 'n320':
			if month <= 6:
				start_year = 2007
			else:
				start_year = 2006
			if month <= 3:
				years = range(start_year, 2011)
			else:
				years = range(start_year, 2010)
		elif res == 'n512':
			#start_year = 2006
			if month <= 2:
				years = range(2011, 2015)
			elif month <= 7:
				years = range(2010, 2015)
			else:
				years = range(2010, 2014)
				
		elif res == 'n768':

				
			if month <= 6:
				years = range(2015, 2018)
			elif month == 7:
				years = range(2014, 2018)
			else:
				years = range(2014, 2017)
				
				
		elif res == 'n1280':
			#start_year = [2006, 2008][lead==6]
			
			end_year=2018 # THIS WILL BE 2019 ONCE HAVE THE 2018-2019 IBTRACS ETC (FULL DATASET)
			
			if month <= 6:
				years=range(2018,end_year+1)
			elif month == 7:
				years=range(2018,end_year+1)
			else:
				years=range(2017,2018)


		print years
		infiles = []
		for y in years:
			tot_file = COMP_PCP_TOT_FILE_TRMM % (y, month)
			if (y, month) in [(2010, 3), (2014, 7), (2017,7)]:
				tot_file = tot_file.replace('.nc', '.%s.nc' % res)
			#print tot_file
			infiles.append(os.path.join(COMP_PCP_TOT_DIR_TRMM, str(y), tot_file))

		print infiles

		outfile = infiles[0].replace('/%d' % years[0], '')
		outfile = outfile.replace('%d%02d' % (years[0], month), '%d_%d.%02d' % (years[0], years[-1], month))
		if res not in outfile:
			outfile = outfile.replace('.nc', '.%s.nc' % res)
		print outfile
		pl.add_files(infiles, outfile, deal_with_masks=True, cube_count=len(years))
		
def composite_pcp_tc_month_separate_resolutions_trmm_analysis(month):
	"""Computes the total precipitation at a particular forecast lead time (in
	days) for a particular month in the given years."""

	# Iterate for each resolution
	#y1 = [2006, 2008][lead==6]
	#y2 = [2016, 2017][month <= 7]
	#y2=2016
	#res_list = {2006: ['n320', 'n512', 'n768'], 2008: ['n320']}[y1]
	res_list = ['n320', 'n512', 'n768'] #,'n1280'
	#res_list = ['n320', 'n768']
	for res in res_list:
		print res
		# Iterate for each year
		if res == 'n320':
			if month <= 6:
				start_year = 2007
			else:
				start_year = 2006
			if month <= 3:
				years = range(start_year, 2011)
			else:
				years = range(start_year, 2010)
		elif res == 'n512':
			#start_year = 2006
			if month <= 2:
				years = range(2011, 2015)
			elif month <= 7:
				years = range(2010, 2015)
			else:
				years = range(2010, 2014)
				
		elif res == 'n768':

			#when precip composites are updates for ukmo nwp, the numbers in brackets should be the ones commented out after, rather than ending in 2016
			if month <= 6:
				years = range(2015, 2016) #2015, 2018
			elif month == 7:
				years = range(2014, 2016) #2014, 2018
			else:
				years = range(2014, 2016) #2014, 2017
				
				
		#elif res == 'n1280':
			###start_year = [2006, 2008][lead==6]
			
			#end_year=2018 # THIS WILL BE 2019 ONCE HAVE THE 2018-2019 IBTRACS ETC (FULL DATASET)
			
			#if month <= 6:
				#years=range(2018,end_year+1)
			#elif month == 7:
				#years=range(2018,end_year+1)
			#else:
				#years=range(2017,2018)


		print years
		infiles = []
		for y in years:
			tot_file = COMP_PCP_TC_FILE_TRMM_ANALYSIS % (y, month)
			if (y, month) in [(2010, 3), (2014, 7), (2017,7)]:
				tot_file = tot_file.replace('.nc', '.%s.nc' % res)
			#print tot_file
			infiles.append(os.path.join(COMP_PCP_TC_DIR_TRMM_ANALYSIS, str(y), tot_file))

		print infiles

		outfile = infiles[0].replace('/%d' % years[0], '')
		outfile = outfile.replace('%d%02d' % (years[0], month), '%d_%d.%02d' % (years[0], years[-1], month))
		if res not in outfile:
			outfile = outfile.replace('.nc', '.%s.nc' % res)
		print outfile
		pl.add_files(infiles, outfile, deal_with_masks=True, cube_count=len(years))

def composite_pcp_tot_all_separate_resolutions(lead):
	"""Computes the total precipitation at a particular forecast lead time (in
	days) in the given range of years."""
	"""Requires monthly totals first"""

	# Iterate for each resolution
	for res in ['n320', 'n512', 'n768','n1280']:
	#for res in ['n768']:

		if res == 'n320':
			y1 = '072006'
			y2 = '032010'
		elif res == 'n512':
			y1 = '032010'
			y2 = '072014'
		elif res == 'n768':
			y1 = '072014'
			y2 = '072017'
		elif res == 'n1280':
			y1 = '072017'
			y2 = '062018' #EVENTUALLY WILL BE 062019 
		

		infile = []

		for m in xrange(1, 13):

			if res == 'n320':
				if lead == 6:
					start_year = 2009
				elif m <= 6:
					start_year = 2007
				else:
					start_year = 2006
				if m <= 3:
					years = range(start_year, 2011)
				else:
					years = range(start_year, 2010)

			elif res == 'n512':
				#start_year = 2006
				if m <= 2:
					years = range(2011, 2015)
				elif m <= 7:
					years = range(2010, 2015)
				else:
					years = range(2010, 2014)

			elif res == 'n768':
					
				if month <= 6:
					years = range(2015, 2018)
				elif month == 7:
					years = range(2014, 2018)
				else:
					years = range(2014, 2017)
				
				
			elif res == 'n1280':
			
				end_year=2018 # THIS WILL BE 2019 ONCE HAVE THE 2018-2019 IBTRACS ETC (FULL DATASET_
			
				if month <= 6:
					years=range(2018,end_year+1)
				elif month == 7:
					years=range(2018,end_year+1)
				else:
					years=range(2017,2018) #THIS IS CURRENTLY 2017-2017 SO THAT WE JUST HAVE 072017-062018; EVENTUALLY THIS WILL ALSO BE "end_year+1" as well
					
					

			if lead is None:
				tot_file = COMP_PCP_TOT_FILE.replace('%d%02d', '%d_%d.%02d')
				tot_file = tot_file % (years[0], years[-1], m)
			else:
				tot_file = COMP_PCP_TOT_FILE.replace('%d%02d', '%d_%d.%02d')
				tot_file = tot_file % (lead, years[0], years[-1], m)
			tot_file = tot_file.replace('.nc', '.%s.nc' % res)
			infile.append(os.path.join(COMP_PCP_TOT_DIR, [str(lead), 'all'][lead is None], tot_file))

		print infile
		outfilename= 'ukmo_nwp.comp_pcp_tc.%d_days.%s-%s.%s.nc' % (lead,y1 ,y2 , res)
		outfile = os.path.join(COMP_PCP_TOT_DIR, str(lead), outfilename)
		pl.add_files(infile, outfile, deal_with_masks=True, cube_count=12)
		print '################'
		print outfile
		print '################'


def composite_pcp_tot_all_separate_resolutions_trmm():
	"""Computes the total precipitation at a particular forecast lead time (in
	days) in the given range of years."""
	"""Requires monthly totals first"""

	# Iterate for each resolution
	for res in ['n320', 'n512', 'n768','n1280']:
	#for res in ['n768']:

		if res == 'n320':
			y1 = '072006'
			y2 = '032010'
		elif res == 'n512':
			y1 = '032010'
			y2 = '072014'
		elif res == 'n768':
			y1 = '072014'
			y2 = '072017'
		elif res == 'n1280':
			y1 = '072017'
			y2 = '062018'

		infile = []

		for m in xrange(1, 13):

			if res == 'n320':
				#if lead == 6:
					#start_year = 2009
				if m <= 6:
					start_year = 2007
				else:
					start_year = 2006
				if m <= 3:
					years = range(start_year, 2011)
				else:
					years = range(start_year, 2010)

			elif res == 'n512':
				#start_year = 2006
				if m <= 2:
					years = range(2011, 2015)
				elif m <= 7:
					years = range(2010, 2015)
				else:
					years = range(2010, 2014)

			elif res == 'n768':
					
				if month <= 6:
					years = range(2015, 2018)
				elif month == 7:
					years = range(2014, 2018)
				else:
					years = range(2014, 2017)
				
				
			elif res == 'n1280':
			
				end_year=2018 # THIS WILL BE 2019 ONCE HAVE THE 2018-2019 IBTRACS ETC (FULL DATASET_
			
				if month <= 6:
					years=range(2018,end_year+1)
				elif month == 7:
					years=range(2018,end_year+1)
				else:
					years=range(2017,2018) #THIS IS CURRENTLY 2017-2017 SO THAT WE JUST HAVE 072017-062018; EVENTUALLY THIS WILL ALSO BE "end_year+1" as well

			#if lead is None:
				#tot_file = COMP_PCP_TOT_FILE_TRMM.replace('%d%02d', '%d_%d.%02d')
				#tot_file = tot_file % (years[0], years[-1], m)
			#else:
			tot_file = COMP_PCP_TOT_FILE_TRMM.replace('%d%02d', '%d_%d.%02d')
			tot_file = tot_file % ( years[0], years[-1], m)
			tot_file = tot_file.replace('.nc', '.%s.nc' % res)
			infile.append(os.path.join(COMP_PCP_TOT_DIR_TRMM, tot_file))

		print infile
		outfilename= 'trmm.comp_pcp_tc.%s-%s.%s.nc' % (y1 ,y2 , res)
		outfile = os.path.join(COMP_PCP_TOT_DIR_TRMM, outfilename)
		pl.add_files(infile, outfile, deal_with_masks=True, cube_count=12)
		print '################'
		print outfile
		print '################'
		
def composite_pcp_tc_all_separate_resolutions_trmm_analysis():
	"""Computes the total precipitation at a particular forecast lead time (in
	days) in the given range of years."""
	"""Requires monthly totals first"""

	# Iterate for each resolution
	for res in ['n320','n512','n768','n1280']: #
	#for res in ['n768']:

		if res == 'n320':
			y1 = '072006'
			y2 = '032010'
		elif res == 'n512':
			y1 = '032010'
			y2 = '072014'
		elif res == 'n768':
			y1 = '072014'
			y2 = '072017'
		elif res == 'n1280':
			y1 = '072017'
			y2 = '062018'

		infile = []

		for m in xrange(1, 13):

			if res == 'n320':
				#if lead == 6:
					#start_year = 2009
				if m <= 6:
					start_year = 2007
				else:
					start_year = 2006
				if m <= 3:
					years = range(start_year, 2011)
				else:
					years = range(start_year, 2010)

			elif res == 'n512':
				#start_year = 2006
				if m <= 2:
					years = range(2011, 2015)
				elif m <= 7:
					years = range(2010, 2015)
				else:
					years = range(2010, 2014)

			elif res == 'n768':
					
				if m <= 6:
					years = range(2015, 2018) #2015,2018
				elif m == 7:
					years = range(2014, 2018) #2014,2018
				else:
					years = range(2014, 2017) #2014,2017
				
				
			elif res == 'n1280':
			
				end_year=2018 # THIS WILL BE 2019 ONCE HAVE THE 2018-2019 IBTRACS ETC (FULL DATASET_
			
				if m <= 6:
					years=range(2018,end_year+1)
				elif m == 7:
					years=range(2018,end_year+1)
				else:
					years=range(2017,2018) #THIS IS CURRENTLY 2017-2017 SO THAT WE JUST HAVE 072017-062018; EVENTUALLY THIS WILL ALSO BE "end_year+1" as well

			#if lead is None:
				#tot_file = COMP_PCP_TOT_FILE_TRMM.replace('%d%02d', '%d_%d.%02d')
				#tot_file = tot_file % (years[0], years[-1], m)
			#else:
			tot_file = COMP_PCP_TC_FILE_TRMM_ANALYSIS.replace('%d%02d', '%d_%d.%02d')
			tot_file = tot_file % ( years[0], years[-1], m)
			tot_file = tot_file.replace('.nc', '.%s.nc' % res)
			infile.append(os.path.join(COMP_PCP_TC_DIR_TRMM_ANALYSIS, tot_file))

		print infile
		outfilename= 'trmm.comp_pcp_tc.analysis.%s-%s.%s.nc' % (y1 ,y2 , res)
		outfile = os.path.join(COMP_PCP_TC_DIR_TRMM_ANALYSIS, outfilename)
		pl.add_files(infile, outfile, deal_with_masks=True, cube_count=12)
		print '################'
		print outfile
		print '################'

def composite_pcp_tot_seasons_res(lead):
	"""Computes the total precipitation at a particular forecast lead time (in
	days) in the given range of years for each season JJASON and DJFMAM."""

	#first_year = [2006, 2009][lead == 6]
	months_dict = {'NDJFMA': [11,12,1,2,3,4], 'MJJASO': [5,6,7,8,9,10], 'DJFM': [12,1,2,3]}
	for res in ['n320', 'n512', 'n768','n1280']:

		if res == 'n320':
			y1 = '072006'
			y2 = '032010'
		elif res == 'n512':
			y1 = '032010'
			y2 = '072014'
		elif res == 'n768':
			y1 = '072014'
			y2 = '072017'
		elif res == 'n1280':
			y1 = '072017'
			y2 = '062018'

		for season, months in months_dict.iteritems():
			if season in ['NDJFMA', 'MJJASO']:
				cube_count=6
			elif season == 'DJFM':
				cube_count=4
			infile = []
			for m in months:

				if res == 'n320':
					if lead == 6:
						start_year = 2009
					elif m <= 6:
						start_year = 2007
					else:
						start_year = 2006
					if m <= 3:
						years = range(start_year, 2011)
					else:
						years = range(start_year, 2010)

				elif res == 'n512':
					start_year = 2006
					if m <= 2:
						years = range(2011, 2015)
					elif m <= 7:
						years = range(2010, 2015)
					else:
						years = range(2010, 2014)

				elif res == 'n768':
					
					if month <= 6:
						years = range(2015, 2018)
					elif month == 7:
						years = range(2014, 2018)
					else:
						years = range(2014, 2017)
				
				
				elif res == 'n1280':
			
					end_year=2018 # THIS WILL BE 2019 ONCE HAVE THE 2018-2019 IBTRACS ETC (FULL DATASET_
			
					if month <= 6:
						years=range(2018,end_year+1)
					elif month == 7:
						years=range(2018,end_year+1)
					else:
						years=range(2017,2018) #THIS IS CURRENTLY 2017-2017 SO THAT WE JUST HAVE 072017-062018; EVENTUALLY THIS WILL ALSO BE "end_year+1" as well


				#if m <= 7:
					#last_year = 2017
				#else:
					#last_year = 2016
				tot_file = COMP_PCP_TOT_FILE.replace('%d%02d', '%d_%d.%02d.%s')
				tot_file = tot_file % (lead, years[0], years[-1], m, res)
				infile.append(os.path.join(COMP_PCP_TOT_DIR, str(lead),tot_file))

			outfilename = 'ukmo_nwp.comp_pcp_tc.%d_days.%s-%s.%s.%s.nc' % (lead, y1, y2, res, season)
			outfile = os.path.join(COMP_PCP_TOT_DIR, str(lead), outfilename)

			pl.add_files(infile, outfile, deal_with_masks=True, cube_count=cube_count)


def composite_pcp_tot_seasons_res_trmm():
	"""Computes the total precipitation at a particular forecast lead time (in
	days) in the given range of years for each season JJASON and DJFMAM."""

	#first_year = [2006, 2009][lead == 6]
	months_dict = {'NDJFMA': [11,12,1,2,3,4], 'MJJASO': [5,6,7,8,9,10], 'DJFM': [12,1,2,3]}
	for res in ['n320', 'n512', 'n768','n1280']:

		if res == 'n320':
			y1 = '072006'
			y2 = '032010'
		elif res == 'n512':
			y1 = '032010'
			y2 = '072014'
		elif res == 'n768':
			y1 = '072014'
			y2 = '072017'
		elif res == 'n1280':
			y1 = '072017'
			y2 = '062018' 

		for season, months in months_dict.iteritems():
			if season in ['NDJFMA', 'MJJASO']:
				cube_count=6
			elif season == 'DJFM':
				cube_count=4
			infile = []
			for m in months:

				if res == 'n320':
					#if lead == 6:
						#start_year = 2009
					if m <= 6:
						start_year = 2007
					else:
						start_year = 2006
					if m <= 3:
						years = range(start_year, 2011)
					else:
						years = range(start_year, 2010)

				elif res == 'n512':
					#start_year = 2006
					if m <= 2:
						years = range(2011, 2015)
					elif m <= 7:
						years = range(2010, 2015)
					else:
						years = range(2010, 2014)

				elif res == 'n768':
					
					if month <= 6:
						years = range(2015, 2018)
					elif month == 7:
						years = range(2014, 2018)
					else:
						years = range(2014, 2017)
				
				
				elif res == 'n1280':
			
					end_year=2018 # THIS WILL BE 2019 ONCE HAVE THE 2018-2019 IBTRACS ETC (FULL DATASET_
			
					if month <= 6:
						years=range(2018,end_year+1)
					elif month == 7:
						years=range(2018,end_year+1)
					else:
						years=range(2017,2018) #THIS IS CURRENTLY 2017-2017 SO THAT WE JUST HAVE 072017-062018; EVENTUALLY THIS WILL ALSO BE "end_year+1" as well


				#if m <= 7:
					#last_year = 2017
				#else:
					#last_year = 2016
				tot_file = COMP_PCP_TOT_FILE_TRMM.replace('%d%02d', '%d_%d.%02d.%s')
				tot_file = tot_file % (years[0], years[-1], m, res)
				infile.append(os.path.join(COMP_PCP_TOT_DIR_TRMM,tot_file))

			outfilename = 'trmm.comp_pcp_tc.%s-%s.%s.%s.nc' % (y1, y2, res, season)
			outfile = os.path.join(COMP_PCP_TOT_DIR_TRMM, outfilename)

			pl.add_files(infile, outfile, deal_with_masks=True, cube_count=cube_count)
