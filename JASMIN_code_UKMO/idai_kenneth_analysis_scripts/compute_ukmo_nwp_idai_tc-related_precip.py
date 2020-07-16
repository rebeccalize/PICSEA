import picsea_library as pl
import composite_functions_ukmo_nwp as cf
import datetime
import iris
import itertools
import numpy as np
import os
import warnings
import sys


lead=int(sys.argv[1])
lead=lead-1

#indata
TRACK_DIR_3WAY = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/IDAI/reformatted_idai_forecast_tracks/"
TRACK_FILE='idai_UKMO_%Y%m%d%H_interpolated.txt' #ukmo_nwp.2006111400.tcident.interp

#outdata
COMP_PCP_TC_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/IDAI/precip_composites/"
COMP_PCP_TC_FILE = 'idai_ukmo_nwp.comp_pcp_tc.%d_days.%d%02d.nc'

COMP_PCP_TOT_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/IDAI/precip_composites/"
COMP_PCP_TOT_FILE = 'idai_ukmo_nwp.comp_pcp_tot.%d_days.%d%02d.nc'


#ukmo precip data
MODEL_DIR = '/gws/nopw/j04/klingaman/fascinate/ukmo_nwp/precip/%Y/%Y%m%d%H'
MODEL_FILE = 'prods_op_gl-mn_%Y%m%d_%H_???.nc'
#MODEL_FILE = 'prods_op_gl-mn_%%Y%%m%%d_%%H_%%s.nc'
#MODEL_FILE = 'prods_op_gl-mn_{0}{1}{2}_{3}_{4}.nc'
#before some date in 2007, the data looks different (stored in one file rather than split into many)
SINGLE_MODEL_FILE = 'prods_op_gl-mn_%Y%m%d_%H.nc'

if not os.path.exists(COMP_PCP_TC_DIR):
	os.makedirs(COMP_PCP_TC_DIR)


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
	#if forecast_date >= datetime.datetime(2007, 11, 28, 0):
	model_path = forecast_date.strftime(os.path.join(MODEL_DIR, MODEL_FILE))
	print model_path
	#else:
		#model_path = forecast_date.strftime(os.path.join(MODEL_DIR, SINGLE_MODEL_FILE))
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		#cubes = iris.cube.CubeList([c for c in iris.load(model_path, callback=cb) if c.var_name.startswith('UM_m10s20i013_vn') and not c.var_name.endswith('_1')])
		cubes = iris.cube.CubeList([c for c in iris.load(model_path) if c.var_name.startswith('UM_m01s05i216_vn1009') and not c.var_name.endswith('_1')])
		print cubes[0]
		print len(cubes)
	if not len(cubes):
		raise IOError('No precip data for %s T+%03d' % (forecast_date.strftime('%Y%m%d %H:00'), forecast_hours))

	# Some Cubes have scalar time coord; promote to dim coord
	#for ii, icube in enumerate(cubes):
	if cubes[0].ndim == 2:
		cubes[0] = iris.util.new_axis(cubes[0], cubes[0].coord(axis='T'))
	iris.util.unify_time_units(cubes)

	print cubes
	print cubes[0]
	# Get rid of duplicate times (at the moment, works only if they appeared in
	# a Cube with scalar time dimension)
	for ii, icube in enumerate(cubes):
		if len(icube.coord(axis='T').points) == 1:
			for jcube in (cubes[:ii] + cubes[ii+1:]):
				if icube.coord(axis='T').points[0] in \
						jcube.coord(axis='T').points:
					cubes.remove(icube)
					break
	cube = cubes.concatenate_cube()

	print cube
	#print cube
	#cube = cubes[0]

	# Get correct time
	tcoord = cube.coord(axis='T')
	print tcoord
	#con_func = lambda tt: tt.point == (gt + datetime.timedelta(hours=3))
	#print gt + datetime.timedelta(hours=3)
	#con_func = lambda tt: tt.bound == (gt + datetime.timedelta(hours=3))
	#con = iris.Constraint(**{tcoord.name(): con_func})
	#print con
	#with iris.FUTURE.context(cell_datetime_objects=True):
		#pcp = cube.extract(con)

	#myDate = gt + datetime.timedelta(hours=3)
	#myConstraint = iris.Constraint(time = myDate, cube.coord('time').units, cube.coord('time'.units.calendar))
	#pcp = cube.extract(myConstraint)

	time_units = cube.coord('time').units
	val = time_units.date2num(gt + datetime.timedelta(hours=3))
	constraint = iris.Constraint(time=val)
	pcp = cube.extract(constraint)

	#pcp = cubes[0].extract(con)
	print "PCP"
	print pcp
	if pcp is None:
		raise IOError('###### 2: No precip data for %s T+%03d' % (forecast_date.strftime('%Y%m%d %H:00'), forecast_hours))

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


def composite_pcp_tc_year_month(year, month, lead):
	"""Computes a composite of the precipitation due to all TCs at a particular
	forecast lead time (in days) for a particular year and month

	Total is divided by 2 at the end as each day is composited from both the
	00Z and the 12Z forecast.  (Of course there may be no tracks at certain
	times, in which case that time just contributes 0 to the total.)

	**Arguments**

	*year*, *month*
		`int`s, year and month of validity times for which to calculate the
		composite

	*lead*
		`int`, length of time after forecast initialization
		
	"""
	
		# If no Cubes, create a dummy one with zeros


	# Check lead time is in available range
	if not 0 <= lead <= 6:
		raise ValueError('lead=%s; must be 0 <= lead <= 6' % lead)
	#if year == 2017:
		#if not 1 <= month <= 7:
			#raise ValueError('Data in 2017 used up to July only')

	# Check whether output file already exists

	#infilename = TRACK_FILE %
	infile = os.path.join(TRACK_DIR_3WAY, TRACK_FILE)
	print infile

	outdir = os.path.join(COMP_PCP_TC_DIR, str(lead), str(year))
	comp_file = COMP_PCP_TC_FILE % (lead, year, month)
	outfile = os.path.join(outdir, comp_file)

	print outfile

	if os.path.isfile(outfile):
		raise ValueError('Output file %s already exists' % outfile)
	if not os.path.isdir(outdir):
		os.makedirs(outdir)

	# Iterate for every time available in this year and month
	t1 = datetime.datetime(year, month, 4, 12)
	print "t1:", t1
	t2 = datetime.datetime(year, month, 14, 0)
	print "t2:", t2
	dt = datetime.timedelta(days=1)
	#if (year, month) == (2017, 7):
		#t2 = datetime.datetime(2017, 7, 11)
	#elif month == 12:
		#t2 = datetime.datetime(year+1, 1, 1) - dt
	#else:
		#t2 = datetime.datetime(year, month+1, 1) - dt
	pcp_cubes = iris.cube.CubeList()
	exclude = pl.exclude_days(lead)
	count_days = 0
	#if (year, month) == (2010, 3):
		#count_days_res = {'n320': 0.0, 'n512': 0.0}
		#start_n512 = datetime.datetime(2010, 3, 9, 12)
	#elif (year, month) == (2014, 7):
		#count_days_res = {'n512': 0.0, 'n768': 0.0}
		#start_n768 = datetime.datetime(2014, 7, 15, 12)
	vt_list = []
	print "EXCLUDED DATES:", exclude
	while t1 <= t2:

		# Check whether this day is on the list of those to exclude
		if t1 in exclude:
			print t1, '- EXCLUDE'
			t1 += dt
			continue
		#if t1.timetuple()[:3] == (2017, 7, 11):
			#count_days += 0.5
		#else:
		count_days += 1
		print t1.strftime('%Y/%m/%d')

		# Get list of forecast and validity times for the three forecasts to be
		# used
		ftime_deltas = np.arange(-12, 13, 12) - lead*24
		print "ftime_deltas:", ftime_deltas
		ftimes=[]
		for hh in ftime_deltas:
				ftimes.append(t1 + datetime.timedelta(hours=hh))
		#ftimes = (t1 + datetime.timedelta(hours=hh) for hh in ftime_deltas)
		print "ftimes:", ftimes
		vtimes = (np.array([15, 21]) + lead*24, np.arange(3, 22, 6) + lead*24, np.array([3, 9]) + lead*24)
		print "vtimes:", vtimes




		# Iterate for each of the three forecasts
		for ff, vv in itertools.izip(ftimes, vtimes):
			print "vv:", vv

			# If on or after 2017/07/11 12:00, skip
			#if ff >= datetime.datetime(2017, 7, 11, 12):
				#continue

			# Get year, month, day, hour, lon, lat from file
			this_infile = ff.strftime(infile)
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				ain = np.genfromtxt(this_infile, dtype=float, skip_header=1, usecols=range(3, 9))

			# Count days for each resolution
			#for v in vv:
				#vt = ff + datetime.timedelta(hours=v)
				#if vt not in vt_list:
					#vt_list.append(vt)
					#if (year, month) == (2010, 3):
						#if vt < start_n512:
							#count_days_res['n320'] += 0.25
						#else:
							#count_days_res['n512'] += 0.25
					#elif (year, month) == (2014, 7):
						#if vt < start_n768:
							#count_days_res['n512'] += 0.25
						#else:
							#count_days_res['n768'] += 0.25

			# If no tracks in this forecast, skip it
			if not ain.size:
				print '   ', ff, '- no tracks'
				continue

			# Iterate for every validity time required from this forecast
			for v in vv:
				print "v:", v

				# Get track(s) with point(s) this time
				gd = ff + datetime.timedelta(hours=v)
				print "gd:", gd
				aint = ain[np.where((ain[:, 0] == gd.year) &\
									(ain[:, 1] == gd.month) &\
									(ain[:, 2] == gd.day) &\
									(ain[:, 3] == gd.hour))]

				if not aint.size:
					print '   ', ff, 'T+%03d' % v, '- no tracks'
					continue
				print '   ', ff, 'T+%03d' % v



				# Iterate for each track
				print "running nwp_pcp_accumulation function for points along track:"
				#print "ff:", ff
				#print "v:", v
				#print "aint:", aint
				for lon, lat in aint[:, [4, 5]]:
					print aint[:,[4,5]]
					print "this_pcp = nwp_pcp_accumulation(",ff, v, lon, lat,")"
					this_pcp = nwp_pcp_accumulation(ff, v, lon, lat)
					print "finished running nwp_pcp_accumulation"
					if this_pcp is None:
						print "no precip, running dummy_cube()"
						this_pcp = iris.cube.CubeList([dummy_cube()])
					else:
						this_pcp.coord(axis='X').var_name = 'longitude'
						this_pcp.coord(axis='Y').var_name = 'latitude'
						this_pcp.coord(axis='X').attributes = {}
						this_pcp.coord(axis='Y').attributes = {}
						pcp_cubes.append(iris.util.squeeze(this_pcp))

		# Increment time
		t1 += dt

	def dummy_cube():
		dummy = None
		dummy_t = datetime.datetime(year, month, 1)
		while dummy is None:
			dummy = this_pcp = cf.nwp_pcp_accumulation(dummy_t, 3)
			dummy_t += dt
		dummy = iris.util.squeeze(dummy)
		dummy.data = np.zeros_like(dummy.data)
		dummy.remove_coord(dummy.coord(axis='T'))
		return dummy
	if not len(pcp_cubes):
		pcp_cubes = iris.cube.CubeList([dummy_cube()])

	if not len(pcp_cubes):
		pcp_cubes = iris.cube.CubeList([dummy_cube()])

	# Sum over Cubes and divide by 2
	pcp = pl.add_cubes(pcp_cubes, deal_with_masks=False, contributing_days=False)/2.

	# Set metadata
	pcp.units = 'mm'
	pcp.standard_name = 'lwe_thickness_of_precipitation_amount'
	pcp.long_name = 'precipitation'
	pcp.var_name = 'pcp'
	pcp.attributes['contributing_days'] = count_days

	# Save
	iris.save(pcp, outfile)
	print "saved: ", outfile

	# For months with more than one resolution, sum separately and divide by 2
	#if (year, month) in [(2010, 3), (2014, 7)]:
		#if year == 2010:
			#res_list= ['n320','n512']
		#elif year == 2014:
			#res_list = ['n512', 'n768']
		#res_list = {2010: ['n320', 'n512'], 2014: ['n512', 'n768']}[year]
		#pcp_sep = pl.add_cubes(pcp_cubes, deal_with_masks=False, separate_resolutions=True, contributing_days=False)
		#file_tot = os.path.join(COMP_PCP_TOT_DIR, str(lead), str(year),COMP_PCP_TOT_FILE % (lead, year, month))
		#for k in pcp_sep.iterkeys():
			#pcp_sep_k = pcp_sep[k]/2.

			# Set metadata
			#pcp_sep_k.units = 'mm'
			#pcp_sep_k.standard_name = 'lwe_thickness_of_precipitation_amount'
			#pcp_sep_k.long_name = 'precipitation'
			#pcp_sep_k.var_name = 'pcp'

			# Number of contributing days is difficult to count so just get the
			# value from the total pcp composites (the number should be the
			# same anyway)
			#res = {640: 'n320', 1024: 'n512', 1536: 'n768'}[k[1]]
			#res_list.remove(res)
			#file_tot_k = file_tot.replace('.nc', '%s.nc' % res)
			#cube_tot_k = iris.load_cube(file_tot_k)
			#pcp_sep_k.attributes['contributing_days'] = \
				#float(cube_tot_k.attributes['contributing_days'])

			# Save
			#outfile_k = outfile.replace('.nc', '.%s.nc' % res)
			#iris.save(pcp_sep_k, outfile_k)
			#print outfile_k

		# If any resolutions are still in res_list it means there were no
		# tracks at that resolution, so save an empty Cube
		#for res in res_list:
			#pcp_sep_k = dummy_cube()
			#pcp_sep_k.units = 'mm'
			#pcp_sep_k.standard_name = 'lwe_thickness_of_precipitation_amount'
			#pcp_sep_k.long_name = 'precipitation'
			#pcp_sep_k.var_name = 'pcp'
			#file_tot_k = file_tot.replace('.nc', '%s.nc' % res)
			#cube_tot_k = iris.load_cube(file_tot_k)
			#pcp_sep_k.attributes['contributing_days'] = float(cube_tot_k.attributes['contributing_days'])
			#outfile_k = outfile.replace('.nc', '.%s.nc' % res)
			#iris.save(pcp_sep_k, outfile_k)
			#print outfile_k

def composite_pcp_tc_year(year1, year2, lead):
	"""Computes the total composite precipitation due to observed TCs for a particular
	TC season (July year1 - June year2)."""

	infiles = []
	for m in [7, 8, 9, 10, 11, 12]:  # 7,
		# m = str(m).zfill(2)
		filename = COMP_PCP_TC_FILE % (lead, year1, m)
		infile = os.path.join(COMP_PCP_TC_DIR, str(lead), str(year1), filename )
		infiles.append(infile)

	for m in [1, 2, 3, 4, 5, 6]:
		# m=str(m).zfill(2)
		# year=year2
		filename = COMP_PCP_TC_FILE % (lead, year2, m)
		infile = os.path.join(COMP_PCP_TC_DIR, str(lead), str(year2), filename)
		infiles.append(infile)


	print infiles

	outfile = "ukmo_nwp.comp_pcp_tc.%d-days.%d-%d.nc" % (lead, year1, year2)
	outpath = os.path.join(COMP_PCP_TC_DIR, str(lead), outfile)

	pl.add_files(infiles, outpath, cube_count=12)

def composite_pcp_tc_year_separate_resolutions(year1, year2, lead, threeway=False):
	"""Computes the total composite precipitation due to TCs at a particular
	forecast lead time (in days) for a particular year."""

	# Check whether separation of resolutions is required
	#if year not in [2010, 2014]:
	#    raise ValueError('All data in %d are at the same resolution' % year)

	# Get list of infiles
	infiles = []
	for m in [7, 8, 9, 10, 11, 12]:  # 7,
		if (year1, m) == (2014, 7):
			res = ['n512', 'n768']
			infile = [os.path.join(COMP_PCP_TC_DIR, str(lead), str(year1), COMP_PCP_TC_FILE % (lead, year1, m))] * 2
			infile = [i.replace('.nc', '%s.nc' % r) for i, r in zip(infile, res)]
			for a in infile:
				infiles.append(a)
		else:
			# m = str(m).zfill(2)
			filename = COMP_PCP_TC_FILE % (lead, year1, m)
			infile = os.path.join(COMP_PCP_TC_DIR, str(lead), str(year1), filename)
			infiles.append(infile)

	for m in [1, 2, 3, 4, 5, 6]:
		if (year2, m) == (2010, 3):
			res = ['n320', 'n512']
			infile = [os.path.join(COMP_PCP_TC_DIR, str(lead), str(year2), COMP_PCP_TC_FILE % (lead, year2, m))] * 2
			infile = [i.replace('.nc', '%s.nc' % r) for i, r in zip(infile, res)]
			for a in infile:
				infiles.append(a)
		# m=str(m).zfill(2)
		# year=year2
		else:
			filename = COMP_PCP_TC_FILE % (lead, year2, m)
			infile = os.path.join(COMP_PCP_TC_DIR, str(lead), str(year2), filename)
			infiles.append(infile)

	print infiles

	outfile = "ukmo_nwp.comp_pcp_tc.%d-days.%d-%d.nc" % (lead, year1, year2)
	outpath = os.path.join(COMP_PCP_TC_DIR, str(lead), outfile)

	pl.add_files(infiles, outpath, separate_resolutions=True, deal_with_masks=False, cube_count=13)




composite_pcp_tc_year_month(2019,3,lead)

#composite_pcp_tc_year_month(2006, 12, 0)

#Need contributing days from total (not tc-related) UKMO NWP precip composites, for resolution changes
#So need to produce the total precip composites before I can run 2010 and 2014, and 2009 is still missing

#if lead == 6:
	#no lead time 7 days ahead for 2006-2008
	#year1s = [2010, 2011, 2012, 2013, 2015] #2006, 2007,  #will need to add 2008, 2009
	#year2s = [2011, 2012, 2013, 2014, 2016] #2007, 2008,  #will need to add 2009, 2010
#else:
	#year1s = [2006, 2007, 2008, 2010, 2011, 2012, 2013, 2015] #2006, 2007,  #will need to add 2008, 2009
	#year2s = [2007, 2008, 2009, 2011, 2012, 2013, 2014, 2016] #2007, 2008,  #will need to add 2009, 2010

#year1s_sep_res = [2014] #2009,
#year2s_sep_res = [2015] #2010,
#year1s = [2006]
#year2s = [2007]


#run monthly precip totals:
#for year in [2014]: #
	#no lead time 6 in 2006, and also only run from July onwards
	#if year == 2006:
		#if lead in [0,6]:
			#continue
		#else:
			##for month in [7,8,9,10,11,12]:
				#composite_pcp_tc_year_month(year, month, lead)

	#no lead time 6 in 2007 and 2008
	#elif year in [2007, 2008]:
		#if lead == 6:
			#continue
		#else:
			#for month in [1,2,3,4,5,6,7,8,9,10,11,12]:
				#composite_pcp_tc_year_month(year, month, lead)

	#if year == 2014:
		#for month in [7]:
			#composite_pcp_tc_year_month(year, month, lead)

	#else:
		#for month in [1,2,3,4,5,6,7,8,9,10,11,12]:
			#composite_pcp_tc_year_month(year, month, lead)

#run yearly (tc season) precip totals:
#for year1, year2 in zip(year1s, year2s):
	#composite_pcp_tc_year(year1, year2, lead)

#for year1, year2 in zip(year1s_sep_res, year2s_sep_res):
	#composite_pcp_tc_year_separate_resolutions(year1,year2,lead)
