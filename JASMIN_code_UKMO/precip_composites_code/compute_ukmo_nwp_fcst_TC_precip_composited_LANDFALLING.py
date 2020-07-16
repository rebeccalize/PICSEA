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
TRACK_DIR_3WAY = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/seychelles_landfalling/"
TRACK_FILE='ukmo_nwp.%Y%m%d%H.tcident.interp' #ukmo_nwp.2006111400.tcident.interp

#outdata
#COMP_PCP_TC_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_precip_composites/"
COMP_PCP_TC_DIR = './'
COMP_PCP_TC_FILE = 'ukmo_nwp.comp_pcp_tc.%d_days.%d%02d.nc'

COMP_PCP_TOT_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_TOTAL_precip_composites/"
COMP_PCP_TOT_FILE = 'ukmo_nwp.comp_pcp_tot.%d_days.%d%02d.nc'


#ukmo precip data
#Simon used "test" directory instead of "precip" - data looks the same, what's going on?
MODEL_DIR = '/gws/nopw/j04/klingaman/fascinate/ukmo_nwp/precip/%Y/%Y%m%d%H'
MODEL_FILE = 'prods_op_gl-mn_%Y%m%d_%H_???.nc'
#before some date in 2007, the data looks different (stored in one file rather than split into many)
SINGLE_MODEL_FILE = 'prods_op_gl-mn_%Y%m%d_%H.nc'

if not os.path.exists(COMP_PCP_TC_DIR):
	os.makedirs(COMP_PCP_TC_DIR)


def composite_pcp_tc_year_month(year, month, lead):
	"""Computes a composite of the precipitation due to all TCs at a particular
	forecast lead time (in days) for a particular year and month.

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

	# Check lead time is in available range
	if not 0 <= lead <= 6:
		raise ValueError('lead=%s; must be 0 <= lead <= 6' % lead)
	#if year == 2017:
		#if not 1 <= month <= 7:
			#raise ValueError('Data in 2017 used up to July only')

	# Check whether output file already exists

	#infilename = TRACK_FILE %
	infile = os.path.join(TRACK_DIR_3WAY, '%Y', TRACK_FILE)
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
	t1 = datetime.datetime(year, month, 1, 0)
	dt = datetime.timedelta(days=1)
	if (year, month) == (2017, 7):
		t2 = datetime.datetime(2017, 7, 11)
	elif month == 12:
		t2 = datetime.datetime(year+1, 1, 1) - dt
	else:
		t2 = datetime.datetime(year, month+1, 1) - dt
	pcp_cubes = iris.cube.CubeList()
	exclude = pl.exclude_days(lead)
	count_days = 0
	if (year, month) == (2010, 3):
		count_days_res = {'n320': 0.0, 'n512': 0.0}
		start_n512 = datetime.datetime(2010, 3, 9, 12)
	elif (year, month) == (2014, 7):
		count_days_res = {'n512': 0.0, 'n768': 0.0}
		start_n768 = datetime.datetime(2014, 7, 15, 12)
	vt_list = []
	while t1 <= t2:

		# Check whether this day is on the list of those to exclude
		if t1 in exclude:
			print t1, '- EXCLUDE'
			t1 += dt
			continue
		if t1.timetuple()[:3] == (2017, 7, 11):
			count_days += 0.5
		else:
			count_days += 1
		print t1.strftime('%Y/%m/%d')

		# Get list of forecast and validity times for the three forecasts to be
		# used
		ftime_deltas = np.arange(-12, 13, 12) - lead*24
		ftimes = (t1 + datetime.timedelta(hours=hh) for hh in ftime_deltas)
		vtimes = (np.array([15, 21]) + lead*24, np.arange(3, 22, 6) + lead*24,
				  np.array([3, 9]) + lead*24)

		# Iterate for each of the three forecasts
		for ff, vv in itertools.izip(ftimes, vtimes):

			# If on or after 2017/07/11 12:00, skip
			#if ff >= datetime.datetime(2017, 7, 11, 12):
				#continue

			# Get year, month, day, hour, lon, lat from file
			this_infile = ff.strftime(infile)
			with warnings.catch_warnings():
				warnings.simplefilter('ignore')
				ain = np.genfromtxt(this_infile, dtype=float, skip_header=1, usecols=range(3, 9))

			# Count days for each resolution
			for v in vv:
				vt = ff + datetime.timedelta(hours=v)
				if vt not in vt_list:
					vt_list.append(vt)
					if (year, month) == (2010, 3):
						if vt < start_n512:
							count_days_res['n320'] += 0.25
						else:
							count_days_res['n512'] += 0.25
					elif (year, month) == (2014, 7):
						if vt < start_n768:
							count_days_res['n512'] += 0.25
						else:
							count_days_res['n768'] += 0.25

			# If no tracks in this forecast, skip it
			if not ain.size:
				print '   ', ff, '- no tracks'
				continue

			# Iterate for every validity time required from this forecast
			for v in vv:

				# Get track(s) with point(s) this time
				gd = ff + datetime.timedelta(hours=v)
				aint = ain[np.where((ain[:, 0] == gd.year) &\
									(ain[:, 1] == gd.month) &\
									(ain[:, 2] == gd.day) &\
									(ain[:, 3] == gd.hour))]
				if not aint.size:
					print '   ', ff, 'T+%03d' % v, '- no tracks'
					continue
				print '   ', ff, 'T+%03d' % v

				# Iterate for each track
				for lon, lat in aint[:, [4, 5]]:
					this_pcp = cf.nwp_pcp_accumulation(ff, v, lon, lat)
					this_pcp.coord(axis='X').var_name = 'longitude'
					this_pcp.coord(axis='Y').var_name = 'latitude'
					this_pcp.coord(axis='X').attributes = {}
					this_pcp.coord(axis='Y').attributes = {}
					pcp_cubes.append(iris.util.squeeze(this_pcp))

		# Increment time
		t1 += dt

	# If no Cubes, create a dummy one with zeros
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
	print outfile

	# For months with more than one resolution, sum separately and divide by 2
	if (year, month) in [(2010, 3), (2014, 7)]:
		if year == 2010:
			res_list= ['n320','n512']
		elif year == 2014:
			res_list = ['n512', 'n768']
		#res_list = {2010: ['n320', 'n512'], 2014: ['n512', 'n768']}[year]
		pcp_sep = pl.add_cubes(pcp_cubes, deal_with_masks=False, separate_resolutions=True, contributing_days=False)
		file_tot = os.path.join(COMP_PCP_TOT_DIR, str(lead), str(year),COMP_PCP_TOT_FILE % (lead, year, month))
		for k in pcp_sep.iterkeys():
			pcp_sep_k = pcp_sep[k]/2.

			# Set metadata
			pcp_sep_k.units = 'mm'
			pcp_sep_k.standard_name = 'lwe_thickness_of_precipitation_amount'
			pcp_sep_k.long_name = 'precipitation'
			pcp_sep_k.var_name = 'pcp'

			# Number of contributing days is difficult to count so just get the
			# value from the total pcp composites (the number should be the
			# same anyway)
			res = {640: 'n320', 1024: 'n512', 1536: 'n768'}[k[1]]
			res_list.remove(res)
			file_tot_k = file_tot.replace('.nc', '%s.nc' % res)
			cube_tot_k = iris.load_cube(file_tot_k)
			pcp_sep_k.attributes['contributing_days'] = \
				float(cube_tot_k.attributes['contributing_days'])

			# Save
			outfile_k = outfile.replace('.nc', '.%s.nc' % res)
			iris.save(pcp_sep_k, outfile_k)
			print outfile_k

		# If any resolutions are still in res_list it means there were no
		# tracks at that resolution, so save an empty Cube
		for res in res_list:
			pcp_sep_k = dummy_cube()
			pcp_sep_k.units = 'mm'
			pcp_sep_k.standard_name = 'lwe_thickness_of_precipitation_amount'
			pcp_sep_k.long_name = 'precipitation'
			pcp_sep_k.var_name = 'pcp'
			file_tot_k = file_tot.replace('.nc', '%s.nc' % res)
			cube_tot_k = iris.load_cube(file_tot_k)
			pcp_sep_k.attributes['contributing_days'] = float(cube_tot_k.attributes['contributing_days'])
			outfile_k = outfile.replace('.nc', '.%s.nc' % res)
			iris.save(pcp_sep_k, outfile_k)
			print outfile_k

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


#composite_pcp_tc_year_month(2006, 12, 0)

#Need contributing days from total (not tc-related) UKMO NWP precip composites, for resolution changes
#So need to produce the total precip composites before I can run 2010 and 2014, and 2009 is still missing

if lead == 6:
	#no lead time 7 days ahead for 2006-2008
	year1s = [2010, 2011, 2012, 2013, 2015] #2006, 2007,  #will need to add 2008, 2009
	year2s = [2011, 2012, 2013, 2014, 2016] #2007, 2008,  #will need to add 2009, 2010
else:
	year1s = [2006, 2007, 2008, 2010, 2011, 2012, 2013, 2015] #2006, 2007,  #will need to add 2008, 2009
	year2s = [2007, 2008, 2009, 2011, 2012, 2013, 2014, 2016] #2007, 2008,  #will need to add 2009, 2010

year1s_sep_res = [2014] #2009,
year2s_sep_res = [2015] #2010,
#year1s = [2006]
#year2s = [2007]


#run monthly precip totals:
for year in [2014]: #
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

	if year == 2014:
		for month in [1]:
			composite_pcp_tc_year_month(year, month, lead)

	#else:
		#for month in [1,2,3,4,5,6,7,8,9,10,11,12]:
			#composite_pcp_tc_year_month(year, month, lead)

#run yearly (tc season) precip totals:
#for year1, year2 in zip(year1s, year2s):
	#composite_pcp_tc_year(year1, year2, lead)

#for year1, year2 in zip(year1s_sep_res, year2s_sep_res):
	#composite_pcp_tc_year_separate_resolutions(year1,year2,lead)
