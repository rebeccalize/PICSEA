import picsea_library as pl
import composite_functions_ukmo_nwp as cf
import datetime
import iris
import itertools
import numpy as np
import os

#IBTRACS_DIR = "/gws/nopw/j04/klingaman/emerton/ibtracs_test_data_maps/"
#IBTRACS_FILE = "storms_SH%d%d_interpolated.txt"


COMP_PCP_TOT_FILE = "trmm.comp_pcp_total.%d%02d.nc"

COMP_PCP_TOT_DIR = "/gws/nopw/j04/klingaman/emerton/total_trmm_precip_composites/"
if not os.path.exists(COMP_PCP_TOT_DIR):
	os.makedirs(COMP_PCP_TOT_DIR)


# pcp = pl.read_trmm_3b42(2012, 11, 10, 3)
# print pcp


def composite_pcp_tot_year_month(year, month, lead=None, mjo=None, y1=2006, y2=2017):

	"""Computes the total precipitation from TRMM 3B42 for a particular year and month. For 2017, uses data up to 2017/07/11 00:00 (inclusive) only."""
	print "starting"
	# Check values are in valid ranges
	if lead is not None:
		if not 0 <= lead <= 6:
			raise ValueError('lead=%s; must be 0 <= lead <= 6' % lead)
	if mjo is not None:
		if not 0 <= mjo <= 8:
			raise ValueError('mjo=%s; must be 0 <= mjo <= 8' % mjo)
	# if year == 2017:
	#     if not 1 <= month <= 7:
	#         raise ValueError('Data in 2017 used up to July only')

	# Compute file names and get list of dates to exclude due to missing
	# ukmo_nwp forecasts
	if lead is None:
		if mjo is None:
			outfile = COMP_PCP_TOT_FILE % (year, month)
		else:
			outfile = COMP_PCP_TOT_FILE_MJO % (year, month, mjo)
		exclude = []
	else:
		if not 0 <= lead <= 6:
			raise ValueError('lead=%s; must be None or 0 <= lead <= 6' % lead)
		if mjo is None:
			outfile = COMP_PCP_TOT_FILE_LEAD % (lead, year, month)
		else:
			outfile = COMP_PCP_TOT_FILE_LEAD_MJO % (lead, year, month, mjo)
		exclude = pl.exclude_days(lead)
	outpath = os.path.join(COMP_PCP_TOT_DIR,str(year), outfile)
	if not os.path.exists(os.path.join(COMP_PCP_TOT_DIR,str(year))):
		os.makedirs(os.path.join(COMP_PCP_TOT_DIR,str(year)))
	if os.path.isfile(outpath):
		raise ValueError('Output file %s already exists' % outpath)

	# For a given MJO phase, get all cyclone tracks for this year and month,
	# and exclude those whose genesis was in the wrong month; then get list of
	# times for remaining track points
	if mjo is not None:
		tc_file = os.path.join(TRACK_DIR, TRACK_FILE % (1979, 2016))
		ain = np.genfromtxt(tc_file, dtype=float, skip_header=1,
							usecols=[0]+range(3, 9))
		ain_y = ain[np.where(ain[:, 1] == year)]
		ain_ym = ain_y[np.where(ain_y[:, 2] == month)]
		fn = os.path.join(TRACK_DIR, TRACK_FILE_0 % (1979, 2016) + \
							  '.genesis.mjo%d_trackIDs'% mjo)
		ids_list = np.genfromtxt(fn)
		ain_use = ain_ym[np.in1d(ain_ym[:, 0], ids_list)]
		ain_dates = ain_use[:, range(1, 5)]
		dates = [datetime.datetime(*d) for d in ain_dates.astype(int)]


	# Iterate for every day available in this year and month
	t1 = datetime.datetime(year, month, 1)
	dt = datetime.timedelta(days=1)
	if (year, month) == (2017, 7):
		t2 = datetime.datetime(2017, 7, 11)
	elif month == 12:
		t2 = datetime.datetime(year, 12, 31)
	else:
		t2 = datetime.datetime(year, month+1, 1) - dt
	pcp_cubes = iris.cube.CubeList()
	count_days = 0
	if (year, month) == (2010, 3):
		pcp_cubes_res = {'n320': iris.cube.CubeList(),
						 'n512': iris.cube.CubeList()}
		count_days_res = {'n320': 0.0, 'n512': 0.0}
		start_n512 = datetime.datetime(2010, 3, 9, 12)
	elif (year, month) == (2014, 7):
		pcp_cubes_res = {'n512': iris.cube.CubeList(),
						 'n768': iris.cube.CubeList()}
		count_days_res = {'n512': 0.0, 'n768': 0.0}
		start_n768 = datetime.datetime(2014, 7, 15, 12)
	while t1 <= t2:

		# If on list of dates to exclude, move on
		y, m, d = [getattr(t1, a) for a in ['year', 'month', 'day']]
		if t1 in exclude:
			print t1.strftime('%Y/%m/%d'), '-- EXCLUDE'
			t1 += dt
			continue

		# Iterate for each 6 hour period in day
		# (Remember, h=3 means 00-06UTC, etc.)
		if t1.timetuple()[:3] == (2017, 7, 11):
			count_days += 0.5
		else:
			count_days += 1
		for h in xrange(3, 22, 6):
			this_t_h = t1 + datetime.timedelta(hours=h)
			if this_t_h >= datetime.datetime(2017, 7, 11, 12):
				t1 = t2 + dt
				break
			if (year, month) == (2010, 3):
				if (t1 + datetime.timedelta(hours=h)) < start_n512:
					count_days_res['n320'] += 0.25
				else:
					count_days_res['n512'] += 0.25
			elif (year, month) == (2014, 7):
				if (t1 + datetime.timedelta(hours=h)) < start_n768:
					count_days_res['n512'] += 0.25
				else:
					count_days_res['n768'] += 0.25

			# For a given MJO phase, skip if this time is not on the list found
			# above (note that we are doing this AFTER counting the
			# contributing days above)
			if mjo is not None:
				if (t1 + datetime.timedelta(hours=h)) not in dates:
					continue
			print t1.strftime('%Y/%m/%d'), '%02d:00' % h

			# Get data
			this_pcp = cf.trmm_3b42_pcp_accumulation(y, m, d, h)
			if this_pcp is not None:
				if this_pcp.data.min() < 0:
					raise ValueError('Found negative precip value(s)')
				pcp_cubes.append(this_pcp)

				# For months including a change of resolution, append Cube to
				# the appropriate CubeList
				if (year, month) == (2010, 3):
					if datetime.datetime(y, m, d, h) < start_n512:
						pcp_cubes_res['n320'].append(this_pcp)
					else:
						pcp_cubes_res['n512'].append(this_pcp)
				elif (year, month) == (2014, 7):
					if datetime.datetime(y, m, d, h) < start_n768:
						pcp_cubes_res['n512'].append(this_pcp)
					else:
						pcp_cubes_res['n768'].append(this_pcp)

		# Increment day
		t1 += dt

	# Sum over Cubes
	if len(pcp_cubes):
		pcp_sum = pl.add_cubes(pcp_cubes, contributing_days=False)
	else:
		print 'No data - creating dummy map'
		dummy = cf.trmm_3b42_pcp_accumulation(1998, 1, 1, 3)
		pcp_sum = pl.add_cubes([dummy, dummy], contributing_days=False)
		pcp_sum.data = np.zeros_like(pcp_sum.data)

	# Set metadata
	pcp_sum.units = 'mm'
	pcp_sum.standard_name = 'lwe_thickness_of_precipitation_amount'
	pcp_sum.long_name = 'precipitation'
	pcp_sum.var_name = 'pcp'
	pcp_sum.attributes['contributing_days'] = float(count_days)

	# Save
	print outpath
	iris.save(pcp_sum, outpath)


	# For months with more than one resolution, sum separately and divide by 2
	if (year, month) in [(2010, 3), (2014, 7)]:
		for k in pcp_cubes_res.iterkeys():
			if len(pcp_cubes_res[k]):
				pcp_k = pl.add_cubes(pcp_cubes_res[k],
										   deal_with_masks=False,
										   contributing_days=False)/2.
			else:
				print 'No data - creating dummy map'
				dummy = cf.trmm_3b42_pcp_accumulation(1998, 1, 1, 3)
				pcp_k = pl.add_cubes([dummy, dummy],
											 contributing_days=False)
				pcp_k.data = np.zeros_like(pcp_k.data)

			# Set metadata
			pcp_k.units = 'mm'
			pcp_k.standard_name = 'lwe_thickness_of_precipitation_amount'
			pcp_k.long_name = 'precipitation'
			pcp_k.var_name = 'pcp'
			pcp_k.attributes['contributing_days'] = count_days_res[k]

			# Save
			outpath_k = outpath.replace('.nc', '.%s.nc' % k)
			iris.save(pcp_k, outpath_k)
			#print outpath_k
		


def composite_pcp_tot_year(year1, year2, lead=None):
	"""Computes the total composite precipitation due to observed TCs for a particular
	TC season (July year1 - June year2)."""

	infiles = []
	for m in [7, 8, 9, 10, 11, 12]:  # 7,
		# m = str(m).zfill(2)
		filename = COMP_PCP_TOT_FILE % (year1, m)
		infile = os.path.join(COMP_PCP_TOT_DIR, str(year1), filename )
		infiles.append(infile)

	for m in [1, 2, 3, 4, 5, 6]:
		# m=str(m).zfill(2)
		# year=year2
		filename = COMP_PCP_TOT_FILE % (year2, m)
		infile = os.path.join(COMP_PCP_TOT_DIR, str(year2), filename)
		infiles.append(infile)


	print infiles

	outfile = "trmm.comp_pcp_total.%d-%d.nc" % (year1, year2)
	outpath = os.path.join(COMP_PCP_TOT_DIR, outfile)

	pl.add_files(infiles, outpath, cube_count=12)

for year in [2016, 2017, 2018]: #2006, 2007, 2008, 2009, 2010, 2011,2012, 2013, 2014, 2015, 2016
	if year == 2006:
		for month in range(7,13):
			composite_pcp_tot_year_month(year, month)
	elif year == 2016:
		for month in range(7,13):
			composite_pcp_tot_year_month(year, month)
	elif year == 2018:
		for month in range(1,7):
			composite_pcp_tot_year_month(year, month)
	else:
		for month in range(1,13):
			composite_pcp_tot_year_month(year, month)


#year1s = [2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015]
#year2s = [ 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016]
year1s = [2016, 2017]
year2s = [2017, 2018]

for year1, year2 in zip(year1s, year2s):
	composite_pcp_tot_year(year1, year2, lead=None)