import picsea_library as pl
import composite_functions_ukmo_nwp as cf
import datetime
import iris
import itertools
import numpy as np
import os
import sys

lead=int(sys.argv[1])
lead=lead-1

#indata
TRACK_DIR_3WAY = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_year/"
TRACK_FILE='ukmo_nwp.%Y%m%d%H.tcident.interp' #ukmo_nwp.2006111400.tcident.interp

#outdata
COMP_PCP_TOT_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_TOTAL_precip_composites/"
if not os.path.exists(COMP_PCP_TOT_DIR):
	os.makedirs(COMP_PCP_TOT_DIR)

COMP_PCP_TOT_FILE = 'ukmo_nwp.comp_pcp_tot.%d_days.%d%02d.nc'

def composite_pcp_tot_year_month(year, month, lead):
	"""Computes a composite of all precipitation at a particular forecast lead
	time (in days) for a particular year and month.

	Total is divided by 2 at the end as each day is composited from both the
	00Z and the 12Z forecast.

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

	# Check whether output directory and file already exist
	outfile = os.path.join(COMP_PCP_TOT_DIR, str(lead), str(year),
						   COMP_PCP_TOT_FILE % (lead, year, month))
	if os.path.isfile(outfile):
		raise ValueError('Output file %s already exists' % outfile)

	# Iterate for every time available in this year and month
	########################################################################################
	# MODIFIED BELOW TO TEST SCRIPT, SHOULD BE t1 = datetime.datetime(year, month, 1, 0)
	########################################################################################
	
	#t1 = datetime.datetime(year, month, 1, 0)
	
	t1 = datetime.datetime(year, month, 25, 0) #testing
	
	dt = datetime.timedelta(days=1)
	#if (year, month) == (2017, 7):
		#t2 = datetime.datetime(2017, 7, 11, 12)
	if month == 12:
		t2 = datetime.datetime(year, 12, 31, 12)
	else:
		t2 = datetime.datetime(year, month+1, 1, 12) - dt
	pcp_cubes = iris.cube.CubeList()
	exclude = pl.exclude_days(lead=lead)
	count_days = 0
	if (year, month) == (2010, 3):
		count_days_res = {'n320': 0.0, 'n512': 0.0}
		start_n512 = datetime.datetime(2010, 3, 9, 12)
	elif (year, month) == (2014, 7):
		count_days_res = {'n512': 0.0, 'n768': 0.0}
		start_n768 = datetime.datetime(2014, 7, 15, 12)
	elif (year, month) == (2017,7):
		count_days_res = {'n768': 0.0, 'n1280': 0.0}
		start_n1280 = datetime.datetime(2017, 7, 11, 12)
		
	vt_list = []
	while t1 <= t2:

		# Check whether this day is on the list of those to exclude
		if t1 in exclude:
			print t1.strftime('%Y/%m/%d'), '- EXCLUDE'
			t1 += dt
			continue
		#if t1.timetuple()[:3] == (2017, 7, 11):
			#count_days += 0.5
		else:
			count_days += 1
		print t1.strftime('%Y/%m/%d')

		# Get list of forecast and validity times for the three forecasts to be
		# used
		ftime_deltas = np.arange(-12, 13, 12) - lead*24
		ftimes = [t1 + datetime.timedelta(hours=hh) for hh in ftime_deltas]
		vtimes = (np.array([15, 21]) + lead*24, np.arange(3, 22, 6) + lead*24,
				  np.array([3, 9]) + lead*24)

		# Iterate for each of the three forecasts
		for ff, vv in itertools.izip(ftimes, vtimes):
			#if ff >= datetime.datetime(2017, 7, 11, 12):
				#continue

			# Iterate for every validity time required from this forecast
			for v in vv:
				print '   ', ff, 'T+%03d' % v

				# Count days for each resolution
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
					elif (year, month) == (2017,7):
						if vt < start_n1280:
							count_days_res['n768'] += 0.25
						else:
							count_days_res['n1280'] += 0.25
							

				# Get precipitation data
				this_pcp = cf.nwp_pcp_accumulation(ff, v)
				this_pcp.coord(axis='X').var_name = 'longitude'
				this_pcp.coord(axis='Y').var_name = 'latitude'
				this_pcp.coord(axis='X').attributes = {}
				this_pcp.coord(axis='Y').attributes = {}
				pcp_cubes.append(iris.util.squeeze(this_pcp))

		# Increment time
		t1 += dt

	# Sum over Cubes and divide by 2
	pcp = pl.add_cubes(pcp_cubes, deal_with_masks=False, contributing_days=False)/2.

	# Set metadata
	pcp.units = 'mm'
	pcp.standard_name = 'lwe_thickness_of_precipitation_amount'
	pcp.long_name = 'precipitation'
	pcp.var_name = 'pcp'
	pcp.attributes['contributing_days'] = float(count_days)
	# Save
	print outfile
	iris.save(pcp, outfile)
	#print outfile
	# For months with more than one resolution, sum separately and divide by 2
	if (year, month) in [(2010, 3), (2014, 7), (2017,7)]:
		print "pcp_cubes:", pcp_cubes
		pcp_sep = pl.add_cubes_sep_res(pcp_cubes, year, deal_with_masks=False, separate_resolutions=True, contributing_days=False)

		for k in pcp_sep.iterkeys():
			print k
			print k[1]

		for k in pcp_sep.iterkeys():
			pcp_sep_k = pcp_sep[k]/2.
			#res = {640: 'n320', 1024: 'n512', 1536: 'n768'}[k[1]]
			res=str(k)

			# Set metadata
			pcp_sep_k.units = 'mm'
			pcp_sep_k.standard_name = 'lwe_thickness_of_precipitation_amount'
			pcp_sep_k.long_name = 'precipitation'
			pcp_sep_k.var_name = 'pcp'
			pcp_sep_k.attributes['contributing_days'] = count_days_res[res]

			# Save
			outfile_k = outfile.replace('.nc', str(res)+'.nc')
			iris.save(pcp_sep_k, outfile_k)
			print outfile_k



def composite_pcp_tot_year(year1, year2, lead):
	"""Computes the total composite precipitation due to observed TCs for a particular
	TC season (July year1 - June year2)."""

	infiles = []
	for m in [7, 8, 9, 10, 11, 12]:  # 7,
		# m = str(m).zfill(2)
		filename = COMP_PCP_TOT_FILE % (lead, year1, m)
		infile = os.path.join(COMP_PCP_TOT_DIR, str(lead), str(year1), filename )
		infiles.append(infile)

	for m in [1, 2, 3, 4, 5, 6]:
		# m=str(m).zfill(2)
		# year=year2
		filename = COMP_PCP_TOT_FILE % (lead, year2, m)
		infile = os.path.join(COMP_PCP_TOT_DIR, str(lead), str(year2), filename)
		infiles.append(infile)


	print infiles

	outfile = "ukmo_nwp.comp_pcp_total.%d-days.%d-%d.nc" % (lead, year1, year2)
	outpath = os.path.join(COMP_PCP_TOT_DIR, str(lead), outfile)

	pl.add_files(infiles, outpath, cube_count=12)

def composite_pcp_tot_year_separate_resolutions(year1,year2, lead):
	"""Computes the total precipitation at a particular forecast lead time (in
	days) for a particular year."""

	# Check whether separation of resolutions is required
	#if year not in [2010, 2014]:
		#raise ValueError('All data in %d are at the same resolution' % year)

	# Get list of infiles
	infiles = []

	infiles = []
	for m in [7, 8, 9, 10, 11, 12]:  # 7,
		if (year1, m) == (2014,7):
			res = ['n512', 'n768']
			infile = [os.path.join(COMP_PCP_TOT_DIR, str(lead), str(year1), COMP_PCP_TOT_FILE % (lead, year1, m))] * 2
			infile = [i.replace('.nc', '%s.nc' % r) for i, r in zip(infile, res)]
			for a in infile:
				infiles.append(a)
				
		elif (year1, m) == (2017,7):
			res = ['n768', 'n1280']
			infile = [os.path.join(COMP_PCP_TOT_DIR, str(lead), str(year1), COMP_PCP_TOT_FILE % (lead, year1, m))] * 2
			infile = [i.replace('.nc', '%s.nc' % r) for i, r in zip(infile, res)]
			for a in infile:
				infiles.append(a)
		else:
			# m = str(m).zfill(2)
			filename = COMP_PCP_TOT_FILE % (lead, year1, m)
			infile = os.path.join(COMP_PCP_TOT_DIR, str(lead), str(year1), filename)
			infiles.append(infile)

	for m in [1, 2, 3, 4, 5, 6]:
		if (year2,m) == (2010,3):
			res = ['n320', 'n512']
			infile = [os.path.join(COMP_PCP_TOT_DIR, str(lead), str(year2), COMP_PCP_TOT_FILE % (lead, year2, m))] * 2
			infile = [i.replace('.nc', '%s.nc' % r) for i, r in zip(infile, res)]
			for a in infile:
				infiles.append(a)
		# m=str(m).zfill(2)
		# year=year2
		else:
			filename = COMP_PCP_TOT_FILE % (lead, year2, m)
			infile = os.path.join(COMP_PCP_TOT_DIR, str(lead), str(year2), filename)
			infiles.append(infile)

	print infiles

	outfile = "ukmo_nwp.comp_pcp_total.%d-days.%d-%d.nc" % (lead, year1, year2)
	outpath = os.path.join(COMP_PCP_TOT_DIR, str(lead), outfile)

	pl.add_files(infiles, outpath, separate_resolutions=True,deal_with_masks=False, cube_count=13)


#if lead == 6:
	#no lead time 7 days ahead for 2006-2008
	#year1s = [2010, 2011, 2012, 2013, 2015] #2006, 2007,  #will need to add 2008, 2009
	#year2s = [2011, 2012, 2013, 2014, 2016] #2007, 2008,  #will need to add 2009, 2010
#else:
	#year1s = [2006, 2007, 2008, 2010, 2011, 2012, 2013, 2015] #2006, 2007,  #will need to add 2008, 2009
	#year2s = [2007, 2008, 2009, 2011, 2012, 2013, 2014, 2016] #2007, 2008,  #will need to add 2009, 2010

#year1s_sep_res = [2009,2014]
#year2s_sep_res = [2010,2015]
#year1s = [2006]
#year2s = [2007]

year1s = [2016]
year2s = [2017]

year1s_sep_res = [2017]
year2s_sep_res = [2018]


#run monthly precip totals:
#for lead time 6, 2006 has no data - so can't run lead time 6 for 2006 or January 2007
for year in [2018]: #2011, 2012, 2013, 2014, 2015, 2016  

	for month in [9,10,11,12]: #,3,4,5,6,7,8,9,10,11,12
	
		composite_pcp_tot_year_month(year, month, lead)

#once 2018 issues are fixed, run months 2-12 above, and then run the following:
#run yearly (tc season) precip totals:
for year1, year2 in zip(year1s, year2s):
	composite_pcp_tot_year(year1, year2, lead)

for year1, year2 in zip(year1s_sep_res, year2s_sep_res):
	composite_pcp_tot_year_separate_resolutions(year1,year2,lead)




