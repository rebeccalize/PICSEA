import picsea_library as pl
import composite_functions_ukmo_nwp as cf
import sys
import os
import warnings
import itertools
import numpy as np
import iris
import datetime
from calendar import monthrange
import cf_units
from iris.time import PartialDateTime
import iris.coord_categorisation

#

#Require 2 years, because for the SH we are doing analysis from July - June across 2 years
y1=int(sys.argv[1])
y2=int(sys.argv[2])


UKMO_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis"
MATCH_DIR = os.path.join("/gws/nopw/j04/klingaman/emerton/TRACK/UKMO/Y{0}/UKMO_{1}{2}{3}{4}_VOR_VERTAVG_T63_DET/MATCH-ECMWF-ANALYSIS-IBT-3WAY/")
#MATCH_DIR = os.path.join("/group_workspaces/jasmin2/klingaman/emerton/TRACK/UKMO/Y%Y/UKMO_%Y%m%d%H_VOR_VERTAVG_T63_DET/MATCH-ECMWF-ANALYSIS-IBT-3WAY/")
REFORMAT_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_year/"
TRACK_FILE_PREFIX = 'trmatch_cntl_tr'
REFORMAT_FILENAME = 'ukmo_nwp.{0}{1}{2}{3}.tcident'
INTERPOLATE_FILENAME_2 = 'ukmo_nwp.{0}{1}{2}{3}.tcident.interp'
INTERPOLATE_FILENAME ='ukmo_nwp.%Y%m%d%H.tcident.interp'
DENSITY_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/track_density/"
DENSITY_FILENAME = "ukmo_nwp.density.lt{0}_days.{1}{2}.nc"


# To compute track density per year, first have to do per yearmonth
# Before running this, need a list of excluded days (missing NWP data dates to ignore)
# Calls exclude_days and track_density functions from picsea_library (well, currently from here)
def track_density_year_month(lead, y, year1, year2, month, box_size=2.5, threeway=True):
	"""Calculates the TC track density as the number of tracks per grid
	box (n.b., NOT per unit time), with the given grid box size in degrees, at a
	particular forecast lead time (in days).

	y is the year the month in question belongs to, to get the right data and correctly name the output file
	year1 is the first year of the 2-year period, for the directory name
	year2 is the second year of the 2-year period, for the directory name"""

	#print lead
	# Check lead time is in available range (7 days)
	if not 0 <= lead <= 6:
		raise ValueError('lead=%s; must be 0 <= lead <= 6' % lead)


	# Compute track density

	#track_file = INTERPOLATE_FILENAME.format(y,m,d,h)
	infile = os.path.join(REFORMAT_DIR, '%Y', INTERPOLATE_FILENAME)

	outfile = os.path.join(DENSITY_DIR, str(lead), str(year1)+str(year2), DENSITY_FILENAME.format(l,y,m)) # % (lead, year, month)
	exclude = pl.exclude_days(lead=lead)

	#print exclude

	pl.track_density(infile, outfile, 'ukmo_nwp',year=int(y),month=int(month), box_size=2.5, lead=int(lead), exclude=exclude, contributing_days=True)

def track_density_year(lead, year1, year2, box_size=2.5, threeway=True):
	"""Sums the TC track density for the given lead time and year (TC season)."""

	# Check lead time is in available range
	if not 0 <= lead <= 6:
		raise ValueError('lead=%s; must be 0 <= lead <= 6' % lead)

	# Get list of files
	#dens_dir = DENS_DIR_3WAY if threeway else DENS_DIR
	#m2 = 12
	infiles=[]
	for m in [7,8,9,10,11,12]: #7,
		#m = str(m).zfill(2)
		infiles.append(DENSITY_FILENAME.format(lead,year1,str(m).zfill(2)))

	for m in [1,2,3,4,5,6]:
		#m=str(m).zfill(2)
		#year=year2
		infiles.append(DENSITY_FILENAME.format(lead,year2,str(m).zfill(2)))


	#infiles = [DENS_FILE % (lead, year, month) for month in xrange(1, m2+1)]
	inpaths = [os.path.join(DENSITY_DIR, str(lead), str(year1)+str(year2), ff) for ff in infiles]
	print inpaths

	outfile = "ukmo_nwp.density.lt{0}_days.{1}{2}.nc".format(lead,year1,year2)
	outpath = os.path.join(DENSITY_DIR, str(lead), str(year1)+str(year2), outfile)

	pl.add_files(inpaths, outpath, cube_count=12)


##################################################################################
#RUN FUNCTIONS
##################################################################################

#for each lead time, calculate the yearmonth track density for each month July - June the following year
for l in range(7):
	print "Calculating yearmonth track densities, lead time ",l
	for m in [7,8,9,10,11,12]: #7,
		m=str(m).zfill(2)
		print m
		y=y1
		print y
		track_density_year_month(l,y,y1,y2,m)

	for m in [1,2,3,4,5,6]:
		m=str(m).zfill(2)
		print m
		y=y2
		print y
		track_density_year_month(l, y, y1, y2, m)

	print "Calculating yearly track density, lead time", l, "year ", y1,y2
	track_density_year(l,y1,y2)
