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



REFORMAT_DIR = "/gws/nopw/j04/klingaman/emerton/ibtracs_reformatted_data_and_track_maps/"

INTERPOLATE_FILENAME = 'storms_SH{0}{1}_interpolated.txt'
 
DENSITY_DIR = "/gws/nopw/j04/klingaman/emerton/ibtracs_track_density/"

DENSITY_FILENAME = "ibtracs.density.{0}{1}.nc"


# To compute track density per year, first have to do per yearmonth
# Before running this, need a list of excluded days (missing NWP data dates to ignore)
# Calls exclude_days and track_density functions from picsea_library (well, currently from here)
def track_density_year_month(y, year1, year2, month, box_size=2.5, threeway=True):
	"""Calculates the TC track density as the number of tracks per grid
	box (n.b., NOT per unit time), with the given grid box size in degrees, at a
	particular forecast lead time (in days).

	y is the year the month in question belongs to, to get the right data and correctly name the output file
	year1 is the first year of the 2-year period, for the directory name
	year2 is the second year of the 2-year period, for the directory name"""

	#print lead
	# Check lead time is in available range (7 days)
	#if not 0 <= lead <= 6:
		#raise ValueError('lead=%s; must be 0 <= lead <= 6' % lead)


	# Compute track density

	#track_file = INTERPOLATE_FILENAME.format(y,m,d,h)
	infile = os.path.join(REFORMAT_DIR, INTERPOLATE_FILENAME.format(year1,year2))

	outfile = os.path.join(DENSITY_DIR, DENSITY_FILENAME.format(y,m)) 
	exclude = pl.exclude_days()

	#print exclude

	pl.track_density(infile, outfile, 'ibtracs',year=int(y),month=int(month), box_size=2.5, exclude=exclude, contributing_days=True)

def track_density_year(year1, year2, box_size=2.5, threeway=True):
	"""Sums the TC track density for the given lead time and year (TC season)."""

	# Check lead time is in available range
	#if not 0 <= lead <= 6:
		#raise ValueError('lead=%s; must be 0 <= lead <= 6' % lead)

	# Get list of files
	#dens_dir = DENS_DIR_3WAY if threeway else DENS_DIR
	#m2 = 12
	infiles=[]
	for m in [7,8,9,10,11,12]: #7,
		#m = str(m).zfill(2)
		infiles.append(DENSITY_FILENAME.format(year1,str(m).zfill(2)))

	for m in [1,2,3,4,5,6]:
		#m=str(m).zfill(2)
		#year=year2
		infiles.append(DENSITY_FILENAME.format(year2,str(m).zfill(2)))


	#infiles = [DENS_FILE % (lead, year, month) for month in xrange(1, m2+1)]
	inpaths = [os.path.join(DENSITY_DIR, ff) for ff in infiles]
	print inpaths

	outfile = "ibtracs.density.{0}{1}.nc".format(year1,year2)
	outpath = os.path.join(DENSITY_DIR, outfile)

	pl.add_files(inpaths, outpath, cube_count=12)


##################################################################################
#RUN FUNCTIONS
##################################################################################

#for each lead time, calculate the yearmonth track density for each month July - June the following year
#for l in range(7):
print "Calculating yearmonth track densities:"
for m in [7,8,9,10,11,12]: #7,
	m=str(m).zfill(2)
	print m
	y=y1
	print y
	track_density_year_month(y,y1,y2,m)

for m in [1,2,3,4,5,6]:
	m=str(m).zfill(2)
	print m
	y=y2
	print y
	track_density_year_month(y, y1, y2, m)

print "Calculating yearly track density, year ", y1,y2
track_density_year(y1,y2)
