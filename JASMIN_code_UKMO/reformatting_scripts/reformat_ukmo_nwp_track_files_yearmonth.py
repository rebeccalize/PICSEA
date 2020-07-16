import picsea_library as pl
import sys
import os
from calendar import monthrange

y=sys.argv[1]
y=int(y)
m=sys.argv[2]

UKMO_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis"
MATCH_DIR = os.path.join("/gws/nopw/j04/klingaman/emerton/TRACK/UKMO/Y{0}/UKMO_{1}{2}{3}{4}_VOR_VERTAVG_T63_DET/MATCH-ECMWF-ANALYSIS-IBT-3WAY/")
#MATCH_DIR = os.path.join("/group_workspaces/jasmin2/klingaman/emerton/TRACK/UKMO/Y%Y/UKMO_%Y%m%d%H_VOR_VERTAVG_T63_DET/MATCH-ECMWF-ANALYSIS-IBT-3WAY/")
REFORMAT_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_year/"
TRACK_FILE_PREFIX = 'trmatch_cntl_tr'
REFORMAT_FILENAME = 'ukmo_nwp.{0}{1}{2}{3}.tcident'
INTERPOLATE_FILENAME_2 = 'ukmo_nwp.{0}{1}{2}{3}.tcident.interp'
INTERPOLATE_FILENAME ='ukmo_nwp.%Y%m%d%H.tcident.interp'

def run_rewrite_track_files(y, m, d, h):
	"""Finds the files containing track info (after 3way matching), and rewrites in a more convenient format.
    In the resulting file, track_id and track_number will always be the same for these NWP forecasts"""

	# get the input directory with matched track files for this date
	indir = MATCH_DIR.format(y,y,m,d,h) #% (y, y, m, d, h)
	if os.path.exists(indir):

		# make a list of all the files that start with "trmatch_cntl_tr" - which means there is a matched track
		infiles = [f for f in os.listdir(indir) if f.startswith(TRACK_FILE_PREFIX)]
		infiles.sort()
		# create output filename for this date
		outfile = os.path.join(REFORMAT_DIR, str(y), REFORMAT_FILENAME.format(y,m,d,h)) # % (y, m, d, h)

		# run the function to rewrite these track files in a more convenient format for python
		pl.rewrite_list_of_track_files(indir, infiles, outfile, 2)


def run_interpolate_track_files(y, m, d, h):
	"""Takes the track files just generated using run_rewrite_track_files, and interpolates to halfway between
    the existing times."""

	indir = MATCH_DIR.format(y,y,m,d,h) #% (y, y, m, d, h)
	if os.path.exists(indir):

		infile = os.path.join(REFORMAT_DIR, str(y), REFORMAT_FILENAME.format(y,m,d,h)) # % (y, m, d, h)
		outfile = os.path.join(REFORMAT_DIR, str(y), INTERPOLATE_FILENAME_2.format(y,m,d,h))# % (y, m, d, h))

		pl.interpolate_track_file(infile, outfile)

#monthrange returns (weekday of first day of month, number of days in month)
#we want index [1], number of days in the month, to iterate over
ndays = monthrange(y,int(m))[1]
m=m.zfill(2)

for d in range(1,ndays+1):
	print "Reformatting ",y,m,d
	run_rewrite_track_files(y,m,str(d).zfill(2),str(0).zfill(2))
	run_rewrite_track_files(y, m, str(d).zfill(2), str(12))
	run_interpolate_track_files(y,m,str(d).zfill(2),str(0).zfill(2))
	run_interpolate_track_files(y, m, str(d).zfill(2), str(12))