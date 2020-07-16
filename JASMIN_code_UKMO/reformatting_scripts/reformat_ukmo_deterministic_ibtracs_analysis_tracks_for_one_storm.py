import sys
sys.path.insert(1, '/home/users/emerton/analysis_scripts')
import picsea_library as pl
import os
from calendar import monthrange

y1=int(sys.argv[1])
y2=int(sys.argv[2])
track=sys.argv[3]
#m=sys.argv[2]


UKMO_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis"

MATCH_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/{}_{}/tr{}/"

REFORMAT_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/{}_{}/tr{}/"

TRACK_FILE_PREFIX = 'trmatch_cntl_tr'
ANALYSIS_REFORMAT_FILENAME = 'analysis_tr{}_reformatted.txt' #track number
IBTRACS_REFORMAT_FILENAME = 'ibtracs_tr{}_reformatted.txt'
FORECAST_REFORMAT_FILENAME = 'ukmo_nwp_{}_tr{}_reformatted.txt'



def rewrite_track_file_rebecca(indir, infiles, outfile, find_id):
	"""Takes several input files containing track information from TRACK and rewrites
	them in a more convenient format.

	find_id tells you which track to use from the 3way matched file, which contains the
	analysis (0), ibtracs (1) and forecast (2) tracks (e.g. find_id = 2 is the forecast track)"""

	# Check whether output file already exists
	# if os.path.isfile(outfile):
	# raise ValueError('Output file %s already exists' % outfile)

	# Create output directory if necessary
	#odir = os.path.split(outfile)[0]
	#if not os.path.isdir(odir):
	#	os.makedirs(odir)

	# Create output file and write out header
	fout = open(outfile, 'w')
	header = 'track_id track_number point_number year month day hour ' + \
			 'longitude latitude mslp?? 10mwind/vorticity??\n'
	fout.write(header)

	# Iterate for each input file
	track_id = 0
	print "INFILES = ", infiles
	for infile in infiles:

		# Open file
		print infile
		fin = open(infile)
		track_number = int(track)

		# Skip as far as 'TRACK_ID  <find_id>' (forecast track) then read rest
		# of file
		l = fin.readline()
		while not l.startswith('TRACK_ID  %d' % find_id):
			l = fin.readline()
		l = fin.readline()  # Skip one more line
		rl = fin.readlines()
		track_id += 1

		# Iterate for each line
		point_number = 0
		skip_next = False
		use_track = True
		for line in rl:

			# Skip if necessary
			if skip_next:
				skip_next = False
				continue

			# Split line on delimiters
			data = [e for e in line.split() if e != '&']

			# If new TC, increment track ID
			if data[0] == 'TRACK_ID':
				track_id += 1
				point_number = 0
				skip_next = True
				use_track = int(data[1]) == find_id
				continue

			# Skip if wrong TRACK_ID
			if not use_track:
				continue

			# Get date and hour
			date = data[0]
			year = date[:4]
			month = date[4:6]
			day = date[6:8]
			hour = date[8:10]

			# Get other data
			lon = data[1]
			lat = data[2]
			vort = data[3]
			others = ' '.join(data[4:])
			point_number += 1

			# Write out data
			write_out = (track_id, track_number, point_number, year, month,
						 day, hour, lon, lat, vort, others)
			fout.write('%04d %03d %03d %s %s %s %s %s %s %s %s\n' % write_out)

		# Close input file
		fin.close()

	# Close output file
	fout.close()
	print outfile

def run_rewrite_track_files(date, track, which_track):
	"""Finds the files containing track info (after 3way matching), and rewrites in a more convenient format.
    In the resulting file, track_id and track_number will always be the same for these NWP forecasts"""

	# get the input directory with matched track files for this date
	indir = MATCH_DIR.format(y1,y2,track) #% (y, y, m, d, h)
	#indir = MATCH_DIR % (y, y, m, d, h)
	# make a list of all the files that start with "trmatch_cntl_tr" - which means there is a matched track
	infiles = [indir+date+"_tr"+track+".txt"]
	print infiles
	
	print track
	print y1
	print y2
	
	if which_track == "analysis":
	    find_id = 0
	    outfile = os.path.join(REFORMAT_DIR.format(y1,y2,track), ANALYSIS_REFORMAT_FILENAME.format(track)) # % (y, m, d, h)
	    
	elif which_track == "ibtracs":
	    find_id = 1
	    outfile = os.path.join(REFORMAT_DIR.format(y1,y2,track), IBTRACS_REFORMAT_FILENAME.format(track))
	    
	elif which_track == "ukmo_nwp":
	    find_id = 2
	    outfile = os.path.join(REFORMAT_DIR.format(y1,y2,track), FORECAST_REFORMAT_FILENAME.format(date,track))
	    
	else:
	    print "something went wrong specifying which track to reformat"

	# run the function to rewrite these track files in a more convenient format for python
	rewrite_track_file_rebecca(indir, infiles, outfile, find_id)


def run_interpolate_track_files(y, m, d, h):
	"""Takes the track files just generated using run_rewrite_track_files, and interpolates to halfway between
    the existing times."""

	infile = os.path.join(REFORMAT_DIR, str(y), REFORMAT_FILENAME.format(y,m,d,h)) # % (y, m, d, h)
	outfile = os.path.join(REFORMAT_DIR, str(y), INTERPOLATE_FILENAME_2.format(y,m,d,h))# % (y, m, d, h))

	pl.interpolate_track_file(infile, outfile)

#get the dates of all the forecasts for this storm
dates=[]

search_dir=MATCH_DIR.format(y1,y2,track)
print search_dir
for filename in os.listdir(search_dir):
    print filename
    date,other = filename.split("_")
    print date
    dates.append(date)

#loop over all the dates and reformat the track files
for d,i in zip(dates,range(len(dates))):
    #only want to run the analysis and ibtracs once, and the forecast for each date
    if i == 0: 
        run_rewrite_track_files(d,track,"analysis")
        run_rewrite_track_files(d,track,"ibtracs")
        run_rewrite_track_files(d,track,"ukmo_nwp")
    elif i > 0:
        run_rewrite_track_files(d,track,"ukmo_nwp")


