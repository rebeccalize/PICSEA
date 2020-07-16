import picsea_library as pl
import sys
import os
from calendar import monthrange

y1=int(sys.argv[1])
y2=int(sys.argv[2])
track=sys.argv[3]
#m=sys.argv[2]

#this script should work regardless of how many ensemble members there were at a given time


#UKMO_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_TIGGE_analysis"

datadir = "/gws/nopw/j04/klingaman/emerton/ukmo_TIGGE_analysis/reformatted_track_files_per_storm/Y{}{}/SIO_storms/tr{}/"

savedir = "/gws/nopw/j04/klingaman/emerton/ukmo_TIGGE_analysis/reformatted_track_files_per_storm/Y{}{}/SIO_storms/tr{}/"

TRACK_FILE_PREFIX = 'trmatch_cntl_tr'
ANALYSIS_REFORMAT_FILENAME = 'analysis_tr{}_reformatted.txt' #track number
IBTRACS_REFORMAT_FILENAME = 'ibtracs_tr{}_reformatted.txt'
MEAN_REFORMAT_FILENAME = 'ukmo_tigge_mean_{}_tr{}_reformatted.txt'
CTRL_REFORMAT_FILENAME = 'ukmo_tigge_ctrl_{}_tr{}_reformatted.txt'
ENS_REFORMAT_FILENAME = 'ukmo_tigge_ens_{}_tr{}_reformatted.txt'


def rewrite_track_file_singletrack(indir, infiles, outfile, find_id):
	"""Takes several input files containing track information from TRACK and rewrites
	them in a more convenient format.

	find_id tells you which track to use from the 3way matched file, which, for files where tracking was run then MSLP & winds added separately, contains:
	
	in the ensemble forecast file: analysis (0), ibtracs (1), control (2) ensembles (3-22) #for 2010-2014, or 13 for 2014-2017, or 
	in the mean forecast file: ensemble mean (0)"""
	
	

	# Check whether output file already exists
	# if os.path.isfile(outfile):
	# raise ValueError('Output file %s already exists' % outfile)

	# Create output directory if necessary
	#odir = os.path.split(outfile)[0]
	#if not os.path.isdir(odir):
	#	os.makedirs(odir)

	# Create output file and write out header
	fout = open(outfile, 'w')
	header = 'track_id track_number point_number year month day hour longitude latitude mslp?? 10mwind/vorticity??\n'
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
	
	
	

def rewrite_track_file_ensemble(indir, infiles, outfile):
	"""Takes several input files containing track information from TRACK and rewrites
	them in a more convenient format.

	find_id tells you which track to use from the 3way matched file, which contains:
	
	in the ensemble forecast file: analysis (0), ibtracs (1), control (2) ensembles (3-22) #or 13, or something else
	in the mean forecast file: ensemble mean (0) """

	# Check whether output file already exists
	# if os.path.isfile(outfile):
	# raise ValueError('Output file %s already exists' % outfile)

	# Create output directory if necessary
	#odir = os.path.split(outfile)[0]
	#if not os.path.isdir(odir):
	#	os.makedirs(odir)

	# Create output file and write out header
	fout = open(outfile, 'w')
	header = 'track_id track_number point_number year month day hour longitude latitude mslp?? 10mwind/vorticity??\n'
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
		while not l.startswith('TRACK_ID'):
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
				track_id = int(data[1])
				point_number = 0
				skip_next = True
				#use_track = int(data[1]) >= find_id
				continue

			# Skip if wrong TRACK_ID
			#if not use_track:
				#continue

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

			if track_id > 2: #(don't want to include the analysis, ibtracs or control,
			#but no way to guarantee the control isn't in there using find_ids, as some ensembles
			# have no TC in the foreacst..., so numbers not usually/always consecutive from 3 to 52...)
				

				# Write out data
				write_out = (track_id, track_number, point_number, year, month, day, hour, lon, lat, vort, others)
				fout.write('%04d %03d %03d %s %s %s %s %s %s %s %s\n' % write_out)
				
			else:
				continue

		# Close input file
		fin.close()

	# Close output file
	fout.close()
	print outfile
	
	
	
	

def run_rewrite_track_files(date, track, which_track):
	"""Finds the files containing track info (after 3way matching), and rewrites in a more convenient format.
   	In the resulting file, track_id and track_number will always be the same for these NWP forecasts
    
    
    	find_id tells you which track to use from the 3way matched file, which contains:
	
	in the ensemble forecast file: analysis (0), ibtracs (1), control (2) ensembles (3-22) #or 13, or something else
	in the mean forecast file: ensemble mean (0)  """

	# get the input directory with matched track files for this date
	indir = datadir.format(y1,y2,track) #% (y, y, m, d, h)
	#indir = MATCH_DIR % (y, y, m, d, h)
	# make a list of all the files that start with "trmatch_cntl_tr" - which means there is a matched track
	#infiles = [indir+date+"_tr"+track+".txt"]
	#print infiles
	
	print track
	print y1
	print y2
	
	if which_track == "analysis":
	    find_id = 0
	    outfile = os.path.join(savedir.format(y1,y2,track), ANALYSIS_REFORMAT_FILENAME.format(track)) # % (y, m, d, h)
	    
	    infiles = [indir+date+"_tr"+track+".txt"]
	    print infiles   
	    
	    rewrite_track_file_singletrack(indir, infiles, outfile, find_id) 
	    
	    
	elif which_track == "ibtracs":
	    find_id = 1
	    outfile = os.path.join(savedir.format(y1,y2,track), IBTRACS_REFORMAT_FILENAME.format(track))
	    
	    infiles = [indir+date+"_tr"+track+".txt"]
	    print infiles 
	    
	    rewrite_track_file_singletrack(indir, infiles, outfile, find_id)
	    
	    
	elif which_track == "mean":
	    find_id = 0
	    outfile = os.path.join(savedir.format(y1,y2,track), MEAN_REFORMAT_FILENAME.format(date,track))
	    
	    infiles = [indir+date+"_tr"+track+"_mean.txt"]
	    print infiles 
	    
	    rewrite_track_file_singletrack(indir, infiles, outfile, find_id)
	    
	elif which_track == "control":
	    find_id = 2
	    outfile = os.path.join(savedir.format(y1,y2,track), CTRL_REFORMAT_FILENAME.format(date,track))
	    
	    infiles = [indir+date+"_tr"+track+".txt"]
	    print infiles 
	    
	    rewrite_track_file_singletrack(indir, infiles, outfile, find_id)
	    
	    
	elif which_track == "ensemble":
	    find_id = 3
	    outfile = os.path.join(savedir.format(y1,y2,track), ENS_REFORMAT_FILENAME.format(date,track))
	    
	    infiles = [indir+date+"_tr"+track+".txt"]
	    print infiles 
	    
	    rewrite_track_file_ensemble(indir, infiles, outfile)
	    
	    
	    
	else:
	    print "something went wrong specifying which track to reformat"

	# run the function to rewrite these track files in a more convenient format for python
#	rewrite_track_file_rebecca(indir, infiles, outfile, find_id)
	
	
	


#def run_interpolate_track_files(y, m, d, h):
	"""Takes the track files just generated using run_rewrite_track_files, and interpolates to halfway between
    the existing times."""

	#infile = os.path.join(savedir, str(y), REFORMAT_FILENAME.format(y,m,d,h)) # % (y, m, d, h)
	#outfile = os.path.join(savedir, str(y), INTERPOLATE_FILENAME_2.format(y,m,d,h))# % (y, m, d, h))

	#pl.interpolate_track_file(infile, outfile)
	
	
	
#calling the function a couple of times, as some dates have either the mean or the det, or the ensemble, file missing	

	
#get the dates of all the forecasts for this storm
dates=[]

search_dir=datadir.format(y1,y2,track)
print search_dir
for filename in os.listdir(search_dir):
    if filename.endswith("_mean.txt"):
    	print filename
    	date,other,fcst = filename.split("_")
    	print date
   	dates.append(date)

#loop over all the dates and reformat the track files
for d,i in zip(dates,range(len(dates))):

	run_rewrite_track_files(d,track,"mean")




#get the dates of all the forecasts for this storm
dates=[]

search_dir=datadir.format(y1,y2,track)
print search_dir
for filename in os.listdir(search_dir):
    if filename.endswith(track+".txt"):
    	print filename
    	date,other = filename.split("_")
    	print date
   	dates.append(date)

#loop over all the dates and reformat the track files
for d,i in zip(dates,range(len(dates))):

	if i == 0:
		print d
		print i
		run_rewrite_track_files(d,track,"analysis")
		run_rewrite_track_files(d,track,"ibtracs")

		run_rewrite_track_files(d,track,"control")
		run_rewrite_track_files(d,track,"ensemble")
		
	elif i > 0:
		print d
		print i
		run_rewrite_track_files(d,track,"control")
		run_rewrite_track_files(d,track,"ensemble")

