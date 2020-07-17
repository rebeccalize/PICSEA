import sys
import os
import os.path
import math
import warnings
import itertools
import numpy as np
import datetime


#Reformat a file containing ibtracs track data so that it's easier for python to read
#Note: for the SH, track files run from July to July, with the filename the name of the start year

datadir = "/perm/mo/more/TIGGE/Y20182019/"
savedir = "/perm/mo/more/picsea/reformatted_idai_kenneth_tracks/"

def rewrite_track_file(infile, outfile, _source):
	"""Takes an input file containing track information from TRACK and rewrites
	it in a more convenient format."""

	# Check whether output file already exists
	# if os.path.isfile(outfile):
	# raise ValueError('Output file %s already exists' % outfile)

	# Open input file and create output file
	fin = open(infile)
	fout = open(outfile, 'w')

	# Write out header
	fout.write('track_id track_number point_number year month day hour '
			   'longitude latitude vorticity ???\n')

	# Iterate for each line of file
	if _source == 'jra55':
		track_number = 0
		last_genesis_year = 1958
	elif _source == 'ukmo_nwp':
		for skip in xrange(3):
			_ = fin.next()
	skip_next = False
	track_id = 0
	for line in fin:

		# Skip if necessary
		if skip_next:
			skip_next = False
			continue

		# Split line on delimiters
		data = [e for e in line.split() if e != '&']

		# If new TC, get track ID
		if data[0] == 'TRACK_ID':
			track_id += 1
			if _source == 'ukmo_nwp':
				track_number = int(data[1])
			point_number = 0
			skip_next = True
			continue

		# Get date and hour
		date = data[0]
		year = date[:4]
		month = date[4:6]
		day = date[6:8]
		hour = date[8:10]

		# Check whether this is new year
		if _source == 'jra55':
			if point_number == 0:
				if year != last_genesis_year:
					last_genesis_year = year
					track_number = 1
				else:
					track_number += 1

		# Get other data
		lon = data[1]
		lat = data[2]
		vort = data[3]
		others = ' '.join(data[4:])
		point_number += 1

		# Write out data
		write_out = (track_id, track_number, point_number, year, month, day,
					 hour, lon, lat, vort, others)
		fout.write('%04d %03d %03d %s %s %s %s %s %s %s %s\n' % write_out)

	# Close files
	fin.close()
	fout.close()
	print outfile
	
	
def interpolate_track_file(infile, outfile):
	"""Takes an input file containing track information from TRACK and
	linearly interpolates track points in time from [00,06,12,18] UTC to
	[03,09,15,21] UTC."""

	# Check whether output file already exists
	# if os.path.isfile(outfile):
	# raise ValueError('Output file %s already exists' % outfile)

	# Create output file
	fout = open(outfile, 'w')

	# Write out header
	fmt = '%04d %03d %03.01f %04d %02d %02d %02d %s %s\n'
	fout.write('track_id track_number point_number year month day hour '
			   'longitude latitude\n')

	# Read in data and get number of tracks
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		ain = np.genfromtxt(infile, dtype=float, skip_header=1,
							usecols=range(9))
	if not ain.size:
		fout.close()
		print outfile
		return
	ntrack = int(ain[-1, 0])
	nrows = ain.shape[0]

	# Iterate for each track
	dt_func = np.vectorize(lambda y, m, d, h: datetime.datetime(y, m, d, h))
	dt_inv_func = np.vectorize(lambda dt, att: getattr(dt, att))
	atts = ['year', 'month', 'day', 'hour']
	count_rows = 0
	for itrack in xrange(1, ntrack + 1):

		# Get data for this track
		ain_t = ain[ain[:, 0] == itrack]
		npoint = int(ain_t[-1, 2])

		# Deal with columns which can safely be linearly interpolated
		# (i.e., everything except date/time)
		ain_t_to_interp = ain_t[:, [0, 1, 2, 7, 8]]
		ain_t_interp0 = (ain_t_to_interp[:-1] + ain_t_to_interp[1:]) / 2.

		# Interpolate date/times
		ain_t_datetime = ain_t[:, 3:7].astype(int)
		datetimes = dt_func(*[ain_t_datetime[:, i] for i in xrange(4)])
		datetimes_interp = datetimes[:-1] + (datetimes[1:] - datetimes[:-1]) / 2
		y_interp, m_interp, d_interp, h_interp = \
			[dt_inv_func(datetimes_interp, att).reshape((npoint - 1, 1)) for att
			 in atts]

		# Create array for all interpolated data for this track
		ain_t_interp = np.concatenate((ain_t_interp0[:, :3], y_interp,
									   m_interp, d_interp, h_interp, ain_t_interp0[:, 3:]), axis=1)

		# Write to file
		for ll in ain_t_interp:
			fout.write(fmt % tuple(ll))
			count_rows += 1

	# Sanity check
	if count_rows != nrows - ntrack:
		raise ValueError('Number of rows of output data is %d; should be %d' % \
						 (count_rows, nrows - ntrack))

	# Close file
	fout.close()
	print outfile





#dates = [2019030100,2019030112,2019030200,2019030212,2019030300,2019030312,2019030400,2019030412, 2019030500, 2019030512,2019030600,2019030612,2019030700,2019030712,2019030800,2019030812,2019030900,2019030912,2019031000,2019031012,2019031100,2019031112,2019031200,2019031212,2019031300,2019031312,2019031400, 2019031412]	

for d in range(1,31): #32

	for h in [0,12]:

		day = str(d).zfill(2)
		print day
	
		hr = str(h).zfill(2)
	
		infile_eps = datadir+"201902"+day+hr+"/MATCH-ECMWF-ANALYSIS-IBT-3WAY/trmatch_cntl_tr0001"
		infile_mean = datadir+"201902"+day+hr+"/MATCH-ECMWF-ANALYSIS-IBT-3WAY/trmatch_cntl_tr0001_mean"
		infile_det = datadir+"201902"+day+hr+"/MATCH-ECMWF-ANALYSIS-IBT-3WAY/DET/trmatch_cntl_tr0001"
	
		outfile_eps = savedir+"201902"+day+hr+"_idai_ECMWF_eps_3way_matched_reformatted.txt"
		outfile_mean = savedir+"201902"+day+hr+"_idai_ECMWF_mean_3way_matched_reformatted.txt"
		outfile_det = savedir+"201902"+day+hr+"_idai_ECMWF_det_3way_matched_reformatted.txt"
		
		
		if os.path.isfile(infile_eps):		
			rewrite_track_file(infile_eps,outfile_eps,'ukmo_nwp')
			
		else:
			continue
			
			
		if os.path.isfile(infile_mean):
			rewrite_track_file(infile_mean,outfile_mean,'ukmo_nwp')
		else:
			continue
			
			
		if os.path.isfile(infile_det):
			rewrite_track_file(infile_det,outfile_det,'ukmo_nwp')
		else:
			continue
		
		
	
		infile_eps = savedir+"201902"+day+hr+"_idai_ECMWF_eps_3way_matched_reformatted.txt"
		infile_mean = savedir+"201902"+day+hr+"_idai_ECMWF_mean_3way_matched_reformatted.txt"
		infile_det = savedir+"201902"+day+hr+"_idai_ECMWF_det_3way_matched_reformatted.txt"
		
		outfile_eps = savedir+"201902"+day+hr+"_idai_ECMWF_eps_3way_matched_interpolated.txt"
		outfile_mean = savedir+"201902"+day+hr+"_idai_ECMWF_mean_3way_matched_interpolated.txt"
		outfile_det = savedir+"201902"+day+hr+"_idai_ECMWF_det_3way_matched_interpolated.txt"
		
		
		if os.path.isfile(infile_eps):		
			interpolate_track_file(infile_eps,outfile_eps)
		else:
			continue
			
			
		if os.path.isfile(infile_mean):
			interpolate_track_file(infile_mean,outfile_mean)
		else:
			continue
			
			
		if os.path.isfile(infile_det):
			interpolate_track_file(infile_det,outfile_det)
		else:
			continue
		
		
	
		
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

