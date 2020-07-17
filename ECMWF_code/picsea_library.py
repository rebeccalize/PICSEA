# A collection of functions used for PICSEA
# Use "import picsea_library as pl" from other scripts...

import sys
import os
import math
import warnings
import itertools
import numpy as np
#import iris
import datetime
#import iris.coord_categorisation
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.cm as cm
import numpy as np
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import matplotlib.colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl

# DIRECTORIES & FILENAMES
TRMM_DIR = "/gws/nopw/j04/klingaman/datasets/TRMM_3B42/V7_NC/"
TRMM_FILE = "3B42.%d%02d%02d.%02d.7.nc"
UKMO_DIR = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis"


# Functions included:

# round_up_to_even
# get_dates
# exclude_days
# add_files
# extrange
# rewrite_track_file
# rewrite_list_of_track_files
# interpolate_track_file

# map_tracks_one_year
# map_tracks_multi_year
# map_nwp_tracks_per_storm
# map_composite_data

# read_trmm_3b42

# add_cubes
# bounds

# track_density

def round_up_to_even(f):
	return math.ceil(f / 2.) * 2

#################################################################################
# UKMO NWP HELPER FUNCTIONS
#################################################################################

def get_dates(data_from_track_file):

	data = data_from_track_file
	year = data[:, 3]
	month = data[:, 4]
	day = data[:, 5]
	hh = data[:, 6]

	# join the year, month, day and hour together to get the whole date
	datelist = []
	for i in range(len(year)):
		date = []
		date.append(str(int(year[i])))
		date.append(str(int(month[i])).zfill(2))
		date.append(str(int(day[i])).zfill(2))
		date.append(str(int(hh[i])).zfill(2))
		wholedate = ''.join(date)
		datelist.append(wholedate)

	return datelist


#################################################################################


def exclude_days(lead=None):
	"""Returns list of datetime objects of days to exclude at the given lead
	time."""

	# Get list of datetime objects of missing forecast (T+0) times
	strs = np.genfromtxt(os.path.join(UKMO_DIR, 'IDAI/MISSING_IDAI.txt'), comments='#', delimiter=13, usecols=0, dtype=str)
	dt_l = lambda ss: datetime.datetime.strptime(ss, '%Y/%m/%d %H')
	#dt_l = lambda ss: datetime.datetime.strptime(ss, '%Y/%m/%d/%H')
	dt_func = np.vectorize(dt_l)
	dts = dt_func(strs)

	# If no lead time, just return list of all datetimes
	if lead is None:
		return dts

	# Get range of validity times which should be excluded
	t1 = dts + datetime.timedelta(hours=24 * lead)
	t2 = dts + datetime.timedelta(hours=24 * lead + 21)
	t_excl_all = np.concatenate((t1, t2))
	t_excl_all.sort()

	# Convert all times to midnight and remove duplicates
	hour0_func = np.vectorize(lambda tt: tt - datetime.timedelta(hours=tt.hour))
	return np.unique(hour0_func(t_excl_all))


#################################################################################

#from iris.experimental.equalise_cubes import equalise_attributes

def add_files(infiles, outfile, separate_resolutions=False,
			  deal_with_masks=True, cube_count=None, contributing_days=True):
	"""Adds together the given files.

	If separate_resolutions is True, adds the resolution name (e.g., n768) onto
	the end of the file name.
	"""

	# Check whether output file already exists
	if os.path.isfile(outfile) and not separate_resolutions:
		raise ValueError('Output file %s already exists' % outfile)

	# Sum data and save
	cubes = iris.load_raw(infiles)
	for i in range(len(cubes)):
		print '----------------'
		print cubes[i].coords
		print '----------------'
	if cube_count is not None:
		if len(cubes) != cube_count:
			raise ValueError('%d Cubes found, not %d' % (len(cubes), cube_count))
	total_cubes = add_cubes(cubes, separate_resolutions=separate_resolutions,
							 deal_with_masks=deal_with_masks,
							 contributing_days=contributing_days)
	if separate_resolutions:
		for k in total_cubes.iterkeys():
			res = {640: 'n320', 1024: 'n512', 1536: 'n768'}[k[1]]
			outfile_k = outfile.replace('.nc', '.%s.nc' % res)
			if os.path.isfile(outfile_k):
				raise ValueError('Output file %s already exists' % outfile_k)
			iris.save(total_cubes[k], outfile_k)
			print outfile_k
	else:
		iris.save(total_cubes, outfile)
		print outfile


#################################################################################

def extrange(*args):
	"""extrange([start,] stop[, step]) -> generator object

	Extends the built-in xrange() function by allowing any data types,
	provided *step* can be added to *start*, and *start* can be compared
	to *stop* using >=.

	For example, *start* and *stop* could be `datetime.datetime`s and
	*step* a `datetime.timedelta`.

	The end point is omitted!

	*start* and *step* may each be omitted but will default to integers
	0 and 1 respectively, which may be inappropriate depending on the
	type(s) of the other argument(s) given.
	"""

	# Parse arguments
	msg = None
	if not args:
		msg = 'extrange expected at least 1 arguments, got 0'
	elif len(args) == 1:
		start, stop, step = 0, args[0], 1
	elif len(args) == 2:
		(start, stop), step = args, 1
	elif len(args) == 3:
		start, stop, step = args
	else:
		msg = 'range expected at most 3 arguments, got %d' % len(args)
	if msg is not None:
		raise TypeError(msg)

	# Check for compatible argument types
	_ = start + step
	_ = start >= stop

	# Generate values
	while True:
		if start >= stop:
			break
		yield start
		start += step


#################################################################################
# TRACK REFORMATTING FUNCTIONS
#################################################################################

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


#################################################################################

def rewrite_list_of_track_files(indir, infiles, outfile, find_id):
	"""Takes several input files containing track information from TRACK and rewrites
	them in a more convenient format.

	find_id tells you which track to use from the 3way matched file, which contains the
	analysis, forecast and ibtracs tracks (e.g. find_id = 2 is the forecast track"""

	# Check whether output file already exists
	# if os.path.isfile(outfile):
	# raise ValueError('Output file %s already exists' % outfile)

	# Create output directory if necessary
	odir = os.path.split(outfile)[0]
	if not os.path.isdir(odir):
		os.makedirs(odir)

	# Create output file and write out header
	fout = open(outfile, 'w')
	header = 'track_id track_number point_number year month day hour ' + \
			 'longitude latitude vorticity ???\n'
	fout.write(header)

	# Iterate for each input file
	track_id = 0
	for infile in infiles:

		# Open file
		fin = open(os.path.join(indir, infile))
		track_number = int(infile[-4:])

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


#################################################################################    

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


#################################################################################    
# MAPPING FUNCTIONS
#################################################################################

def draw_rectangle( lats, lons, m):
	x, y = m( lons, lats )
	xy = zip(x,y)
	poly = Polygon(xy, facecolor='None', edgecolor='darkgrey', linewidth=1,alpha=0.75)
	plt.gca().add_patch(poly)

def map_tracks_one_year(infile, outfile, region):
	"""Plots a map of storm tracks in the Southern Hemisphere
	from a TRACK file that has been reformatted to suit Python
	Plots the tracks so that they change colour from lighter to darker with time
	Note: the colormap used is hardcoded in the function (where 'colors' is assigned)"""

	# load in data from a specified TRACK file (reformatted & interpolated)
	data = np.genfromtxt(infile, dtype=float, skip_header=1)

	# set up map of specified region
	if region == "SH":
		lat1 = 30
		lat2 = -60
		lon1 = -25
		lon2 = 335

	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')

	m.drawcoastlines(linewidth=0.4, color='lightgray')
	m.drawcountries(linewidth=0.4, color='lightgray')

	no_tracks = int(np.max(data[:, 0]))  # finds the max number in the track_id column i.e. total no. tracks in file

	# iterate over each track in the file
	for i_track in range(1, no_tracks + 1):
		# get the numbers of the rows in the file that belong to this track (no. points can vary per storm)
		data_indices = np.where(data[:, 0] == i_track)
		# count how many track points this storm has
		no_points = len(data_indices[0])
		# set up an array to hold the lat and lon of each point of this storm
		track_coords = np.zeros((no_points, 2))
		# get the lat and lons from the track file for this storm only
		track_coords[:, 0] = data[np.min(data_indices[0]):np.max(data_indices[0]) + 1, 7]  # longitudes
		track_coords[:, 1] = data[np.min(data_indices[0]):np.max(data_indices[0]) + 1, 8]  # latitudes

		# To have the line changing colour from start to finish, split the coordinates into segments
		points = np.array([track_coords[:, 0], track_coords[:, 1]]).T.reshape(-1, 1, 2)
		segments = np.concatenate([points[:-1], points[1:]], axis=1)

		colors = cm.BuPu(np.linspace(0.1, 1, no_points))

		# iterate over each segment of this track
		for p, cl in zip(range(no_points - 1), colors):
			# need to separate the x and y values for each segment, and plot each segment on the map
			xarr = []
			yarr = []
			xarr.append(segments[p][0][0])
			xarr.append(segments[p][1][0])
			yarr.append(segments[p][0][1])
			yarr.append(segments[p][1][1])
			x, y = m(xarr, yarr)
			m.plot(x, y, linewidth=0.5, color=cl)

	fig.subplots_adjust(wspace=0)
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.05, dpi=300)


#################################################################################


def map_tracks_multi_year(infile_list, outfile, region):
	"""Plots a map of storm tracks in the Southern Hemisphere
	from a list of several TRACK files (e.g. one per year) that has been reformatted to suit Python
	Plots the tracks so that they change colour from lighter to darker with time
	Note: the colormap used is hardcoded in the function (where 'colors' is assigned)
	Function takes a list of filenames as  input argument"""

	# set up map of region
	if region == "SH":
		lat1 = 30
		lat2 = -60
		lon1 = -25
		lon2 = 335

	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')

	m.drawcoastlines(linewidth=0.4, color='lightgray')
	m.drawcountries(linewidth=0.4, color='lightgray')

	for file_in_list in infile_list:

		# load in data from a specified TRACK file (reformatted & interpolated)
		data = np.genfromtxt(file_in_list, dtype=float, skip_header=1)
		no_tracks = int(np.max(data[:, 0]))  # finds the max number in the track_id column i.e. total no. tracks in file

		# iterate over each track in the file
		for i_track in range(1, no_tracks + 1):
			# get the numbers of the rows in the file that belong to this track (no. points can vary per storm)
			data_indices = np.where(data[:, 0] == i_track)
			# count how many track points this storm has
			no_points = len(data_indices[0])
			# set up an array to hold the lat and lon of each point of this storm
			track_coords = np.zeros((no_points, 2))
			# get the lat and lons from the track file for this storm only
			track_coords[:, 0] = data[np.min(data_indices[0]):np.max(data_indices[0]) + 1, 7]  # longitudes
			track_coords[:, 1] = data[np.min(data_indices[0]):np.max(data_indices[0]) + 1, 8]  # latitudes

			# To have the line changing colour from start to finish, split the coordinates into segments
			points = np.array([track_coords[:, 0], track_coords[:, 1]]).T.reshape(-1, 1, 2)
			segments = np.concatenate([points[:-1], points[1:]], axis=1)

			colors = cm.BuPu(np.linspace(0.1, 1, no_points))

			# iterate over each segment of this track
			for p, cl in zip(range(no_points - 1), colors):
				# need to separate the x and y values for each segment, and plot each segment on the map
				xarr = []
				yarr = []
				xarr.append(segments[p][0][0])
				xarr.append(segments[p][1][0])
				yarr.append(segments[p][0][1])
				yarr.append(segments[p][1][1])
				x, y = m(xarr, yarr)
				m.plot(x, y, linewidth=0.5, color=cl)

	fig.subplots_adjust(wspace=0)
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.05, dpi=300)

#################################################################################

def map_nwp_tracks_per_storm(nwp_infile_list, analysis_file,ibtracs_file, outfile, region,cat,name):
	"""Plots a map of storm tracks in the Southern Hemisphere
	Plots a map of storm tracks for one storm in the Southern Hemisphere
	from a list of several TRACK files (e.g. one per year) that have been reformatted to suit Python
	Plots all forecasts of one storm, and the reanalysis and ibtracs tracks.
	Plots the tracks so that they change colour from lighter to darker with time
	Note: the colormap used is hardcoded in the function (where 'colors' is assigned)
	Function takes a list of filenames as  input argument"""

	# set up map of region
	if region == "SH":
		lat1 = 30
		lat2 = -60
		lon1 = -25
		lon2 = 335

	elif region == "SIO":
		lat1=30
		lat2=-55
		lon1=-10
		lon2=130

	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')

	RSMClats = [-40, 0, 0, -40]
	RSMClons = [30, 30, 90, 90]

	draw_rectangle(RSMClats,RSMClons,m)

	#m.bluemarble(alpha=0.75)
	m.drawcoastlines(linewidth=0.4, color='lightgray')
	m.drawcountries(linewidth=0.4, color='lightgray')
	m.fillcontinents(color='white')


	no_tracks=len(nwp_infile_list)
	colors = cm.BuPu(np.linspace(0.2, 1, no_tracks)) #BuPu

	#loop over the files containing the forecast tracks, and over the colours to plot each track a different colour
	for file_in_list,c in zip(nwp_infile_list,colors):

		# load in data from a specified TRACK file (reformatted & interpolated)
		data = np.genfromtxt(file_in_list, dtype=float, skip_header=1)

		#plot this forecast track - loops over colours and plots each track a different colour
		x,y = m(data[:,7], data[:,8])
		m.plot(x,y,linewidth=0.5,color=c)

		#OR the following plots each track so it changes colour along its length

		# count how many track points this storm has
		#no_points = len(data[:,0])
		# set up an array to hold the lat and lon of each point of this storm
		#track_coords = np.zeros((no_points, 2))
		# get the lat and lons from the track file for this storm only
		#track_coords[:, 0] = data[:, 7]  # longitudes
		#track_coords[:, 1] = data[:, 8]  # latitudes

		# To have the line changing colour from start to finish, split the coordinates into segments
		#points = np.array([track_coords[:, 0], track_coords[:, 1]]).T.reshape(-1, 1, 2)
		#segments = np.concatenate([points[:-1], points[1:]], axis=1)

		#colors = cm.BuPu(np.linspace(0.2, 1, no_points))

		# iterate over each segment of this track
		#for p, cl in zip(range(no_points - 1), colors):
			# need to separate the x and y values for each segment, and plot each segment on the map
			#xarr = []
			#yarr = []
			#xarr.append(segments[p][0][0])
			#xarr.append(segments[p][1][0])
			#yarr.append(segments[p][0][1])
			#yarr.append(segments[p][1][1])
			#x, y = m(xarr, yarr)
			#m.plot(x, y, linewidth=0.5, color=cl)


	#plot the analysis and ibtracs tracks
	analysis_data = np.genfromtxt(analysis_file, dtype=float, skip_header=1)
	syear=int(analysis_data[0,3])
	smonth=int(analysis_data[0,4])
	sday=int(analysis_data[0,5])
	eyear=int(analysis_data[-1,3])
	emonth=int(analysis_data[-1,4])
	eday=int(analysis_data[-1,5])
	storm_dates=str(sday)+'/'+str(smonth)+'/'+str(syear)+' - '+str(eday)+'/'+str(emonth)+'/'+str(eyear)
	print storm_dates
	x,y = m(analysis_data[:,7], analysis_data[:,8])
	m.plot(x,y,linewidth=0.75, color='k')

	ibtracs_data = np.genfromtxt(ibtracs_file, dtype=float, skip_header=1)
	x,y = m(ibtracs_data[:,7], ibtracs_data[:,8])
	m.plot(x,y, linewidth=0.75, color='k',linestyle='--')

	#legend
	if cat == "TS":
		title = cat+" "+str(name)+"\n"+storm_dates
	else:
		title="Category "+str(cat)+" "+str(name)+"\n"+storm_dates
	an = plt.Line2D((0, 1), (0, 0), color='k',linewidth=0.5)
	ib = plt.Line2D((0, 1), (0, 0), color='k', linestyle='--',linewidth=0.5)
	nwp = plt.Line2D((0, 1), (0, 0), color='darkorchid',linewidth=0.5)
	legend = ax.legend((an, ib, nwp), ['Analysis', 'IBTrACS', 'UKMO Forecasts'],title=title, fontsize=5, loc='lower left')
	plt.setp(legend.get_title(), fontsize='5')
	legend._legend_box.align = "left"

	#save and close the plot
	fig.subplots_adjust(wspace=0)
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.05, dpi=500)
	plt.close()

#################################################################################

class MidpointNormalize(matplotlib.colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

	#def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		#x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		#return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def map_composite_data(data, lats, lons, outfile, region, cmap, label, bounds,dtype=None):
	"""Plots a map of data (e.g. track density, precip) in the Southern Hemisphere"""

	# set up map of region
	if region == "SH":
		lat1 = 30
		lat2 = -60
		lon1 = -25
		lon2 = 335

	elif region == "SIO":
		lat1 = 30
		lat2 = -55
		lon1 = -10
		lon2 = 130

	elif region == "SIO2":
		lat1 = 0
		lat2 = -30
		lon1 = 20
		lon2 = 60

	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])


	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')

	RSMClats = [-40, 0, 0, -40]
	RSMClons = [30, 30, 90, 90]

	# pl.draw_rectangle(RSMClats,RSMClons,m)

	m.drawcoastlines(linewidth=0.6, color='k') #gray
	m.drawcountries(linewidth=0.6, color='k') #gray

	cmap = plt.get_cmap(cmap)
	colors = cmap(np.linspace(0.1, 1, cmap.N))
	cmap2 = LinearSegmentedColormap.from_list(cmap, colors)

	max = round_up_to_even(np.max(data))
	print max
	#bounds = np.linspace(0, max, 11)
	#bounds = np.linspace(0,4000,9)
	print bounds
	#bounds = np.linspace(0, np.max(data))
	norm = BoundaryNorm(bounds, ncolors=cmap2.N, clip=True)

	# mask 0 values so that we can set them to white/transparent rather than pale blue
	data = np.ma.masked_where(data == 0, data)
	lons, lats = np.meshgrid(lons, lats)

	#im = m.pcolormesh(lons, lats, data, shading='flat', cmap=cmap, latlon=True, norm=norm)
	#im.cmap.set_bad('white', alpha=0)  # set missing values to white/transparent

	if dtype == 'bias':
		im = m.pcolormesh(lons, lats, data, shading='flat', cmap=cmap, latlon=True, norm = MidpointNormalize(midpoint=0., vmin=bounds[0], vmax=bounds[-1]) ) #vmin=bounds[0], vmax=bounds[-1])
		im.cmap.set_bad('white', alpha=0)  # set missing values to white/transparent
		cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="both")

	elif dtype == "perc":
		im = m.pcolormesh(lons, lats, data, shading='flat', cmap=cmap2, latlon=True, norm=norm)
		im.cmap.set_bad('white', alpha=0)  # set missing values to white/transparent
		cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, extend="max")

	elif dtype == None:
		im = m.pcolormesh(lons, lats, data, shading='flat', cmap=cmap2, latlon=True, norm=norm)
		im.cmap.set_bad('white', alpha=0)  # set missing values to white/transparent
		cbar = fig.colorbar(im, ax=ax,fraction=0.046, pad=0.04, extend="max")  # auto colorbar  #, extend="max"

	cbar.ax.tick_params(labelsize=8)  # colorbar font size
	plt.text(0.01, 0.02, label, transform=ax.transAxes, fontsize=6)
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.5, dpi=400)
	plt.close()

################################################################################# 
# TRMM DATA FUNCTIONS / CUBE FUNCTIONS 
#################################################################################   

def read_trmm_3b42(year, month, day, hour):
	"""Returns the precip data from a TRMM 3B42 file as a Cube (Python package Iris,
	which uses Cubes to hold data and metadata)."""

	# Convert arguments to integers
	year, month, day, hour = [int(x) for x in (year, month, day, hour)]

	# Callback function
	def cb(cube, field, filename):
		try:
			del cube.calendar
		except AttributeError:
			pass

	# Read data
	# iris.util.squeeze removes any dimensions with length 1 (e.g. a cube with shape (1, 360, 360) is changed to shape (360, 360)
	pcp = iris.util.squeeze(iris.load_cube(os.path.join(TRMM_DIR, str(year),
														TRMM_FILE % (year, month, day, hour)), 'precipitation (mm/hr)',
										   callback=cb))
	pcp.long_name = 'precipitation rate'
	pcp.standard_name = 'lwe_precipitation_rate'
	pcp.var_name = 'pcp'
	pcp.units = 'mm hr-1'

	# Return
	return pcp


#################################################################################

def add_cubes_sep_res(cubes, year, contributing_days=True, separate_resolutions=True,deal_with_masks=True):
	"""Returns the sum of the given Cubes.

	If deal_with_masks is True (default), treats missing data as 0 unless it is
	missing in the same place in all Cubes.  This behaviour must be turned off
	if the Cubes are of differing resolutions and separate_resolutions is False
	(which is presumably safe for model or reanalysis data, but possibly not
	for observations as they may contain missing values).

	If separate_resolutions is True, returns a dictionary with Cube shapes as
	the keys and Cubes as the values.
	"""

	# Set bounds and remove lat/lon coord attributes
	c0 = cubes[0]
	c1 = cubes[-1]
	for c in cubes[1:]:
		for ax in 'XY':
			if c.shape == c0.shape:
				if not c.coord(axis=ax) == c0.coord(axis=ax):
					c.replace_coord(c0.coord(axis=ax).copy())
			elif c.shape == c1.shape:
				if not c.coord(axis=ax) == c1.coord(axis=ax):
					c.replace_coord(c1.coord(axis=ax).copy())
	bounds(cubes)
	for c in cubes:
		for ax in 'XY':
			c.coord(axis=ax).attributes = {}

	# Deal with masked data
	if deal_with_masks:

		# If contributing_days is in all Cubes, add it up
		if contributing_days:
			count_days = sum([c.attributes['contributing_days'] for c in cubes])

		# Get masks
		masks = [c.data.mask if isinstance(c.data, np.ma.MaskedArray) else np.zeros_like(c.data).astype(bool) for c in cubes]

		# Set masked data to 0.0 and sum
		for c, m in itertools.izip(cubes, masks):
			c.data = np.where(m, 0.0, c.data)
		total = sum(cubes)

		# Mask data if missing in all Cubes
		mask = sum([arr.astype(int) for arr in masks]) == len(masks)
		total.data = np.ma.MaskedArray(total.data, mask=mask)
		total.metadata = cubes[0].metadata
		if contributing_days:
			total.attributes['contributing_days'] = float(count_days)
		return total

	# Ignore masks
	# Allow for different resolutions here
	# (What if we need to handle masked data at multiple resolutions?  We'll
	# worry about that if it ever actually happens...)
	else:

		# Sum over Cubes, doing each shape separately due to different resolutions
		totals_dict = {}

		if year == 2010:
			totals_dict['n320'] = {}
			totals_dict['n512'] = {}
			n320_count=0
			n512_count = 0

		elif year == 2014:
			totals_dict['n768'] = {}
			totals_dict['n512'] = {}
			n512_count = 0
			n768_count = 0


		for c in cubes:

			if year == 2010:

				if c.shape == (481, 640):
					if n320_count == 0:
						totals_dict['n320'] = c
					else:
						totals_dict['n320'] += c
					n320_count+=1

				elif c.shape == (769,1024):
					if n512_count == 0:
						pass
					elif n512_count == 1:
						totals_dict['n512'] = c
					else:
						totals_dict['n512'] += c
					n512_count += 1

			elif year == 2014:

				if c.shape == (769,1024):
					#print c
					if n512_count==0:
						#pass
						totals_dict['n512'] = c
					#elif n512_count==1:
						#totals_dict['n512'] = c
					else:
						totals_dict['n512'] += c
					n512_count+=1

				elif c.shape == (1152,1536):
					print c
					if n768_count==0:
						pass
						#totals_dict['n768'] = c
					elif n768_count==1:
						totals_dict['n768'] = c
					else:
						totals_dict['n768'] += c
					n768_count+=1


		if separate_resolutions:

			for c in totals_dict.itervalues():
				c.metadata = cubes[0].metadata

			return totals_dict



############################################################################################

def add_cubes(cubes, contributing_days=True, separate_resolutions=False,deal_with_masks=True):
	"""Returns the sum of the given Cubes.

	If deal_with_masks is True (default), treats missing data as 0 unless it is
	missing in the same place in all Cubes.  This behaviour must be turned off
	if the Cubes are of differing resolutions and separate_resolutions is False
	(which is presumably safe for model or reanalysis data, but possibly not
	for observations as they may contain missing values).

	If separate_resolutions is True, returns a dictionary with Cube shapes as
	the keys and Cubes as the values.
	"""

	# Set bounds and remove lat/lon coord attributes
	c0 = cubes[0]
	c1 = cubes[-1]
	#for coord in c0.coords():
		#print coord.long_name
		#print '##############'
	for c in cubes[1:]:
		#added this when had an error about longitude coordinate being missing in the regridded files
		#it's not missing, but the longname was missing, so take this from the one that has it!
		for coord,c0coord in zip(c.coords(), c0.coords()):
			print coord.long_name
			coord.long_name = c0coord.long_name
		for ax in 'XY':
			if c.shape == c0.shape:
				if not c.coord(axis=ax) == c0.coord(axis=ax):
					c.replace_coord(c0.coord(axis=ax).copy()) #original
			elif c.shape == c1.shape:
				if not c.coord(axis=ax) == c1.coord(axis=ax):
					c.replace_coord(c1.coord(axis=ax).copy())
	bounds(cubes)
	for c in cubes:
		for ax in 'XY':
			c.coord(axis=ax).attributes = {}

	# Deal with masked data
	if deal_with_masks:

		# If contributing_days is in all Cubes, add it up
		if contributing_days:
			count_days = sum([c.attributes['contributing_days'] for c in cubes])

		# Get masks
		masks = [c.data.mask if isinstance(c.data, np.ma.MaskedArray) else np.zeros_like(c.data).astype(bool) for c in cubes]

		# Set masked data to 0.0 and sum
		for c, m in itertools.izip(cubes, masks):
			c.data = np.where(m, 0.0, c.data)
		total = sum(cubes)

		# Mask data if missing in all Cubes
		mask = sum([arr.astype(int) for arr in masks]) == len(masks)
		total.data = np.ma.MaskedArray(total.data, mask=mask)
		total.metadata = cubes[0].metadata
		if contributing_days:
			total.attributes['contributing_days'] = float(count_days)
		return total

	# Ignore masks
	# Allow for different resolutions here
	# (What if we need to handle masked data at multiple resolutions?  We'll
	# worry about that if it ever actually happens...)
	else:

		# Sum over Cubes, doing each shape separately (in case we have
		# differing resolutions)
		totals_dict = {}

		#print "totals_dict.keys() : ", totals_dict.keys()
		if contributing_days:
			contrib_dict = {}
		for c in cubes:
			#print c.shape
			#print type(c)
			#print c
			if c.shape in totals_dict.keys():
				#try:
					#if c.shape == (769, 1024):
						#trycount+=1
				totals_dict[c.shape] += c
				#print totals_dict[c.shape]
				#print "totals_dict.keys() : ", totals_dict.keys()
				if contributing_days:
					contrib_dict[c.shape] += c.attributes['contributing_days']
				#except:
					#print c
					#if c.shape == (769, 1024):
						#errcount+=1

			else:
				#if c.shape == (769, 1024):
					#trycount += 1
				totals_dict[c.shape] = c
				#print "totals_dict.keys() : ", totals_dict.keys()
				if contributing_days:
					contrib_dict[c.shape] = c.attributes['contributing_days']

		#print "number of times it tried to add a file:", trycount
		#print "number of times files were passed over due to an error:", errcount
		# Return each resolution separately
		if separate_resolutions:
			#print totals_dict.iterkeys()
			#print totals_dict

			# Apply metadata
			for c in totals_dict.itervalues():
				c.metadata = cubes[0].metadata

			# Apply contributing_days attribute
			for k in totals_dict.iterkeys():
				if contributing_days:
					totals_dict[k].attributes['contributing_days'] = float(contrib_dict[k])
			return totals_dict


		# Return single Cube
		else:

			# Regrid data to lowest resolution if necessary and sum
			totals = iris.cube.CubeList(totals_dict.values())
			if len(totals) > 1:
				size = lambda c: np.array(c.shape).prod()
				sizes = np.array([size(c) for c in totals])
				smallest = np.argwhere(sizes == min(sizes)).min()
				for i, c in enumerate(totals):
					if i != smallest:
						totals[i] = c.regrid(totals[smallest], iris.analysis.AreaWeighted())
			total = sum(totals)

			# Apply metadata
			total.metadata = cubes[0].metadata

			# Apply contributing_days attribute
			if contributing_days:
				total.attributes['contributing_days'] = float(sum(contrib_dict.values()))
			return total

#################################################################################

def bounds(cubes, dims=None, **kwargs):
    """Applies contiguous bounds to given dimensions, or all dimensions if None
    given.

    **Arguments**

    *cubes*
        `iris.cube.Cube` or `iris.cube.CubeList`

    **Optional arguments**

    *dims*=None
        `str` or `list` of dimension names (uses all dimensions if None)

    **kwargs
        Will be passed to guess_bounds (see iris.coords.DimCoord.guess_bounds
        docstring for details)
    """

    # Make sure we have a CubeList
    if isinstance(cubes, iris.cube.Cube):
        cubes = iris.cube.CubeList([cubes])

    # Parse *dims*
    if isinstance(dims, str):
        dims = [dims]

    # Iterate for each Cube
    for cube in cubes:

        # List of dimensions
        if dims is None:
            dims_use = [c.name() for c in cube.dim_coords]
        else:
            dims_use = [d for d in dims]

        # Iterate for each dimension
        for dim in dims_use:
            coord = cube.coord(dim)

            # Co-ordinate length must be at least 2 to have bounds
            if len(coord.points) < 2:
                continue

            # If no bounds, set them
            if not coord.has_bounds():
                coord.guess_bounds(**kwargs)

            # If existing bounds are not contiguous, remake them
            elif not coord.is_contiguous():
                coord.bounds = None
                coord.guess_bounds(**kwargs)

            # If bounds are still not contiguous, adjust them manually
            if not coord.is_contiguous():
                b = coord._bounds.copy()
                for (ii, (b1, b2)) in enumerate(zip(b[:-1], b[1:])):
                    if b1[1] != b2[0]:
                        b[ii][1] = b[ii+1][0]
                coord._bounds = b


#################################################################################
# TRACK DENSITY
#################################################################################  

def track_density(infile, outfile, _source, year, month, box_size=2.5,
				  y1=None, y2=None, lead='all', exclude=[],
				  contributing_days=True):
	"""Calculates the TC track density as the number of tracks per grid box per
	year, with the given grid box size in degrees."""

	# Check whether output file already exists
	if os.path.isfile(outfile):
		raise ValueError('Output file %s already exists' % outfile)
	outdir = os.path.split(outfile)[0]
	if not os.path.isdir(outdir):
		os.makedirs(outdir)

	# Check arguments
	if (_source, lead) == ('ukmo_nwp', 'all'):
		msg = 'Must specify lead time for NWP'
		raise ValueError(msg)

	# Create array
	lats = np.arange(-90 + box_size / 2., 90., box_size)
	lons = np.arange(box_size / 2., 360., box_size)
	lat = iris.coords.DimCoord(lats, units='degrees_north', standard_name='latitude')
	lon = iris.coords.DimCoord(lons, units='degrees_east', standard_name='longitude')
	lat.guess_bounds()
	lon.guess_bounds()
	lat_bounds = lat.bounds
	lon_bounds = lon.bounds
	aout = np.zeros((len(lats), len(lons)), dtype=float)

	# Function for determining grid box position
	def box(bounds, value):
		match_box = ((bounds <= value) == [True, False])
		return int(np.argwhere(match_box.astype(int).sum(axis=1) == 2))

	# JRA55
	if _source == 'jra55':

		# Get track data, and select specified year and month
		ain = np.genfromtxt(infile, dtype=float, skip_header=1,
							usecols=range(3, 9))
		aint = ain[np.where((ain[:, 0] == year) & (ain[:, 1] == month))]

		# Count days
		t1 = datetime.datetime(year, month, 1)
		t2 = datetime.datetime(year + 1 if month == 12 else year,
							   1 if month == 12 else month + 1, 1)
		dt = datetime.timedelta(days=1)
		count_days = 0
		for tt in extrange(t1, t2, dt):
			if tt in exclude:
				continue
			if tt == datetime.datetime(2017, 7, 11):
				count_days += 0.5
			else:
				count_days += 1

		# Iterate for each point in track file
		end_time = datetime.datetime(2017, 7, 11, 12)
		for day, hour, lon_value, lat_value in aint[:, 2:]:
			tt = datetime.datetime(year, month, int(day))
			if tt in exclude:
				continue
			if tt + datetime.timedelta(hours=int(hour)) >= end_time:
				continue
			lon_index = box(lon_bounds, lon_value)
			lat_index = box(lat_bounds, lat_value)
			aout[lat_index, lon_index] += 1.

	# UKMO NWP
	elif _source == 'ukmo_nwp':

		# Iterate for each time available in this year and month
		# print year
		# print type(year)
		t1 = datetime.datetime(year, month, 1, 0)
		if month == 12:
			t2 = datetime.datetime(year + 1, 1, 1, 0)
		else:
			t2 = datetime.datetime(year, month + 1, 1, 0)
		dt = datetime.timedelta(days=1)
		count_days = 0
		for tt in extrange(t1, t2, dt):

			# Check whether to exclude this time
			if tt in exclude:
				print tt, '- EXCLUDE'
				continue
			if tt.timetuple()[:3] == (2017, 7, 11):
				count_days += 0.5
			else:
				count_days += 1

			# Get list of forecast and validity times for the three forecasts to be
			# used - three forecasts???
			ftime_deltas = np.arange(-12, 13, 12) - lead * 24
			ftimes = (tt + datetime.timedelta(hours=hh) for hh in ftime_deltas)
			vtimes = (np.array([15, 21]) + lead * 24, np.arange(3, 22, 6) + lead * 24,
					  np.array([3, 9]) + lead * 24)

			# Iterate for each of the three forecasts
			for ff, vv in itertools.izip(ftimes, vtimes):

				# Check track data exist for current forecast
				this_infile = ff.strftime(infile)
				if not os.path.isfile(this_infile):
					raise ValueError('%s - %s does not exist' % (tt, this_infile))

				# Get ID, year, month, day, hour, lon, lat from file
				with warnings.catch_warnings():
					warnings.simplefilter('ignore')
					ain = np.genfromtxt(this_infile, dtype=float, skip_header=1,
										usecols=[0] + range(3, 9))
				if not ain.size:
					continue

				# Iterate for each time
				for v in vv:
					vt = ff + datetime.timedelta(hours=v)
					aint = ain[np.where((ain[:, 1] == vt.year) & \
										(ain[:, 2] == vt.month) & \
										(ain[:, 3] == vt.day) & \
										(ain[:, 4] == vt.hour))]
					if not aint.size:
						continue

					# Add each track to aout
					for track_id, lon_value, lat_value in aint[:, [0, 5, 6]]:
						lon_index = box(lon_bounds, lon_value)
						lat_index = box(lat_bounds, lat_value)
						aout[lat_index, lon_index] += 1.

		# Divide by 2 (because there are 2 forecasts per day)
		aout /= 2.

	# Divide by number of years in file and save
	density = iris.cube.Cube(aout, var_name='density')
	density.add_dim_coord(lat, 0)
	density.add_dim_coord(lon, 1)
	if contributing_days:
		density.attributes['contributing_days'] = float(count_days)
	iris.save(density, outfile)

