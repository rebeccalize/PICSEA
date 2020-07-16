import sys
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure as fig
import matplotlib.cm as cm
from netCDF4 import Dataset
import os
import fnmatch
from math import *
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime, timedelta
from matplotlib.ticker import MaxNLocator
#from mpl_toolkits.basemap import Basemap




# Creates a dictionary of empty arrays which are filled with statistics later
# Each array can be called separately, but everything is stored within varstruct
def initialise_variable_struct(NT):
	varstruct = {}
	varstruct['det'] = {}
	varstruct['det']['bias'] = {}
	varstruct['det']['bias']['sum'] = np.zeros(NT)
	varstruct['det']['bias']['wgt'] = np.zeros(NT)
	varstruct['det']['bias']['avg'] = np.zeros(NT)
	varstruct['det']['trerr'] = {}
	varstruct['det']['trerr']['sum'] = np.zeros(NT)
	varstruct['det']['trerr']['wgt'] = np.zeros(NT)
	varstruct['det']['trerr']['avg'] = np.zeros(NT)
	return varstruct
	
	
def init_speed_struct(NT, nens): #nens is number of ensemble members - varies depending on the year for UKMO ensemble
	speedstruct = {}
	speedstruct['analysis'] = np.zeros(NT)
	speedstruct['ibtracs'] = np.zeros(NT)
	speedstruct['mean'] = np.zeros(NT)
	speedstruct['det'] = np.zeros(NT)
	speedstruct['eps'] = {}	
	for i in range(5,nens+5):
		speedstruct['eps']['i'] = np.zeros(NT)
	return speedstruct

#operator functions (basic calcs of bias and track error)
def bias(obs_data,fcst_data):
	#first, catch whether the forecast or observed data have missing values (1e+25)
	if fcst_data > 1000000:
		bias = np.ma.masked
	elif obs_data > 1000000:
		bias = np.ma.masked
	#then, if the value isn't masked, go ahead and work out the bias!
	elif not fcst_data == np.ma.masked:
		bias = fcst_data-obs_data
	else:
		bias = np.ma.masked
	return bias


def trerr(obs_data,fcst_data):
	"""Returns distance in km (great circle distance) between 2 points"""

	#print "obs data for trerr:"
	#print obs_data
	#make sure data read in is in the format [lon,lat]

	#trerr=np.ma.masked_all(len(obs_data))
	R = 6378.1  # Earth's equatorial radius

	#calculate great circle distance between the 2 (haversine formula)
	#if necessary, could instead use Vincenty formula (see bookmarked blog post with code for this),
	#which treats earth as ellipsoid rather than sphere, and is more accurate over long distances

	#print "### Calculating track error ###"
	#print "Obs coordinates, lon-lat:", obs_data
	#print "Fcst coordinates, lon-lat:", fcst_data

	if not fcst_data[0] == np.ma.masked:
		# convert coordinates to radians (required for math functions)
		obs_lat=np.deg2rad(obs_data[1])
		obs_lon=np.deg2rad(obs_data[0])
		fcst_lat=np.deg2rad(fcst_data[1])
		fcst_lon=np.deg2rad(fcst_data[0])

		dlon = fcst_lat - obs_lat
		dlat = fcst_lon - obs_lon

		a = sin(dlat / 2) ** 2 + cos(obs_lat) * cos(fcst_lat) * sin(dlon / 2) ** 2
		c = 2 * atan2(sqrt(a), sqrt(1 - a))

		err = R * c
		#print err
	else:
		err = np.ma.masked

	return err
	
	
	
# Compute the propagation speed of storms in forecasts and analysis / ibtracs
# Using Haversine Formula
#This bit just does the speed from one point to the next; "prop_speed_vals" function computes along length of track
def compute_prop_speed(lon1, lat1, lon2, lat2):

    deg2rad = pi / 180.0
    kmh2ms = 1.0/3.6
    earthrad = 6371

    # Convert decimal degrees to radians
    # Enforce that ALL inputs are masked_arrays
    lon1 = np.ma.masked_array(lon1) * deg2rad
    lat1 = np.ma.masked_array(lat1) * deg2rad
    lon2 = np.ma.masked_array(lon2) * deg2rad
    lat2 = np.ma.masked_array(lat2) * deg2rad

    # Haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    km = earthrad * c 
    
    # Calculate speed (6 hour time-steps)
    # Convert to useful units of m/s
    speedkmh = km/6              
    speedms = speedkmh * kmh2ms
    return speedms
    
    
# Compute prop speed values throughout length of track
# I think lon and lat should be arrays; anf fcst should be the type of track we're computing the speed for???
def prop_speed_vals(lons, lats): #, speed_struct, fcst
	speed_array = []
	for i in range(len(lons)-1):  # loop until len(time)-1 as we're taking a forward difference
		speed_array.append(compute_prop_speed(lons[i], lats[i], lons[i+1], lats[i+1]))
	return speed_array
	

def prop_speed_bias(obsdata, fcstdata):
	speed_bias = fcstdata-obsdata #forecast value - analysis or ibtracs value
	return speed_bias

def prop_speed_abs_err(obsdata, fcstdata):
	speed_err = abs(fcstdata-obsdata) #forecast value - analysis or ibtracs value
	return speed_err
	
	

	
	
	
	

# Computes the spread in propagation speed
# Method is slightly different to other statistics due to different data structure
# Fills the sum and wgt arrays for the prop_speed spread
def compute_speed_spread(speed_struct, stat_struct_spread):
	diffsum = np.zeros(NT)
	for i in range(4,24):
		diff=abs(speed_struct['mean'][:] - speed_struct['eps']['i'][:])
		for j in range(NT):
			if not (diff[j] is np.ma.masked):
				diffsum[j] += diff[j] 	
	spread = np.ma.fix_invalid(diffsum/20)
	for k in range(NT):
		if not (spread[k] is np.ma.masked):
			stat_struct_spread['sum'][k] += spread[k]
			stat_struct_spread['wgt'][k] += 1
	return stat_struct_spread
    
    
    
	







# Computes error on variable (using operator function) for specified forecast
# Appends these errors to a 1-D masked array
# If this array contains a value, adds it to array 'sum' in dictionary - sums errors across storm
# Adds 1 to the weight for each file - later divide by file number to compute average
# Note: the last parameter, error_struct, is a dictionary with a ['wgt'] and a ['sum'] member
def sum_errors(obs_data, fcst_data, operator, error_struct):
	#print "obs_data in sum_errors"
	#print obs_data
	var_err = np.ma.array(operator(obs_data, fcst_data))
	#print var_err#
	for i in range(len(var_err)):
		#print var_err[i]
		if var_err[i] == np.ma.masked:
			continue
		elif var_err[i] == np.nan:
			continue
		else:
			error_struct['sum'][i] += var_err[i]
			error_struct['wgt'][i] += 1
	return error_struct

# Instructs sum_errors to compute statistics for each forecast and error type combination
# Fills the empty arrays in the relevant dictionary
def compute_storm_stats(obs_data, fcst_data, stat, stat_struct):
	if stat == 'trerr':
		sum_errors(obs_data, fcst_data, trerr, stat_struct['det']['trerr'])
	elif stat == 'bias':
		sum_errors(obs_data, fcst_data, bias, stat_struct['det']['bias'])
	return stat_struct


# Computes average error per storm for all statistic combinations
# Fills the empty 'avg' arrays in the relevant dictionary
# Returns an array which has masked values rather than NaN as NaN cannot be processed later
def compute_storm_average(err_struct):
	err_struct['det']['bias']['avg'] = np.ma.fix_invalid(err_struct['det']['bias']['sum'] / err_struct['det']['bias']['wgt'])
	err_struct['det']['err']['avg'] = np.ma.fix_invalid(err_struct['det']['err']['sum'] / err_struct['det']['err']['wgt'])
	return err_struct


# -------------------------------------------------------------------------------------

# HELPER FUNCTIONS - COMPUTE STATISTICS FOR ENTIRE DATASET (AVERAGED OVER ALL STORMS)

# -------------------------------------------------------------------------------------


# Creates a dictionary filled with empty arrays
# These are later filled with statistics averaged over all storms
def init_all_storms_struct(NT):
	avgstruct = {}
	avgstruct['det'] = {}
	avgstruct['det']['bias'] = {}
	avgstruct['det']['bias']['sum'] = np.zeros(NT)
	avgstruct['det']['bias']['wgt'] = np.zeros(NT)
	avgstruct['det']['bias']['avg'] = np.zeros(NT)
	avgstruct['det']['trerr'] = {}
	avgstruct['det']['trerr']['sum'] = np.zeros(NT)
	avgstruct['det']['trerr']['wgt'] = np.zeros(NT)
	avgstruct['det']['trerr']['avg'] = np.zeros(NT)

	return avgstruct


# Fills the dictionary 'sum' arrays with the sum of statistics over all storms
# Fills the dictionary 'wgt' arrays with the number of files contributing to 'sum'
def sum_all_storms_stats(error_struct, storm_avg_struct):
	for i in range(NT):
		if not(error_struct['avg'][i] is np.ma.masked):
			storm_avg_struct['sum'][i] += error_struct['avg'][i]
			storm_avg_struct['wgt'][i] += 1
	return storm_avg_struct



# Applies the sum_all_storms_stats function to all combinations of forecast type and error type
def compute_all_storms_stats(error_struct, avg_struct):
	sum_all_storms_stats(error_struct['det']['err'], avg_struct['det']['err'])
	sum_all_storms_stats(error_struct['det']['bias'], avg_struct['det']['bias'])
	return avg_struct, error_struct

# Averages statistics across all storms (entire dataset) using 'sum' and 'wgt' arrays
# Fills the empty dictionary 'avg' arrays with these results
# Converts any NaN values in the array to masked
def compute_all_storms_average(avg_struct):
	avg_struct['det']['err']['avg'] = np.ma.fix_invalid(avg_struct['det']['err']['sum'] / avg_struct['det']['err']['wgt'])
	avg_struct['det']['bias']['avg'] = np.ma.fix_invalid(avg_struct['det']['bias']['sum'] / avg_struct['det']['bias']['wgt'])

















