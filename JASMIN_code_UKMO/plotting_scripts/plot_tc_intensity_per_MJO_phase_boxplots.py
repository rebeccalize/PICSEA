from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import os
import fnmatch
import numpy as np
import matplotlib.dates as mdates
import matplotlib.cbook as cbook
from scipy.interpolate import spline
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch
import matplotlib.patches as mpatches



fig, ax = plt.subplots()
fig.set_size_inches(16, 10)


box12data = np.genfromtxt('intensity_of_TCs_072006-062018_MJO_phase_12.txt', dtype=float)
box12winds = box12data[:,0]
box12winds = box12winds[np.logical_not(np.isnan(box12winds))]
box12mslp = box12data[:,1]
box12mslp = box12mslp[np.logical_not(np.isnan(box12mslp))]

box23data = np.genfromtxt('intensity_of_TCs_072006-062018_MJO_phase_23.txt', dtype=float)
box23winds = box23data[:,0]
box23winds = box23winds[np.logical_not(np.isnan(box23winds))]
box23mslp = box23data[:,1]
box23mslp = box23mslp[np.logical_not(np.isnan(box23mslp))]

box34data = np.genfromtxt('intensity_of_TCs_072006-062018_MJO_phase_34.txt', dtype=float)
box34winds = box34data[:,0]
box34winds = box34winds[np.logical_not(np.isnan(box34winds))]
box34mslp = box34data[:,1]
box34mslp = box34mslp[np.logical_not(np.isnan(box34mslp))]

box45data = np.genfromtxt('intensity_of_TCs_072006-062018_MJO_phase_45.txt', dtype=float)
box45winds = box45data[:,0]
box45winds = box45winds[np.logical_not(np.isnan(box45winds))]
box45mslp = box45data[:,1]
box45mslp = box45mslp[np.logical_not(np.isnan(box45mslp))]

box56data = np.genfromtxt('intensity_of_TCs_072006-062018_MJO_phase_56.txt', dtype=float)
box56winds = box56data[:,0]
box56winds = box56winds[np.logical_not(np.isnan(box56winds))]
box56mslp = box56data[:,1]
box56mslp = box56mslp[np.logical_not(np.isnan(box56mslp))]

box67data = np.genfromtxt('intensity_of_TCs_072006-062018_MJO_phase_67.txt', dtype=float)
box67winds = box67data[:,0]
box67winds = box67winds[np.logical_not(np.isnan(box67winds))]
box67mslp = box67data[:,1]
box67mslp = box67mslp[np.logical_not(np.isnan(box67mslp))]

box78data = np.genfromtxt('intensity_of_TCs_072006-062018_MJO_phase_78.txt', dtype=float)
box78winds = box78data[:,0]
box78winds = box78winds[np.logical_not(np.isnan(box78winds))]
box78mslp = box78data[:,1]
box78mslp = box78mslp[np.logical_not(np.isnan(box78mslp))]

box81data = np.genfromtxt('intensity_of_TCs_072006-062018_MJO_phase_81.txt', dtype=float)
box81winds = box81data[:,0]
box81winds = box81winds[np.logical_not(np.isnan(box81winds))]
box81mslp = box81data[:,1]
box81mslp = box81mslp[np.logical_not(np.isnan(box81mslp))]



data_std_pairs_wind=[box23winds,box45winds,box67winds,box81winds]
data_std_pairs_mslp=[box23mslp,box45mslp,box67mslp,box81mslp]
data_alt_pairs_wind=[box12winds,box34winds,box56winds,box78winds]
data_alt_pairs_mslp=[box12mslp,box34mslp,box56mslp,box78mslp]

print data_std_pairs_wind



#TO USE 5TH AND 95TH PERCENTILES, CHANGE whis=100 to whis=[5,95]
bp = plt.boxplot(data_std_pairs_wind, sym='+', vert=1, whis=100, widths=0.25, patch_artist=True)
plt.setp(bp['boxes'], color='black',linewidth=2)
plt.setp(bp['whiskers'], color='black',linewidth=2,linestyle='solid')
#plt.setp(bp['fliers'], color='red', marker='o',alpha=0.5)
#plt.setp(bp['medians'], color='black',linewidth=1)
plt.setp(bp['caps'],color='black',linewidth=2)

#box = plt.boxplot(data, patch_artist=True)
#colors = ['#1A5899','#1A5899','#9C1127','#529DC8','#BBDAEA','#529DC8','#529DC8','white','#DB6B55','#DB6B55','#529DC8','#9C1127','#9C1127','#DB6B55' ,]
colors = ['#fcc200', '#f05238', '#a1005c', '#08025c'] 
for patch, color,flier,median in zip(bp['boxes'], colors,bp['fliers'],bp['medians']):
    patch.set(facecolor=color,alpha=0.9)
    flier.set(marker='o', color=color, alpha=0.9,markersize=3)
    median.set(color='k',linewidth=2)


#for label positions
xpos=[1,2,3,4]

my_xticks = ['MJO 2-3', 'MJO 4-5','MJO 6-7','MJO 8-1']
#my_yticks = ['0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60',' ','70',' ','80',' ','90',' ','100']
plt.yticks(fontsize=24)
plt.xticks(xpos, my_xticks,fontsize=24)
plt.xlabel('MJO Phase (Standard Pairings)', fontsize = 26)
plt.ylabel(r'Max TC Wind Speed (ms$^{-1}$)', fontsize = 26)


plt.savefig('TC_wind_speed_ibtracs_per_MJO_phase_standard_pairings.png', bbox_inches='tight', dpi = 400)
plt.close()


################################################################################################################

fig, ax = plt.subplots()
fig.set_size_inches(16, 10)

#TO USE 5TH AND 95TH PERCENTILES, CHANGE whis=100 to whis=[5,95]
bp = plt.boxplot(data_alt_pairs_wind, sym='+', vert=1, whis=100, widths=0.25, patch_artist=True)
plt.setp(bp['boxes'], color='black',linewidth=2)
plt.setp(bp['whiskers'], color='black',linewidth=2,linestyle='solid')
#plt.setp(bp['fliers'], color='red', marker='o',alpha=0.5)
#plt.setp(bp['medians'], color='black',linewidth=1)
plt.setp(bp['caps'],color='black',linewidth=2)

#box = plt.boxplot(data, patch_artist=True)
#colors = ['#1A5899','#1A5899','#9C1127','#529DC8','#BBDAEA','#529DC8','#529DC8','white','#DB6B55','#DB6B55','#529DC8','#9C1127','#9C1127','#DB6B55' ,]
colors = ['#fcc200', '#f05238', '#a1005c', '#08025c'] 
for patch, color,flier,median in zip(bp['boxes'], colors,bp['fliers'],bp['medians']):
    patch.set(facecolor=color,alpha=0.9)
    flier.set(marker='o', color=color, alpha=0.9,markersize=3)
    median.set(color='k',linewidth=2)


#for label positions
xpos=[1,2,3,4]

my_xticks = ['MJO 1-2', 'MJO 3-4','MJO 5-6','MJO 7-8']
#my_yticks = ['0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60',' ','70',' ','80',' ','90',' ','100']
plt.yticks(fontsize=24)
plt.xticks(xpos, my_xticks,fontsize=24)
plt.xlabel('MJO Phase (Alternative Pairings)', fontsize = 26)
plt.ylabel(r'Max TC Wind Speed (ms$^{-1}$)', fontsize = 26)


plt.savefig('TC_wind_speed_ibtracs_per_MJO_phase_alternative_pairings.png', bbox_inches='tight', dpi = 400)
plt.close()



###################################################################################################################

fig, ax = plt.subplots()
fig.set_size_inches(16, 10)

#TO USE 5TH AND 95TH PERCENTILES, CHANGE whis=100 to whis=[5,95]
bp = plt.boxplot(data_std_pairs_mslp, sym='+', vert=1, whis=100, widths=0.25, patch_artist=True)
plt.setp(bp['boxes'], color='black',linewidth=2)
plt.setp(bp['whiskers'], color='black',linewidth=2,linestyle='solid')
#plt.setp(bp['fliers'], color='red', marker='o',alpha=0.5)
#plt.setp(bp['medians'], color='black',linewidth=1)
plt.setp(bp['caps'],color='black',linewidth=2)

#box = plt.boxplot(data, patch_artist=True)
#colors = ['#1A5899','#1A5899','#9C1127','#529DC8','#BBDAEA','#529DC8','#529DC8','white','#DB6B55','#DB6B55','#529DC8','#9C1127','#9C1127','#DB6B55' ,]
colors = ['#fcc200', '#f05238', '#a1005c', '#08025c'] 
for patch, color,flier,median in zip(bp['boxes'], colors,bp['fliers'],bp['medians']):
    patch.set(facecolor=color,alpha=0.9)
    flier.set(marker='o', color=color, alpha=0.9,markersize=3)
    median.set(color='k',linewidth=2)


#for label positions
xpos=[1,2,3,4]

my_xticks = ['MJO 2-3', 'MJO 4-5','MJO 6-7','MJO 8-1']
#my_yticks = ['0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60',' ','70',' ','80',' ','90',' ','100']
plt.yticks(fontsize=24)
plt.xticks(xpos, my_xticks,fontsize=24)
plt.xlabel('MJO Phase (Standard Pairings)', fontsize = 26)
plt.ylabel('Minimum Pressure (hPa)', fontsize = 26)


plt.savefig('TC_mslp_ibtracs_per_MJO_phase_standard_pairings.png', bbox_inches='tight', dpi = 400)
plt.close()


################################################################################################################

fig, ax = plt.subplots()
fig.set_size_inches(16, 10)

#TO USE 5TH AND 95TH PERCENTILES, CHANGE whis=100 to whis=[5,95]
bp = plt.boxplot(data_alt_pairs_mslp, sym='+', vert=1, whis=100, widths=0.25, patch_artist=True)
plt.setp(bp['boxes'], color='black',linewidth=2)
plt.setp(bp['whiskers'], color='black',linewidth=2,linestyle='solid')
#plt.setp(bp['fliers'], color='red', marker='o',alpha=0.5)
#plt.setp(bp['medians'], color='black',linewidth=1)
plt.setp(bp['caps'],color='black',linewidth=2)

#box = plt.boxplot(data, patch_artist=True)
#colors = ['#1A5899','#1A5899','#9C1127','#529DC8','#BBDAEA','#529DC8','#529DC8','white','#DB6B55','#DB6B55','#529DC8','#9C1127','#9C1127','#DB6B55' ,]
colors = ['#fcc200', '#f05238', '#a1005c', '#08025c'] 
for patch, color,flier,median in zip(bp['boxes'], colors,bp['fliers'],bp['medians']):
    patch.set(facecolor=color,alpha=0.9)
    flier.set(marker='o', color=color, alpha=0.9,markersize=3)
    median.set(color='k',linewidth=2)


#for label positions
xpos=[1,2,3,4]

my_xticks = ['MJO 1-2', 'MJO 3-4','MJO 5-6','MJO 7-8']
#my_yticks = ['0',' ','10',' ','20',' ','30',' ','40',' ','50',' ','60',' ','70',' ','80',' ','90',' ','100']
plt.yticks(fontsize=24)
plt.xticks(xpos, my_xticks,fontsize=24)
plt.xlabel('MJO Phase (Alternative Pairings)', fontsize = 26)
plt.ylabel('Minimum Pressure (hPa)', fontsize = 26)


plt.savefig('TC_mslp_ibtracs_per_MJO_phase_alternative_pairings.png', bbox_inches='tight', dpi = 400)
plt.close()
