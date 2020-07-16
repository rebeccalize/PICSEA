import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.cm as cm



datadir="/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"

#the following lists contain the number of cyclones in each intensity category, that made landfall in each of mozambique, madagascar and the seychelles


#July 2006 - June 2018:
# No Data, TD, MTS, STS, TC, ITC, VITC
data=[8,9,21,18,17,21]
data = np.array(data)
colours = ['silver','darkorange','black','orangered','firebrick','black']
hatches=[None, '/', '/',None, None, None]
total = 94

perc = np.zeros(6)

for i in range(6):
	p = (int(data[i])/94.)*100.
	print p
	perc[i] = p
	

print perc

fig, ax = plt.subplots()
fig.set_size_inches(5,4)

#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")

	

x = np.arange(6)

for i in range(6):
	
	b = plt.bar(x[i], data[i], color=colours[i], width=0.75, hatch=hatches[i], edgecolor='white')


#plt.xlabel('Cyclone Intensity', fontsize = 9)
plt.xticks(x, ('No Data','Moderate\nTropical Storm', 'Strong\nTropical Storm', 'Tropical\nCyclone', 'Intense\nTropical Cyclone', 'Very Intense\nTropical Cyclone'), fontsize=6, rotation=0)



plt.ylabel('Number of Cyclones (July 2010 - May 2020)', fontsize=9)
plt.ylim(0,22)
plt.yticks(np.arange(0,23,2), fontsize=6)

#ax.tick_params(axis='both', which='major', labelsize=6)

#plt.legend((b1, b2, b3), ('Madagascar', 'Mozambique', 'Seychelles'), fontsize=8)

plt.tight_layout()

plt.savefig("number_of_cyclones_per_category_2010-2020.png", dpi=400)

plt.close()



fig, ax = plt.subplots()
fig.set_size_inches(5,4)

#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")

	

x = np.arange(6)

for i in range(6):
	
	b = plt.bar(x[i], perc[i], color=colours[i], width=0.75, hatch=hatches[i], edgecolor='white')
	

#plt.xlabel('Cyclone Intensity', fontsize = 9)
plt.xticks(x, ('No Data','Moderate\nTropical Storm', 'Strong\nTropical Storm', 'Tropical\nCyclone', 'Intense\nTropical Cyclone', 'Very Intense\nTropical Cyclone'), fontsize=6, rotation=0)



plt.ylabel('Percentage of Cyclones (July 2010 - May 2020)', fontsize=9)
plt.ylim(0,30)
plt.yticks(np.arange(0,31,5), fontsize=6)

#ax.tick_params(axis='both', which='major', labelsize=6)

#plt.legend((b1, b2, b3), ('Madagascar', 'Mozambique', 'Seychelles'), fontsize=8)

plt.tight_layout()

plt.savefig("percentage_of_cyclones_per_category_2010-2020.png", dpi=400)

plt.close()











	



