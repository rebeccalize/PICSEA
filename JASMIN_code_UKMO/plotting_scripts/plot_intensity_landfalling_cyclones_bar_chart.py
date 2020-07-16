import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.cm as cm



datadir="/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_track_files_by_storm/"

#the following lists contain the number of cyclones in each intensity category, that made landfall in each of mozambique, madagascar and the seychelles


#July 2006 - June 2018:
# TD, MTS, STS, TC, ITC, VITC
moz=[0, 4, 1, 1, 4] #TS, 3, 1, TS, 4, TS, 3, TS, 4, TS, TS, 2, TD, TS, 3, 4
mada= [0, 4, 1, 7, 1] #4, 1, 4, 4, 3, 1, TS, 4, 3, TS, TS, TS, 4, 1, TS, TS, 4, TS, TS, 3, TS, 4, TS, 4, 4, 3, 4, 4, TS, 4, TS, TS, TS, TS, 5, TS, 4, 2, TD, 3, TS, 1, 3, TS, 4, 3 
sey = [0, 1, 0, 2, 2] #4, TS, TS, 5, 4


#And as a percentage of the total number of storms...

moz_p =[None] * 6 #16
mada_p = [None] * 6 #46
sey_p = [None] * 6 #5


#calculate percentage of total storms in each category
for i in range(5):
	moz_p[i] = (moz[i] / 16.) * 100.
	mada_p[i] = (mada[i] / 46.) * 100.
	sey_p[i] = (sey[i] / 5.) * 100.
	
print moz_p
print mada_p
print sey_p




fig, ax = plt.subplots()
fig.set_size_inches(5,4)

#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")

	
data = [mada, moz, sey]
x = np.arange(5)

b1 = plt.bar(x+0.00, data[0], color='#bc5090', width=0.25)
b2 = plt.bar(x+0.25, data[1], color='#003f5c', width=0.25)
b3 = plt.bar(x+0.50, data[2], color='#ffa600', width=0.25)

#plt.xlabel('Cyclone Intensity', fontsize = 9)
plt.xticks(x+0.25, ('Moderate\nTropical Storm', 'Strong\nTropical Storm', 'Tropical\nCyclone', 'Intense\nTropical Cyclone', 'Very Intense\nTropical Cyclone'), fontsize=6.5, rotation=0)

plt.ylabel('Number of Cyclones (July 2010 - May 2020)', fontsize=9)
plt.ylim(0,10)
plt.yticks(np.arange(0,11,2), fontsize=6.5)

#ax.tick_params(axis='both', which='major', labelsize=6)

plt.legend((b1, b2, b3), ('Madagascar (13)', 'Mozambique (10)', 'Seychelles (5)'), fontsize=8)

plt.tight_layout()

plt.savefig("number_of_cyclones_per_category_each_country_2010-2020.png", dpi=400)

plt.close()



fig, ax = plt.subplots()
fig.set_size_inches(5,4)

#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")


data_p = [mada_p, moz_p, sey_p]
x = np.arange(5)

b1 = plt.bar(x+0.00, data_p[0], color='#bc5090', width=0.25)
b2 = plt.bar(x+0.25, data_p[1], color='#003f5c', width=0.25)
b3 = plt.bar(x+0.50, data_p[2], color='#ffa600', width=0.25)

plt.xlabel('Cyclone Intensity\n(equivalent SSHS category)', fontsize = 9)
plt.xticks(x+0.25, ('Moderate\nTropical Storm', 'Strong\nTropical Storm', 'Tropical\nCyclone', 'Intense\nTropical Cyclone', 'Very Intense\nTropical Cyclone'), fontsize=8,rotation=45)

plt.ylabel('% of Cyclones', fontsize=9)
plt.ylim(0,100)
plt.yticks(np.arange(0,101,10), fontsize=8)

#ax.tick_params(axis='both', which='major', labelsize=6)

plt.legend((b1, b2, b3), ('Madagascar (46)', 'Mozambique (16)', 'Seychelles (5)'), fontsize=8)

plt.tight_layout()

#plt.savefig("percentage_of_cyclones_per_category_each_country_2006-2018.png", dpi=400)

plt.close()



fig, ax = plt.subplots()
fig.set_size_inches(4.5,6)

#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")


TD=np.array([6.25, 2.1739130434782608, 0.0])
TS=np.array([43.75, 39.130434782608695,40.0])
C1=np.array([6.25, 8.695652173913043, 0.0])
C2=np.array([6.25, 2.1739130434782608, 0.0])
C3=np.array([18.75, 15.217391304347828, 0.0])
C4=np.array([18.75, 30.434782608695656, 40.0])
C5=np.array([0.0, 2.1739130434782608, 20.0])

print TD[0]+TS[0]+C1[0]+C2[0]+C3[0]+C4[0]+C5[0]
print TD[1]+TS[1]+C1[1]+C2[1]+C3[1]+C4[1]+C5[1]
print TD[2]+TS[2]+C1[2]+C2[2]+C3[2]+C4[2]+C5[2]

x=np.arange(3)
w=0.4

c = cm.YlGnBu(np.linspace(0.1, 1, 7))

b1 = plt.bar(x, TD,color=c[0], width=w)
b2 = plt.bar(x, TS,color=c[1], width=w, bottom=TD)
b3 = plt.bar(x, C1,color=c[2], width=w, bottom=TD+TS)
b4 = plt.bar(x, C2,color=c[3], width=w, bottom=TD+TS+C1)
b5 = plt.bar(x, C3,color=c[4], width=w, bottom=TD+TS+C1+C2)
b6 = plt.bar(x, C4,color=c[5], width=w, bottom=TD+TS+C1+C2+C3)
b7 = plt.bar(x, C5,color=c[6], width=w, bottom=TD+TS+C1+C2++C3+C4)


#plt.xlabel('Country', fontsize = 9)
plt.xticks(x, ('Madagascar\n(46)', 'Mozambique\n(16)', 'Seychelles\n(5)'), fontsize=11)

plt.ylabel('% of Cyclones\n(equivalent SSHS category)', fontsize=11)
plt.ylim(0,100)
plt.yticks(np.arange(0,101,10), fontsize=10)


#plt.legend((b1[0], b2[0], b3[0], b4[0], b5[0], b6[0], b7[0]), ('TD', 'TS', '1', '2', '3', '4', '5'), bbox_to_anchor=(1.01, 1.0))
plt.legend((b7[0],b6[0], b5[0], b4[0], b3[0], b2[0], b1[0]), ('5','4','3','2','1','TS','TD'), bbox_to_anchor=(1.01, 1.0))

plt.tight_layout()

#plt.savefig("stacked_barchart_percentage_of_cyclones_per_category_each_country_2006-2018.png", dpi=400)

plt.close()









	



