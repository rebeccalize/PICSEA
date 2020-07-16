import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.cm as cm

def bracket(ax, pos=[0,0], scalex=0.8, scaley=1, text="",textkw = {}, linekw = {}):
    x = np.array([0, 0.05, 0.45,0.5])
    y = np.array([0,-0.01,-0.01,-0.02])
    x = np.concatenate((x,x+0.5)) 
    y = np.concatenate((y,y[::-1]))
    ax.plot(x*scalex+pos[0], y*scaley+pos[1], clip_on=False, 
            transform=ax.get_xaxis_transform(), **linekw)
    ax.text(pos[0]+0.5*scalex, (y.min()-0.01)*scaley+pos[1], text, 
                transform=ax.get_xaxis_transform(),
                ha="center", va="top", **textkw)


#total number of TCs per MJO phase, based on first point of the analysis track
p12=26
p23=18
p34=16
p45=15
p56=26
p67=30
p78=19
p81=25

#total number of TCs per MJO phase per year (TC season from 2006-2007 to 2017-2018):

#old totals using first point of the ibtracs track
#p12y=[2,5,1,3,0,4,2,1,1,0,3,2]
#p23y=[2,4,1,3,0,4,0,2,0,1,3,4]
#p34y=[1,2,1,2,0,0,0,3,0,5,0,3]
#p45y=[1,1,1,1,0,1,0,1,2,4,1,1]
#p56y=[2,0,3,2,2,3,1,0,2,0,1,0]
#p67y=[1,2,3,3,2,2,2,2,5,0,0,0]
#p78y=[0,2,3,3,0,0,2,2,7,3,0,3]
#p81y=[1,2,3,3,0,0,3,1,3,3,0,3]

#new totals using first point of the analysis track
p12y=[1,8,0,2,0,4,4,2,0,0,2,3]
p23y=[0,3,0,2,0,2,1,1,1,2,2,3]
p34y=[2,0,0,3,1,1,0,1,2,4,1,1]
p45y=[2,0,2,3,1,3,1,0,1,2,0,0]
p56y=[3,0,5,4,1,4,2,4,3,0,0,0]
p67y=[3,3,4,4,1,2,2,4,5,0,0,2]
p78y=[0,3,3,2,1,0,1,0,3,2,0,4]
p81y=[1,5,2,2,1,2,3,2,1,2,1,3]



fig, ax = plt.subplots()
fig.set_size_inches(5,4)

#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")

	
#data = [p12,p34,p56,p78]
x = np.arange(8)

colors=['#fcc200', '#f05238', '#a1005c', '#08025c']

b1 = plt.bar(x[0], p12, color=colors[0], width=0.25)
b2 = plt.bar(x[1], p23, color=colors[0], width=0.25,hatch="/", edgecolor='white')
b3 = plt.bar(x[2], p34, color=colors[1], width=0.25)
b4 = plt.bar(x[3], p45, color=colors[1], width=0.25,hatch="/", edgecolor='white')
b5 = plt.bar(x[4], p56, color=colors[2], width=0.25)
b6 = plt.bar(x[5], p67, color=colors[2], width=0.25,hatch="/", edgecolor='white')
b7 = plt.bar(x[6], p78, color=colors[3], width=0.25)
b8 = plt.bar(x[7], p81, color=colors[3], width=0.25,hatch="/", edgecolor='white')

#plt.xlabel('MJO Phase (Alternative Pairings)', fontsize = 9)
plt.xticks(x, ('MJO 1-2', 'MJO 2-3', 'MJO 3-4', 'MJO 4-5','MJO 5-6', 'MJO 6-7','MJO 7-8','MJO 8-1'), fontsize=8)

plt.ylabel('Number of Cyclones', fontsize=9)
plt.ylim(0,30)
plt.yticks(np.arange(0,31,2), fontsize=8)

#ax.tick_params(axis='both', which='major', labelsize=6)

#plt.legend((b1, b2, b3), ('Madagascar', 'Mozambique', 'Seychelles'), fontsize=8)

plt.tight_layout()

plt.savefig("number_of_cyclones_per_MJO_phase_both_pairings_analysisstart.png", bbox_inches='tight', dpi=400)

plt.close()


##########################################################################################################################


fig, ax = plt.subplots()
fig.set_size_inches(5,4)

data = [p12y,p34y,p56y,p78y]
x = np.arange(12)

b1 = plt.bar(x-0.4, data[0], color=colors[0], width=0.2,align='center')
b2 = plt.bar(x-0.2, data[1], color=colors[1], width=0.2,align='center')
b3 = plt.bar(x, data[2], color=colors[2], width=0.2,align='center')
b4 = plt.bar(x+0.2, data[3], color=colors[3], width=0.2,align='center')

#plt.xlabel('MJO Phase (Alternative Pairings)', fontsize = 9)
plt.xticks(x-0.1, ('2006-2007', '2007-2008', '2008-2009','2009-2010', '2010-2011','2011-2012','2012-2013','2013-2014','2014-2015','2015-2016','2016-2017','2017-2018'), fontsize=6,rotation=40)

bracket(ax, text='',pos=[-0.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[0.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[1.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[2.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[3.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[4.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[5.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[6.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[7.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[8.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[9.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[10.5,-0.01],linekw=dict(color='k',lw=0.75))



plt.ylabel('Number of Cyclones', fontsize=9)
plt.ylim(0,8)
plt.yticks(np.arange(0,9,2), fontsize=8)

#ax.tick_params(axis='both', which='major', labelsize=6)

plt.legend((b1, b2, b3, b4), ('MJO 1-2', 'MJO 3-4', 'MJO 5-6', 'MJO 7-8'), fontsize=8)

plt.tight_layout()

plt.savefig("number_of_cyclones_per_year_per_MJO_phase_alternative_pairings_analysisstart.png", bbox_inches='tight', dpi=400)

plt.close()




##########################################################################################################################


fig, ax = plt.subplots()
fig.set_size_inches(5,4)

data = [p23y,p45y,p67y,p81y]
x = np.arange(12)

b1 = plt.bar(x-0.4, data[0], color=colors[0], width=0.2,align='center')
b2 = plt.bar(x-0.2, data[1], color=colors[1], width=0.2,align='center')
b3 = plt.bar(x, data[2], color=colors[2], width=0.2,align='center')
b4 = plt.bar(x+0.2, data[3], color=colors[3], width=0.2,align='center')

#plt.xlabel('MJO Phase (Standard Pairings)', fontsize = 9)
plt.xticks(x-0.1, ('2006-2007', '2007-2008', '2008-2009','2009-2010', '2010-2011','2011-2012','2012-2013','2013-2014','2014-2015','2015-2016','2016-2017','2017-2018'), fontsize=6,rotation=40)

bracket(ax, text='',pos=[-0.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[0.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[1.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[2.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[3.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[4.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[5.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[6.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[7.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[8.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[9.5,-0.01],linekw=dict(color='k',lw=0.75))
bracket(ax, text='',pos=[10.5,-0.01],linekw=dict(color='k',lw=0.75))



plt.ylabel('Number of Cyclones', fontsize=9)
plt.ylim(0,8)
plt.yticks(np.arange(0,9,2), fontsize=8)

#ax.tick_params(axis='both', which='major', labelsize=6)

plt.legend((b1, b2, b3, b4), ('MJO 2-3', 'MJO 4-5', 'MJO 6-7', 'MJO 8-1'), fontsize=8)

plt.tight_layout()

plt.savefig("number_of_cyclones_per_year_per_MJO_phase_standard_pairings_analysisstart.png", bbox_inches='tight', dpi=400)

plt.close()
