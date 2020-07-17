import fnmatch
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm
import numpy as np
import matplotlib.path as mpath



datadir = "/perm/mo/more/picsea/IDAI_KENNETH/reformatted_idai_kenneth_tracks/"


def get_hurricane_symbol():
    u = np.array([  [2.444,7.553],
                    [0.513,7.046],
                    [-1.243,5.433],
                    [-2.353,2.975],
                    [-2.578,0.092],
                    [-2.075,-1.795],
                    [-0.336,-2.870],
                    [2.609,-2.016]  ])
    u[:,0] -= 0.098
    codes = [1] + [2]*(len(u)-2) + [2] 
    u = np.append(u, -u[::-1], axis=0)
    codes += codes

    return mpath.Path(3*u, codes, closed=False)
    
    
  
def map_nwp_tracks_per_storm(date, outfile):
	#lat1 = -2
	#lat2 = -40
	#lon1=0
	#lon2=60
	
	#lat1 = 0
	#lat2 = -30
	#lon1=20
	#lon2=60
	
	lat1=-10
	lat2=-30
	lon1=18
	lon2=55
	
	fig = plt.figure(figsize=(6, 3))
	ax = fig.add_axes([0.05, 0.1, 0.9, 0.96])

	m = Basemap(llcrnrlon=lon1, llcrnrlat=lat2, urcrnrlon=lon2, urcrnrlat=lat1, projection='mill', resolution='l')

	RSMClats = [-40, 0, 0, -40]
	RSMClons = [30, 30, 90, 90]
	
	m.drawcoastlines(linewidth=0.5, color='k') #darkgray
	m.drawcountries(linewidth=0.5, color='k') #darkgray
	m.fillcontinents(color='silver',alpha=0.25)
	
	hurricane=get_hurricane_symbol()
	
	eps_data = np.genfromtxt(datadir+date+"_idai_ECMWF_eps_3way_matched_reformatted.txt", dtype=float, skip_header=1, usecols=np.arange(0,11))
	
	#tid = track id (column 1) in the eps reformatted files
	#1 = analysis track, 2 = best track, 3 = control, 4 - 53 = ensemble members 1-50
	
	#total length of the file, to loop over and find the individual tracks
	#probably a better way to do this, but it's friday 
	
	e = range(4,54)
	e.append(3)
	e.append(2)
	e.append(1)
	
	NP = len(eps_data[:,0])	
	
	print e
	
	for tid in e:
	
		
		if tid == 1:
			c = 'k'
			ls = '-'
			lw=0.7
			
		elif tid == 2:
			c = 'None'
			ls = '--'
			lw=0.75
			
		elif tid == 3:
			c = 'None' #'#FF9F1C'
			ls = '-'
			lw=0.75
		
		elif tid >= 4:
			c = 'lightblue'
			ls = '-'
			lw=0.5
			
		else:
			print "something went wrong with the tid"
			
	
		#tid = str(tid).zfill(4)
		#print float(tid)
		
		xyi = []
	
		for i in range(NP):
			#print eps_data[i,0]
			if eps_data[i,0] == float(tid):
				xyi.append(i)
		
		#print xyi	
		xs = []
		ys = []
		
		for i in xyi:
			xs.append(eps_data[i,7])
			ys.append(eps_data[i,8])
		
		
		if tid == 1:
			no_points = len(xs)
			points = np.array([xs, ys]).T.reshape(-1,1,2)
			segments = np.concatenate([points[:-1], points[1:]], axis=1)
			#colors = cm.twilight(np.linspace(0.5,1,no_points))
			cmap = LinearSegmentedColormap.from_list('mycmap',['lightgrey','black','dimgrey'])
			#cmap = LinearSegmentedColormap.from_list('mycmap',['plum','darkmagenta','plum'])
			colors = cmap(np.linspace(0,1,no_points))		

						
			for p, cl in zip(range(no_points - 1), colors):
				#put these here so as to not plot the analysis track before the forecast date, just from the fcst date onwards
				xarr=[]
				yarr=[]
				
				if p < 61: #61 is when it's no longer a TC, 37 is when it became a TC (10th March), I think
					z=10
				elif p > 61:#
					z=9
					
				if p > 37:	
					xarr.append(segments[p][0][0]) #shift these back an indent if plotting the *whole* analysis track
					xarr.append(segments[p][1][0])
					yarr.append(segments[p][0][1])
					yarr.append(segments[p][1][1])
				x,y = m(xarr,yarr)
				m.plot(x,y,linewidth=1.5, color=cl,zorder=z) #zorder=z

			pp=np.linspace(0,92,24)
			#colors2=cm.twilight(np.linspace(0.5,1,24))
			colors2=cmap(np.linspace(0,1,24))
			print pp			
			for p, cl in zip(pp, colors2):
				p=int(p)
				if p < 61:
					z=10
				elif p > 61:
					z=9
				if p < 37:
					xx,yy = m(xs[p], ys[p])
					#m.scatter(xx,yy,marker='o', edgecolors=cl,facecolors=cl, s=8,linewidth=1,zorder=z)
				elif 37 < p < 61:
					xx,yy = m(xs[p], ys[p])
					m.scatter(xx,yy,marker=hurricane, edgecolors=cl, facecolors='None', s=100, linewidth=1.25, zorder=z)
				elif p > 61:
					xx,yy = m(xs[p], ys[p])	
					m.scatter(xx,yy,marker='o',edgecolors=cl,facecolors='None',s=8,linewidth=1,zorder=z)	
			#xx,yy=m(xs[0],ys[0])
			#m.scatter(xx,yy,marker=hurricane, edgecolors=c, facecolors='None',s=25,linewidth=0.6,zorder=10)

		else:		
            		#continue
			x,y = m.shiftdata(xs,ys)
		
			#if tid == 3.:	
				#print xs
				#print ys
				#print x,y
		
			m.plot(x,y,linewidth=lw, color=c, linestyle=ls, latlon=True)
			#if tid == 1.:
				#xx,yy=m(xs[0],ys[0])
				#m.scatter(xx,yy,marker=hurricane, edgecolors=c, facecolors='None',s=25,linewidth=0.6,zorder=10)
			
			#if tid == 3:
				#xx,yy=m(xs[0],ys[0])
				#m.scatter(xx,yy,marker=hurricane, edgecolors=c, facecolors='None',s=25,linewidth=0.6,zorder=10)
		
		
	mean_file = datadir+date+"_idai_ECMWF_mean_3way_matched_reformatted.txt"
	if os.path.isfile(mean_file):	
		mean_data = np.genfromtxt(mean_file, dtype = float, skip_header=1, usecols=np.arange(0,11))
	
		x,y = m.shiftdata(mean_data[:,7],mean_data[:,8])
		m.plot(x,y, linewidth=0.75, color='mediumblue', latlon=True, linestyle='-',zorder=9)
		#xx,yy=m(mean_data[0,7],mean_data[0,8])
		#m.scatter(xx,yy,marker=hurricane, edgecolors='darkturquoise', facecolors='None',s=25,linewidth=0.6,zorder=10)
	
	det_file = datadir+date+"_idai_ECMWF_det_3way_matched_reformatted.txt"
	if os.path.isfile(det_file):
		det_data = np.genfromtxt(det_file, dtype=float, skip_header=1, usecols=np.arange(0,11))
	
		NP = len(det_data[:,0])
		#tid = str(3).zfill(4)
		
		tid = 3
		xyi = []
		for i in range(NP):
			if det_data[i,0] == float(tid):
				xyi.append(i)
	
		xs=[]
		ys=[]
		for i in xyi:
			xs.append(det_data[i,7])
			ys.append(det_data[i,8])
		
		x,y = m.shiftdata(xs,ys)
		m.plot(x,y, linewidth=0.75, color='#F71735', latlon=True, linestyle='-',zorder=9)
		#if len(xs)>0:
			#xx,yy=m(xs[0],ys[0])
			#m.scatter(xx,yy,marker=hurricane, edgecolors='#F71735', facecolors='None',s=25,linewidth=0.6,zorder=10)
		
	#else:
		#continue
	
	
	title="Cyclone Idai ECMWF Forecast\n"+str(date)
	ib = plt.Line2D((0, 1), (0, 0), color='k', linestyle='--',linewidth=0.5)
	an = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-',linewidth=0.5)
	det = plt.Line2D((0, 1), (0, 0), color='#F71735',linewidth=0.5)
	ctrl = plt.Line2D((0, 1), (0, 0), color='#FF9F1C',linewidth=0.5)
	mean = plt.Line2D((0, 1), (0, 0), color='medium',linewidth=0.5)
	eps = plt.Line2D((0, 1), (0, 0), color='lightblue',linewidth=0.5,linestyle=':')
	
	#legend = ax.legend((an, det, ctrl, mean, eps), ['Analysis Track', 'Deterministic', 'Control', 'Ensemble Mean', 'Ensembles'],title=title, fontsize=5, loc='lower left')
	#plt.setp(legend.get_title(), fontsize='5')
	#legend._legend_box.align = "left"

	#save and close the plot
	fig.subplots_adjust(wspace=0)
	plt.savefig(outfile, bbox_inches='tight', pad_inches=0.05, dpi=500)
	plt.close()


map_nwp_tracks_per_storm("2019031000", "idai_forecast_tracks_10March_noctrl.png")
	
#for d in range(22,29): #22,29 #first forecast date that something picked up Idai was 22nd Feb
	#day = str(d).zfill(2)
	
	#for h in [0,12]:
		#hr = str(h).zfill(2)
		
		#date = "201902"+day+hr
		#print date
		
		#outfile = "idai_ecmwf_forecast_tracks_"+date+".png"
		
		#map_nwp_tracks_per_storm(date, outfile)
		
		

#for d in range(1,22): #last forecast date with Idai was 21st March
	#day = str(d).zfill(2)
	
	#for h in [0,12]:
		#hr = str(h).zfill(2)
		
		#date = "201903"+day+hr
		#print date
		
		#outfile = "idai_ecmwf_forecast_tracks_"+date+".png"
		
		#map_nwp_tracks_per_storm(date, outfile)	
	
			
		
			
			
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
