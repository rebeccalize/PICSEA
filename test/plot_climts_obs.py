#!/usr/bin/env python

## To run:
#python plot_climts_obs.py -1 1980 1 2 1989 1 31

import sys
from netCDF4 import Dataset,num2date
import numpy as np
import time
import datetime
import pandas as pd

import matplotlib
matplotlib.use('PS') 

import matplotlib.dates as mdates
import matplotlib.pyplot as plt


def datespan(startDate, endDate, delta=datetime.timedelta(days=1)):
    currentDate = startDate
    while currentDate < endDate:
        yield currentDate
        currentDate += delta

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

st2plot=sys.argv[1]
print "Station=",st2plot


stats='/home/mo/more/GLOFAS/stations_info_MAPPED_050_cal_NEW_man.csv'
print stats

df_sta = pd.read_csv(stats)
statsid=df_sta['id']
print len(statsid)
#print statsid
#print statsid[0]


#print(rundate2plot.strftime("%Y-%m-%d 00:00:00"))
rundate2plot = datetime.datetime(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]))
enddate2plot = datetime.datetime(int(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7]))

daynum=0
for dd in datespan(datetime.datetime(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])),datetime.datetime(int(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7]))+datetime.timedelta(days=1),delta=datetime.timedelta(days=1)):
    daynum=daynum+1
print "daynum=",daynum

#enddate2plot = rundate2plot + datetime.timedelta(hours=int(sys.argv[5]))
print(rundate2plot)
print(enddate2plot)

#fsim1='/scratch/mo/moi/cama_g9y8/outflw_grdc_19790101_19800101.nc'
#fsim1='/scratch/mo/moi/cama_g9y8/outflw_grdc_19790101_19800401.nc'
fsim1='/home/mo/more/GLOFAS/outflw_grdc_19790101_19890101.nc'
fsim2='/home/mo/moi/work/surface_model/TIGGE_verif/rivout_grdc_g4bl_19800101_20140101.nc'
fobs='/home/mo/moi/work/surface_model/GRDC/World_Discharge_JRC/cal_2013_grdcDisData.nc'
#fobsnew='/home/mo/moi/work/surface_model/GRDC/netcdfs/discharge_obs_daily_grdc_mergedwith_oldcama_19800101_20131231.nc'



#load stuff
nc = Dataset(fsim1,'r')
time_sim1=num2date(nc.variables['time'],getattr(nc.variables['time'],'units'))
data_sim1=nc.variables['outflw'][:]
nc.close()

nc = Dataset(fsim2,'r')
time_sim2=num2date(nc.variables['time'],getattr(nc.variables['time'],'units'))
data_sim2=nc.variables['rivout'][:]
nc.close()

nc = Dataset(fobs,'r')
time_obs=num2date(nc.variables['time'],getattr(nc.variables['time'],'units'))
data_obs=nc.variables['dis'][:]
grdcid=nc.variables['grdcID'][:]
river=nc.variables['Criver'][:]
station=nc.variables['Cstation'][:]
country=nc.variables['Ccountry'][:]
lat=nc.variables['lat'][:]
lon=nc.variables['lon'][:]
nc.close()


#FIND STATION IN SIM1 (349 stations)
print len(statsid)
if int(st2plot) > 500:
    for i in range(len(statsid)):
#        print i,st2plot,statsid[i]
        if int(st2plot) == statsid[i]:
           stsim1=i
else:
    stsim1=int(st2plot)
print "StationID in SIM1= ",stsim1 #int(st2plot)-1



#FIND STATION IN SIM2 and OBS (413 stations)
print len(grdcid)
if int(st2plot) > 500:
    for i in range(len(grdcid)):
#        print i,st2plot,grdcid[i]
        if int(st2plot) == grdcid[i]:
           stobs=i
else:
    stobs=int(st2plot)
print "StationID in SIM2/OBS= ",stobs #int(st2plot)-1

if int(st2plot) == -1:
	stsim1=-1
	stobs=-1






yyy=int(sys.argv[2])
mmm=int(sys.argv[3])
ddd=int(sys.argv[4])
ddd1=ddd-1
if ddd1 < 1:
    mmm1=mmm-1
    if mmm1 < 1:
        yyy=yyy-1
        mmm=12
        ddd=31
    else:
        mmm=mmm1
        if mmm == 1 or mmm == 3 or mmm == 5 or mmm == 7 or mmm == 8 or mmm == 10:
            ddd=31
        elif mmm == 4 or mmm == 6 or mmm == 9 or mmm == 11:
            ddd=30
        else:
            if yyy == 2008:
                ddd=29
            else:
                ddd=28
else:
    ddd=ddd1
print int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
print yyy,mmm,ddd

print len(data_obs)
d=datetime.datetime(1979,12,31)
for k in range(len(data_obs)):
    d = d + datetime.timedelta(days=1)
    if d == datetime.datetime(yyy,mmm,ddd):
        obsstart=k
print 'obsstart= ',obsstart



print len(data_sim1)
d=datetime.datetime(1979,1,1)
for k in range(len(data_sim1)):
    d = d + datetime.timedelta(days=1)
    if d == datetime.datetime(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])):
        sim1start=k
print 'sim1start= ',sim1start



print len(data_sim2)
d=datetime.datetime(1980,1,1)
for k in range(len(data_sim2)):
    d = d + datetime.timedelta(days=1)
    if d == datetime.datetime(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])):
        sim2start=k
print 'sim2start= ',sim2start



print "Starts=",obsstart,sim1start,sim2start



#simstart=0
#obsstart=0

dates=[]
num=daynum #3500 #len(time_sim1)
ss1 = np.zeros(num)
ss2 = np.zeros(num)
sob = np.zeros(num)
for t in range(num):
    dates.append(datetime.datetime(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]))+datetime.timedelta(days=t))
    if stobs == -1:
        tmin=0
        tmax=349
    else:
        tmin=stobs
        tmax=stobs+1

    if stobs == -1:
        sm1=0
        sm2=0
        ob=0
        cnt=0
        for i in range(tmin,tmax):
#            print len(grdcid)
            found=0
            for j in range(len(grdcid)):
                #print j,st2plot,grdcid[j]
                if statsid[i] == grdcid[j]:
                    stt=j
                    found=1
            if found == 0:
#                print "ERROR",i,statsid[i]
                continue
#                sys.exit()
#            else:
#                print i,statsid[i],stt
                    
            if data_obs[obsstart+t,stt] >= 0 and data_obs[obsstart+t,stt] < 1000000:
                sm1=sm1+data_sim1[sim1start+t,i]
                sm2=sm2+data_sim2[sim2start+t,stt]
                ob=ob+data_obs[obsstart+t,stt]
                cnt=cnt+1
        print "cnt=",t,cnt
        if cnt == 0:
            ss1[t] = np.nan
            ss2[t] = np.nan
            sob[t] = np.nan
        else:
            ss1[t] = sm1/float(cnt)
            ss2[t] = sm2/float(cnt)
            sob[t] = ob/float(cnt)
    else:
        ob=0
        cnt=0
        for i in range(tmin,tmax):
            if data_obs[obsstart+t,i] >= 0 and data_obs[obsstart+t,i] < 1000000:
                ob=ob+data_obs[obsstart+t,i]
                cnt=cnt+1
        if cnt == 0:
            sob[t] = np.nan
        else:
            sob[t] = ob/float(cnt)

        found=0
        for j in range(len(grdcid)):
            #print j,st2plot,grdcid[j]
            if statsid[i] == grdcid[j]:
                stt=j
                found=1
        if found == 0:
            print "ERROR",i,j
            continue
#            sys.exit()
#        else:
#            print i,j

        sm1=0
        sm2=0
        cnt=0
        for i in range(tmin,tmax):
            sm2=sm2+data_sim2[sim2start+t,i]
            cnt=cnt+1
        if cnt == 0:
            ss2[t] = np.nan
            sob[t] = np.nan
        else:
            ss2[t] = sm2/float(cnt)
            sob[t] = ob/float(cnt)

        for i in range(stt,stt):
            sm1=sm1+data_sim1[sim1start+t,i]
            cnt=cnt+1
        if cnt == 0:
            ss1[t] = np.nan
        else:
            ss1[t] = sm1/float(cnt)


    print t,tmin,tmax,cnt
print "comp done"

print dates

plt.plot(dates,ss1[:],'r',lw=1,label='EI+GPCP2 NEW CAMA') #g35p')
plt.plot(dates,ss2[:],'g',lw=1,label='EI+GPCP2 OLD CAMA') #g4bl')
plt.plot(dates,sob[:],'b',lw=1,label='GRDC')
#plt.plot(time_sim,data_sim[:,stsim],'k',lw=6,ls="-.",label='Long ERA')






#plt.plot(time_sim1,data_sim1[:,stsim1],'r',lw=4,label='EI+GPCP2 NEW CAMA') #g35p')
#plt.plot(time_sim2,data_sim2[:,stobs],'g',lw=4,label='EI+GPCP2 OLD CAMA') #g4bl')
#plt.plot(time_obs,data_obs[:,stobs],'b',lw=4,label='GRDC')
##plt.plot(time_sim,data_sim[:,stsim],'k',lw=6,ls="-.",label='Long ERA')
#
#print(r'$\alpha_i > \beta_i$')
#plt.ylabel("m$\3_i$ /s")
plt.ylabel(r'$m^{3}/s$')

#print(data_sim[:,stsim])


#print(rundate2plot.year)

#plt.xlim(datetime.datetime(1980,9,1),datetime.datetime(1980,9,30))
#plt.xlim(datetime.datetime(2011,3,11),datetime.datetime(2011,3,21))
plt.xlim(datetime.datetime(rundate2plot.year,rundate2plot.month,rundate2plot.day),datetime.datetime(enddate2plot.year,enddate2plot.month,enddate2plot.day))
#plt.xlim(datetime.datetime(1980,1,1),datetime.datetime(2013,12,31))
#plt.xlim(datetime.datetime(2013,5,30),datetime.datetime(2013,6,8))


print station[stobs][0:30]
print river[stobs][0:30]
print country[stobs][0:30]

s=""
for i in station[stobs][0:6]:
    s+=str(i)
print(s)

s+=str("  ")
for i in river[stobs][0:5]:
    s+=str(i)
print(s)

s+=str("  ")
for i in country[stobs][0:2]:
    s+=str(i)
print(s)

s+=str("  ")
s+=str(lat[stobs])
s+=str(" ")
s+=str(lon[stobs])
print(s)

#s+=str("  RT: 20090210 00 UTC")

title=s
plt.title(title, weight = 'bold')



#print len(data_obs)
#d=datetime.datetime(1979,12,31)
##d=datetime.datetime(2008,12,31)
#for k in range(len(data_obs)):
#    d = d + datetime.timedelta(days=1)
#    if d == datetime.datetime(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])):
#        obsstart=k
#print 'obsstart= ',obsstart


print len(data_sim1)
d=datetime.datetime(1978,12,31)
for k in range(len(data_sim1)):
    d = d + datetime.timedelta(days=1)
    if d == datetime.datetime(int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])):
        sim1start=k
print 'sim1start= ',sim1start


#min=1000000
#for x in sob[obsstart:obsstart+daynum]:
#    if x < min:
#        min=x
#print(min)
#for x in ss1[simstart:simstart+daynum]:
#    if x < min:
#        min=x
#for x in ss2[simstart:simstart+daynum]:
#    if x < min:
#        min=x
#print(min)
   

#max=-100000
#for x in sob[obsstart:obsstart+daynum]:
#    if x > max:
#        max=x
#print(max)
#for x in ss1[simstart:simstart+daynum]:
#    if x > max:
#        max=x
#for x in ss2[simstart:simstart+daynum]:
#    if x > max:
#        max=x
#print(max)



#ymin = min([min(y_list) for y_list in data_for1[:,st]])
#ymax = max([max(y_list) for y_list in y_list_of_lists])


#plt.legend(loc='upper left')

#plt.plot(legend=False)

#plt.fmt_xdata = mdates.DateFormatter('%Y-%m-%d')

#plt.autofmt_xdate()

#ymin=min-500
#if ymin < 0:
#	ymin=0
#ymax=max+max/float(4)
#plt.ylim((ymin,ymax))
plt.ylim((0,15000))
plt.legend(loc='upper left')

plt.gcf().autofmt_xdate()

plt.show()
print "StationID in SIM1= ",stsim1 #int(st2plot)-1
print "StationID in SIM2/OBS= ",stobs #int(st2plot)-1


print '/home/mo/more/GLOFAS/comp_' + st2plot + '.ps'
plt.savefig('/home/mo/more/GLOFAS/comp_' + st2plot + '.ps')
#print './comp2.ps'
#plt.savefig('/home/mo/more/GLOFAS/Qcomp3.ps')
