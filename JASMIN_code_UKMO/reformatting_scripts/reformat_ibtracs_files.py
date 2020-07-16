import picsea_library as pl

#Reformat a file containing ibtracs track data so that it's easier for python to read
#Note: for the SH, track files run from July to July, with the filename the name of the start year

datadir = "/gws/nopw/j04/klingaman/emerton/TRACK/ANALYSES/"
savedir = "/gws/nopw/j04/klingaman/emerton/analysis_reformatted_track_files_yearmonth/"


#for year in [2016, 2017]: #2006,2016
	#infile = datadir+"storms_SH"+str(year)
	#outfile = savedir+"storms_SH"+str(year)+"_reformatted.txt"

	#pl.rewrite_track_file(infile, outfile, 'ukmo_nwp')
	

#for year in [2016, 2017]: #2006,2016
	#infile = savedir+"storms_SH"+str(year)+"_reformatted.txt"
	#outfile = savedir+"storms_SH"+str(year)+"_interpolated.txt"
	
	#pl.interpolate_track_file(infile,outfile)
	
year1s = [2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017]
year2s = [2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018]

for y1,y2 in zip(year1s,year2s):

	print y1,y2,' reformatting'
	infile = datadir+"OP_JUL-JUN"+str(y1)+str(y2)+"_VOR_VERTAVG_T63/ff_trs_neg.2day.addT63vor_addwind10m_addmslp_dummy.new.TC.ibtref"
	outfile = savedir+"analysis."+str(y1)+str(y2)+"reformatted.txt"
	
	pl.rewrite_track_file(infile,outfile,'analysis')
	
for y1,y2 in zip(year1s,year2s):

	print y1,y2,' interpolating'
	infile = savedir+"analysis."+str(y1)+str(y2)+"reformatted.txt"
	outfile = savedir+"analysis."+str(y1)+str(y2)+"interpolated.txt"
	
	pl.interpolate_track_file(infile,outfile)


#REFORMAT THE 2018-2019 (JULY - APRIL) ANALYSIS TRACKS (matched to ibtracs) for SEARCH PROPOSAL

#datadir = "/gws/nopw/j04/klingaman/emerton/TRACK/ANALYSES/OP_JUL-APR20182019_VOR_VERTAVG_T63/"
#savedir = "/gws/nopw/j04/klingaman/emerton/TRACK/ANALYSES/OP_JUL-APR20182019_VOR_VERTAVG_T63/"


#infile = datadir + "kenneth_analysis_track.txt"
#outfile = savedir + "kenneth_analysis_track_reformatted.txt"

#pl.rewrite_track_file(infile, outfile, 'ukmo_nwp')


#infile = savedir + "kenneth_analysis_track_reformatted.txt"
#outfile = savedir + "kenneth_analysis_track_interpolated.txt"

#pl.interpolate_track_file(infile, outfile)

