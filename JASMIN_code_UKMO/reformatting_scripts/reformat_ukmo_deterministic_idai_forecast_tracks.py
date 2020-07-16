import picsea_library as pl

#Reformat a file containing ibtracs track data so that it's easier for python to read
#Note: for the SH, track files run from July to July, with the filename the name of the start year

datadir = "/gws/nopw/j04/klingaman/emerton/TRACK/UKMO/Y2019-SH/"
savedir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/reformatted_idai_forecast_tracks/"


#dates = [2019030412, 2019030500, 2019030512,2019030600,2019030612,2019030700,2019030800,2019030812,2019030900,2019030912,2019031000,2019031012,2019031100,2019031112,2019031212,2019031300,2019031312,2019031400]
dates = [2019030712, 2019031200]	
for date in dates:

	infile = datadir+"idai_UKMO_"+str(date)+".txt"
	outfile = savedir+"idai_UKMO_"+str(date)+"_reformatted.txt"
	
	pl.rewrite_track_file(infile,outfile,'ukmo_nwp')
	
	infile = savedir+"idai_UKMO_"+str(date)+"_reformatted.txt"
	outfile = savedir+"idai_UKMO_"+str(date)+"_interpolated.txt"

	pl.interpolate_track_file(infile,outfile)


