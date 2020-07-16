import picsea_library as pl
import composite_functions_ukmo_nwp as cf
import datetime
import iris
import itertools
import numpy as np
import os
import warnings
import sys


#in the composite_functions_ukmo_nwp.py library, the functions are set up to run the following across the entire dataset
#those functions need to be updated with the latest resolutions and years as the precip dataset is updated


#print "starting months - forecasts'
#for month in [1,2,3,4,5,6,7,8,9,10,11,12]: #7 DIDNT WORK - july separate resolutions have same grid, somehow :(, will need to rerun when fixed
	#for lead in [0,1,2,3,4,5,6]:
		#print "month: ", month, "lead: ", lead
		#cf.composite_pcp_tot_month_separate_resolutions(month, lead)


#print "starting all - forecasts"
#for lead in [0,1,2,3,4,5,6]:
	#cf.composite_pcp_tot_all_separate_resolutions(lead)

#for lead in [0,1,2,3,4, 5,6]:
	#cf.composite_pcp_tot_seasons_res(lead)

#print 'starting months - TRMM'
#for month in [1,2,3,4,5,6,7,8,9,10,11,12]:
	#cf.composite_pcp_tot_month_separate_resolutions_trmm(month)

#print 'starting all - TRMM'
#cf.composite_pcp_tot_all_separate_resolutions_trmm()
#cf.composite_pcp_tot_seasons_res_trmm()

#the next couple of lines will calculate the 2006-2016 TRMM precip composites using the analysis tracks rather than ibtracs tracks
#the functions in composite_functions_ukmo_nwp will need to be updated using the comments, when updating analysis to go up to 2018
#all the functions will need updating to bring the analysis up to 2019...
#for month in [1,2,3,4,5,6,7,8,9,10,11,12]:
	#cf.composite_pcp_tc_month_separate_resolutions_trmm_analysis(month)
	
#cf.composite_pcp_tc_all_separate_resolutions_trmm_analysis()



#trmm_tot_dir="/gws/nopw/j04/klingaman/emerton/total_trmm_precip_composites/"
##last input file will eventually be -062019:
#trmm_tot_infiles=[trmm_tot_dir+"trmm.comp_pcp_tot.072006-032010_regridded_n320.nc", trmm_tot_dir+"trmm.comp_pcp_tot.032010-072014_regridded_n320.nc", trmm_tot_dir+"trmm.comp_pcp_tot.072014-072017_regridded_n320.nc", trmm_tot_dir+"trmm.comp_pcp_tot.072017-062018_regridded_n320.nc"]
#trmm_tot_outfile=trmm_tot_dir+"trmm.comp_pcp_tot.072006-062018_n320.nc"
#pl.add_files(trmm_tot_infiles,trmm_tot_outfile)

#trmm_tc_dir="/gws/nopw/j04/klingaman/emerton/ibtracs_precip_composites/"
##last input file will eventually be -062019:
#trmm_tc_infiles=[trmm_tc_dir+"trmm.comp_pcp_tc.072006-032010_regridded_n320.nc", trmm_tc_dir+"trmm.comp_pcp_tc.032010-072014_regridded_n320.nc", trmm_tc_dir+"trmm.comp_pcp_tc.072014-072017_regridded_n320.nc", trmm_tc_dir+"trmm.comp_pcp_tc.072017-062018_regridded_n320.nc"]
#trmm_tc_outfile=trmm_tc_dir+"trmm.comp_pcp_tc.072006-122016_n320.nc"
#pl.add_files(trmm_tc_infiles,trmm_tc_outfile)

#####################################################################################################################
#THIS BIT NEEDS MODIFYING AND RUNNING FOR THE TRMM ANALYSIS TRACK COMPOSITES!
#this, again, will be used to calculate the 2006-2016 TRMM precip composites using the analysis tracks rather than ibtracs tracks
#this will also need updating when eventually update the analysis to 2018 or 2019
trmm_tc_dir="/gws/nopw/j04/klingaman/emerton/analysis_trmm_tc_precip_composites/"
#last input file will eventually be -062019:
trmm_tc_infiles=[trmm_tc_dir+"trmm.comp_pcp_tc.analysis.072006-032010_regridded_n320.nc", trmm_tc_dir+"trmm.comp_pcp_tc.analysis.032010-072014_regridded_n320.nc", trmm_tc_dir+"trmm.comp_pcp_tc.analysis.072014-122016_regridded_n320.nc"] #, trmm_tc_dir+"trmm.comp_pcp_tc.072014-072017_regridded_n320.nc", trmm_tc_dir+"trmm.comp_pcp_tc.072017-062018_regridded_n320.nc"
trmm_tc_outfile=trmm_tc_dir+"trmm.comp_pcp_tc.072006-122016_n320.nc"
pl.add_files(trmm_tc_infiles,trmm_tc_outfile)
#####################################################################################################################


#for lead in [0,1,2,3,4,5,6]:
	#ukmo_tot_dir="/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_TOTAL_precip_composites/"+str(lead)+"/"
	#ukmo_tot_infiles=[ukmo_tot_dir+"ukmo_nwp.comp_pcp_tot."+str(lead)+"_days.072006-032010.n320.nc",ukmo_tot_dir+"ukmo_nwp.comp_pcp_tot."+str(lead)+"_days.032010-072014.regridded_n320.nc", ukmo_tot_dir+"ukmo_nwp.comp_pcp_tot."+str(lead)+"_days.072014-072017.regridded_n320.nc", ukmo_tot_dir+"ukmo_nwp.comp_pcp_tot."+str(lead)+"_days.072017-062018.regridded_n320.nc"]
	#print ukmo_tot_infiles
	#ukmo_tot_outfile=ukmo_tot_dir+"ukmo_nwp.comp_pcp_tot."+str(lead)+"_days.072006-062018.n320.nc"
	#pl.add_files(ukmo_tot_infiles,ukmo_tot_outfile,cube_count=4)

#for lead in [0,1,2,3,4,5,6]:
	#ukmo_tc_dir="/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_precip_composites/"+str(lead)+"/"
	#ukmo_tc_infiles=[ukmo_tc_dir+"ukmo_nwp.comp_pcp_tc."+str(lead)+"_days.072006-032010.n320.nc",ukmo_tc_dir+"ukmo_nwp.comp_pcp_tc."+str(lead)+"_days.032010-072014.regridded_n320.nc", ukmo_tc_dir+"ukmo_nwp.comp_pcp_tc."+str(lead)+"_days.072014-072017.regridded_n320.nc", ukmo_tc_dir+"ukmo_nwp.comp_pcp_tc."+str(lead)+"_days.072017-062018.regridded_n320.nc"]
	#ukmo_tc_outfile = ukmo_tc_dir + "ukmo_nwp.comp_pcp_tc." + str(lead) + "_days.072006-062018.n320.nc"
	#pl.add_files(ukmo_tc_infiles,ukmo_tc_outfile,cube_count=4)

#for season in ['DJFM', 'NDJFMA', 'MJJASO']:
	#trmm_tot_dir = "/gws/nopw/j04/klingaman/emerton/total_trmm_precip_composites/"
	#trmm_tot_infiles=[trmm_tot_dir+"trmm.comp_pcp_tot.072006-032010."+season+"_regridded_n320.nc",trmm_tot_dir+"trmm.comp_pcp_tot.032010-072014."+season+"_regridded_n320.nc",trmm_tot_dir+"trmm.comp_pcp_tot.072014-072017."+season+"_regridded_n320.nc",trmm_tot_dir+"trmm.comp_pcp_tot.072017-062018."+season+"_regridded_n320.nc" ]
	#trmm_tot_outfile=trmm_tot_dir+"trmm.comp_pcp_tot.072006-062018."+season+"_n320.nc"
	#pl.add_files(trmm_tot_infiles,trmm_tot_outfile)

	#trmm_tc_dir = "/gws/nopw/j04/klingaman/emerton/ibtracs_precip_composites/"
	#trmm_tc_infiles = [trmm_tc_dir + "trmm.comp_pcp_tc.072006-032010." + season + "_regridded_n320.nc",
	                  #  trmm_tc_dir + "trmm.comp_pcp_tc.032010-072014." + season + "_regridded_n320.nc",
	                   # trmm_tc_dir + "trmm.comp_pcp_tc.072014-072017." + season + "_regridded_n320.nc",trmm_tc_dir + "trmm.comp_pcp_tc.072017-062018." + season + "_regridded_n320.nc"]
	#trmm_tc_outfile = trmm_tc_dir + "trmm.comp_pcp_tc.072006-062018." + season + "_n320.nc"
	#pl.add_files(trmm_tc_infiles, trmm_tc_outfile)

	#for lead in [0, 1, 2, 3, 4, 5, 6]:
		#ukmo_tot_dir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_TOTAL_precip_composites/" + str(lead) + "/"
		#ukmo_tot_infiles = [ukmo_tot_dir + "ukmo_nwp.comp_pcp_tot." + str(lead) + "_days.072006-032010.n320."+season+".nc",
		                    #ukmo_tot_dir + "ukmo_nwp.comp_pcp_tot." + str(lead) + "_days.032010-072014."+season+".regridded_n320.nc",
		                    #ukmo_tot_dir + "ukmo_nwp.comp_pcp_tot." + str(lead) + "_days.072014-072017."+season+".regridded_n320.nc", ukmo_tot_dir + "ukmo_nwp.comp_pcp_tot." + str(lead) + "_days.072017-062018."+season+".regridded_n320.nc"]
		#print ukmo_tot_infiles
		#ukmo_tot_outfile = ukmo_tot_dir + "ukmo_nwp.comp_pcp_tot." + str(lead) + "_days.072006-062018."+season+".n320.nc"
		#pl.add_files(ukmo_tot_infiles, ukmo_tot_outfile, cube_count=4)

		#ukmo_tc_dir = "/gws/nopw/j04/klingaman/emerton/ukmo_nwp_analysis/ukmo_nwp_precip_composites/" + str(lead) + "/"
		#ukmo_tc_infiles = [ukmo_tc_dir + "ukmo_nwp.comp_pcp_tc." + str(lead) + "_days.072006-032010.n320."+season+".nc",
		                   #ukmo_tc_dir + "ukmo_nwp.comp_pcp_tc." + str(lead) + "_days.032010-072014."+season+".regridded_n320.nc",
		                   #ukmo_tc_dir + "ukmo_nwp.comp_pcp_tc." + str(lead) + "_days.072014-072017."+season+".regridded_n320.nc",ukmo_tc_dir + "ukmo_nwp.comp_pcp_tc." + str(lead) + "_days.072017-062018."+season+".regridded_n320.nc"]
		#ukmo_tc_outfile = ukmo_tc_dir + "ukmo_nwp.comp_pcp_tc." + str(lead) + "_days.072006-062018."+season+".n320.nc"
		#pl.add_files(ukmo_tc_infiles, ukmo_tc_outfile, cube_count=4)















