import picsea_library as pl
import datetime
import iris
import itertools
import numpy as np
import os

IBTRACS_DIR = "/gws/nopw/j04/klingaman/emerton/ibtracs_reformatted_data_and_track_maps/"
IBTRACS_FILE = "storms_SH%d%d_interpolated.txt"
COMP_PCP_TC_IBTRACS_FILE = "ibtracs.comp_pcp_tc.%d%02d.nc"

COMP_PCP_TC_IBTRACS_DIR = "/gws/nopw/j04/klingaman/emerton/ibtracs_precip_composites_wcont_days/"
if not os.path.exists(COMP_PCP_TC_IBTRACS_DIR):
    os.makedirs(COMP_PCP_TC_IBTRACS_DIR)


#pcp = pl.read_trmm_3b42(2012, 11, 10, 3)
#print pcp

def trmm_3b42_pcp_accumulation(year, month, day, hour, lon=None, lat=None,
                           quiet=False):
    """Returns a map of 6-hour accumulated precipitation in mm from TRMM 3B42,
    over the 6-hour period CENTRED ON the given time in UTC.  Hour must be one
    of [3, 9, 15, 21].

    Values more than 5 degrees from the given lat/lon location set to zero (if
    lat/lon are given).

    3B42 is provided every 3 hours as a mean over the 3-hour period CENTERED on
    the given time.  Hence, if hour=9 then we sum:
          TRMM 06UTC file x1.5
        + TRMM 09UTC file x3.0
        + TRMM 12UTC file x1.5
    """

    # Check we have a valid time (same as track points)
    valid_hours = range(3, 22, 6)
    if hour not in valid_hours:
        raise ValueError('Invalid hour %s; must be one of %s' % (hour,
                                                                 valid_hours))

    # Iterate for each time
    t = datetime.datetime(year, month, day, hour)
    dt = datetime.timedelta(hours=3)
    pcp_list = iris.cube.CubeList()
    for mult, time in itertools.izip([1.5, 3., 1.5], [t-dt, t, t+dt]):
        if not quiet:
            print '   ', time.strftime('%Y/%m/%d %H:%M'), 'x%.01f' % mult

        # Open TRMM 3B42 file for given time
        pcp = pl.read_trmm_3b42(time.year, time.month, time.day, time.hour)

        ##################################################################
        # Set to zero outside 5 degree radius
        # Tracks use longitude values 0 -> 360 but TRMM uses -180 -> 180
        ##################################################################
        if lat is not None and lon is not None:
            lats = pcp.coord(axis='Y').points
            lons1 = pcp.coord(axis='X').points
            lons = np.where(lons1<0, lons1+360., lons1)
            lons_, lats_ = np.meshgrid(lons, lats)
            pcp.data = np.ma.where(np.hypot(lons_-lon, lats_-lat) > 5, 0.,
                                   pcp.data*mult)
        else:
            pcp.data = pcp.data*mult
        pcp_list.append(pcp)

    # Return sum
    tot_pcp = pl.add_cubes(pcp_list, contributing_days=False)
    tot_pcp.units = 'mm'
    tot_pcp.standard_name = 'lwe_thickness_of_precipitation_amount'
    tot_pcp.long_name = 'accumulated precipitation'
    tot_pcp.var_name = 'pcp'
    return tot_pcp

def composite_pcp_tc_ibtracs_year_month(yyyy, month):
                                
    """Computes a composite of the precipitation due to (actual) TCs for a particular
    year and month - using ibtracs TCs, NOT NWP!
    
    Only the most basic compositing function - see other functions for including the phase of the MJO (to composite only those TCS whose genesis was in a certain phase, or to select a specific lead time, or to use 3-way matched tracks """

    print yyyy
    print month
    if month in [7,8,9,10,11,12]:
        year = yyyy
    elif month in [1,2,3,4,5,6]:
        year = yyyy-1

    infile = os.path.join(IBTRACS_DIR, IBTRACS_FILE % (year, year+1))
    print "using file: ", infile
    comp_file = COMP_PCP_TC_IBTRACS_FILE % (yyyy, month)

        
    outdir = COMP_PCP_TC_IBTRACS_DIR
    outfile = os.path.join(outdir, str(yyyy), comp_file)
    
    if os.path.isfile(outfile):
        raise ValueError('Output file %s already exists' % outfile)
        
    odir = os.path.split(outfile)[0]
    if not os.path.isdir(odir):
        os.makedirs(odir)

    exclude = []  # if using a lead time, would need to use the pl.exclude_days function here

    # Count number of days "included" (this includes days for which no tracks
    # exist - i.e., count how many days are not positively excluded)
    
    #t1 = first day of the month
    #t2 = last day of the month 
    #(essentially counting the number of days in the month)
    t1 = datetime.datetime(yyyy, month, 1, 0)
    dt = datetime.timedelta(days=1)
    count_days = 0
    #not sure why we care specifically about july 2017? is this when the dataset ended?
    #if the yearmonth is 201707, set t2 to 20170711... 
    #if (year, month) == (2017, 7):
        #t2 = datetime.datetime(2017, 7, 11)
    #account for December needing to use 1st Jan to find the last date
    if month == 12:
        t2 = datetime.datetime(yyyy+1, 1, 1) - dt
    #otherwise, set t2 to the last day of the month
    else:
        t2 = datetime.datetime(yyyy, month+1, 1) - dt

    if (yyyy, month) == (2010, 3):
        count_days_res = {'n320': 0, 'n512': 0}
        start_n512 = datetime.datetime(2010, 3, 9)  # Actually at mid-day
    elif (yyyy, month) == (2014, 7):
        count_days_res = {'n512': 0, 'n768': 0}
        start_n768 = datetime.datetime(2014, 7, 15)  # Actually at mid-day
    while t1 <= t2:
        if t1 not in exclude:
            if t1.timetuple()[:3] == (2017, 7, 11):
                count_days += 0.5
            else:
                count_days += 1
            if (yyyy, month) == (2010, 3):
                if t1 < start_n512:
                    count_days_res['n320'] += 1
                elif t1 == start_n512:
                    count_days_res['n320'] += 0.5
                    count_days_res['n512'] += 0.5
                else:
                    count_days_res['n512'] += 1
            elif (yyyy, month) == (2014, 7):
                if t1 < start_n768:
                    count_days_res['n512'] += 1
                elif t1 == start_n768:
                    count_days_res['n512'] += 0.5
                    count_days_res['n768'] += 0.5
                else:
                    count_days_res['n768'] += 1
        t1 += dt
    # Get year, month, day, hour, lon, lat from file
    ain = np.genfromtxt(infile, dtype=float, skip_header=1, usecols=[0]+range(3, 9))

    # Select track points in the correct year and month
    ain_y = ain[np.where(ain[:, 1] == yyyy)]
    ain_ym = ain_y[np.where(ain_y[:, 2] == month)]
    ain_use = ain_ym
          
    # Iterate for each point
    pcp_list = iris.cube.CubeList()
    if (yyyy, month) == (2010, 3):
        pcp_list_res = {'n320': iris.cube.CubeList(),
                        'n512': iris.cube.CubeList()}
    elif (yyyy, month) == (2014, 7):
        pcp_list_res = {'n512': iris.cube.CubeList(),
                        'n768': iris.cube.CubeList()}

    excluded = []
    for t_id, y, m, d, h, lon, lat in ain_use:
        t_id, y, m, d, h = [int(x) for x in t_id, y, m, d, h]

        this_t_h = datetime.datetime(y, m, d, h)
        this_t = datetime.datetime(y, m, d)

        if this_t in exclude:
            if this_t not in excluded:
                print '%02d/%02d/%02d -- EXCLUDE' % (y, m, d)
                excluded.append(this_t)
            continue

        print '%02d/%02d/%02d %02d:00' % (y, m, d, h)
        


        # Get precipitation
        this_pcp = trmm_3b42_pcp_accumulation(y, m, d, h, lon, lat)
        pcp_list.append(this_pcp)

        # For months including a change of resolution, append Cube to
        # the appropriate CubeList
        if (yyyy, month) == (2010, 3):
            if datetime.datetime(y, m, d, h) < start_n512:
                pcp_list_res['n320'].append(this_pcp)
            else:
                pcp_list_res['n512'].append(this_pcp)
        elif (yyyy, month) == (2014, 7):
            if datetime.datetime(y, m, d, h) < start_n768:
                pcp_list_res['n512'].append(this_pcp)
            else:
                pcp_list_res['n768'].append(this_pcp)

        # Check for negative values
        if pcp_list[-1].data.min() < 0:
            raise ValueError('Precip contains negative values')

    # Save
    if len(pcp_list):
        pcp = pl.add_cubes(pcp_list, contributing_days=False)
    else:
        print 'No data - creating dummy map'
        dummy = trmm_3b42_pcp_accumulation(1998, 1, 1, 3)
        pcp = pl.add_cubes([dummy, dummy], contributing_days=False)
        pcp.data = np.zeros_like(pcp.data)
    pcp.units = 'mm'
    pcp.standard_name = 'lwe_thickness_of_precipitation_amount'
    pcp.long_name = 'accumulated precipitation'
    pcp.var_name = 'pcp'
    pcp.attributes['contributing_days'] = float(count_days)
    iris.save(pcp, outfile)
    print outfile

    # For months with more than one resolution, sum separately
    if (yyyy, month) in [(2010, 3), (2014, 7)]:
        for k in pcp_list_res.iterkeys():
            if len(pcp_list_res[k]):
                pcp_k = pl.add_cubes(pcp_list_res[k], contributing_days=False)
            else:
                print 'No data - creating dummy map'
                dummy = trmm_3b42_pcp_accumulation(1998, 1, 1, 3)
                pcp_k = pl.add_cubes([dummy, dummy], contributing_days=False)
                pcp_k.data = np.zeros_like(pcp.data)

            # Set metadata
            pcp_k.units = 'mm'
            pcp_k.standard_name = 'lwe_thickness_of_precipitation_amount'
            pcp_k.long_name = 'accumulated precipitation'
            pcp_k.var_name = 'pcp'
            pcp_k.attributes['contributing_days'] = float(count_days_res[k])

            # Save
            outfile_k = outfile.replace('.nc', '.%s.nc' % k)
            iris.save(pcp_k, outfile_k)
            print outfile_k

def composite_pcp_tc_year(year1, year2, lead=None):
    """Computes the total composite precipitation due to observed TCs for a particular
    TC season (July year1 - June year2)."""

    infiles = []
    for m in [7, 8, 9, 10, 11, 12]:  # 7,
        # m = str(m).zfill(2)
        filename = COMP_PCP_TC_IBTRACS_FILE % (year1, m)
        infile = os.path.join(COMP_PCP_TC_IBTRACS_DIR, str(year1), filename )
        infiles.append(infile)

    for m in [1, 2, 3, 4, 5, 6]:
        # m=str(m).zfill(2)
        # year=year2
        filename = COMP_PCP_TC_IBTRACS_FILE % (year2, m)
        infile = os.path.join(COMP_PCP_TC_IBTRACS_DIR, str(year2), filename)
        infiles.append(infile)


    print infiles

    outfile = "ibtracs.comp_pcp_tc.%d-%d.nc" % (year1, year2)
    outpath = os.path.join(COMP_PCP_TC_IBTRACS_DIR, outfile)

    pl.add_files(infiles, outpath, cube_count=12)

for year in [2016, 2017, 2018]: #2006, 2007, 2008, 2009,2010, 2011, 2012, 2013, 2014, 2015,
    if year == 2016:
        for month in range(7,13):
            composite_pcp_tc_ibtracs_year_month(year, month)
    elif year == 2018:
        for month in range(1,7):
            composite_pcp_tc_ibtracs_year_month(year, month)
    else:
        for month in range(1,13):
            composite_pcp_tc_ibtracs_year_month(year, month)


#year1s = [2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015]
#year2s = [ 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016]
year1s = [2016, 2017]
year2s = [2017, 2018]

for year1, year2 in zip(year1s, year2s):
    composite_pcp_tc_year(year1, year2, lead=None)



