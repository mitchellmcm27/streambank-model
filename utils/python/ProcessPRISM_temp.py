from dateutil.relativedelta import relativedelta

def basin_average(ras, basin):
    '''
    Average a raster value within a drainage basin. Clip to bounding box, resample, clip to geometry
    '''
    extent = arcpy.Describe(basin).extent
    bbox = [extent.XMin-0.1, extent.YMin-0.1, extent.XMax+0.1, extent.YMax+0.1]
    bbox = ' '.join(map(str, bbox))
    clipped = 'K:/GIS/MODEL/output/temp/clipped.tif'
    arcpy.Clip_management(ras,bbox,clipped)

    resampled = 'K:/GIS/MODEL/output/temp/resampled.tif'
    arcpy.Resample_management(clipped, resampled, res_cellsz, resampling_method)

    clipped_basin = 'K:/GIS/MODEL/output/temp/clipped_basin.tif'
    arcpy.Clip_management(resampled, '#', clipped_basin, basin, '#', 'ClippingGeometry')
    x = arcpy.RasterToNumPyArray(clipped_basin, '#', '#', '#', -9999)
    mx = numpy.ma.masked_values(x, -9999.)
    return mx.mean()

def point_value(ras, lat, lon):
    point = str(lon)+' '+str(lat)
    return arcpy.GetCellValue_management(ras,point)

def calc_PRISM_monthly( sdate, edate, basin_id_site, site_id, lat, lon): # called for each site
    """
    Calculate average monthly precip over a given area
    """
    print 'Starting date: '
    print sdate
    print 'BasinID from site:'
    print basin_id_site

    # Select the drainage basin matching BasinID from the site
    cBasins = arcpy.da.SearchCursor(basins_shp, 'BasinID')
    Basins = [row[0] for row in cBasins]
    for basin in Basins:
        basin_id_basins = int(basin)
        if basin_id_basins == basin_id_site:
            print 'Matched basin ' + str(basin_id_basins) + ' to ' + str(site_id)+'\n'
            # create a feature layer of current basin only
            layer = 'a_basin'
            where = "{0} = '{1}'".format(arcpy.AddFieldDelimiters(basins_shp, 'BasinID'), basin_id_basins)
            arcpy.MakeFeatureLayer_management(basins_shp, layer, where)
            this_basin = 'K:/GIS/MODEL/output/temp/this_basin.shp'
            arcpy.CopyFeatures_management(layer, this_basin)
            del cBasins
            break

    P = []
    Tmean = []
    Tmax = []
    Tmin = []
    f = []
    alpha = []

    date=sdate

    j=0
    while date < edate:
        wildcard = '%04d%02d' % (date.year, date.month)

        print wildcard
##
        ## Monthly average precip
        arcpy.env.workspace = base_workspace + 'PRISM/precip/monthly'
        rlist = arcpy.ListRasters('*_'+wildcard+'*')
        ras = rlist[0]
        P.append(basin_average(ras, this_basin))
        print str(P[j]) + ' mm average precip'

##        ## Monthly average mean temp
##        arcpy.env.workspace = base_workspace + 'PRISM/tmean'
##        rlist = arcpy.ListRasters('*_'+wildcard+'*')
##        ras = rlist[0]
##        Tmean.append(basin_average(ras, this_basin))
##        print str(Tmean[j]) + ' C mean daily average temp'
##
##        ## Monthly average maximum temp
##        arcpy.env.workspace = base_workspace + 'PRISM/tmax'
##        rlist = arcpy.ListRasters('*_'+wildcard+'*')
##        ras = rlist[0]
##        Tmax.append(basin_average(ras, this_basin))
##        print str(Tmax[j]) + ' C mean daily max temp'
##
##        ## Monthly average minimum temp
##        arcpy.env.workspace = base_workspace + 'PRISM/tmin'
##        rlist = arcpy.ListRasters('*_'+wildcard+'*')
##        ras = rlist[0]
##        Tmin.append(basin_average(ras, this_basin))
##        print str(Tmin[j]) + ' C mean daily min temp'

        ## Daily wet or dry days
        arcpy.env.workspace = base_workspace + 'PRISM/precip/daily'
        rlist = arcpy.ListRasters('*_'+wildcard+'*')
        wet = 0
        dry = 0
        for ras in rlist:
            result = point_value(ras, lat, lon)
            precip = float(result.getOutput(0))
            if precip < 1:
                dry += 1
            else:
                wet += 1
        print str(wet) + ' wet days'
        print str(dry) + ' dry days'
        f.append(float(wet)/(float(wet)+float(dry)))
        print str(f[j]) + ' storm frequency'
##        alpha.append(P[j]/float(wet))
##        print str(alpha[j]) + ' mm average precip per wet day\n'
        j+=1
        date = date + relativedelta(months=+1)
    print '\n'
    return (numpy.asarray(P), numpy.asarray(Tmean), numpy.asarray(Tmax), numpy.asarray(Tmin), numpy.asarray(f), numpy.asarray(alpha)) # NP array

def calc_PET(Tmean_monthly, Tmax_monthly, Tmin_monthly, lat, sdate,edate):
    """
    Calulate potential evapotranspiration from temperatures and latitude and dates
    """
    from numpy import sin, cos, tan, arccos, pi, sqrt
    from calendar import monthrange


    mos = []
    yrs = []
    days = []
    date = sdate
    while date < edate:
        mos.append(date.month)
        yrs.append(date.year)
        days.append(monthrange(date.year,date.month)[1])
        date=date+relativedelta(months=+1)
    mos = numpy.asarray(mos)
    yrs = numpy.asarray(yrs)
    print mos
    print yrs
    print days

    lat = numpy.radians(lat)
    JD = numpy.asarray(mos*30 - 15, dtype=float) # days of the year at the 15th of each month
    TD = Tmax_monthly - Tmin_monthly

    KT = 0.00185*TD*TD - 0.0433*TD + 0.4023

    d = 0.409*sin(2*pi*JD/365 - 1.39) # solar declination (rad)
    S = arccos(-tan(lat)*tan(d)) # sunset hour angle (rad)
    G = 0.0820 # Mj/m^2/min
    r = 1 + 0.033*cos(2*pi*JD/365) # inverse rel. dist. from sun
    Ra = 1440/pi * (G*r) * (S*sin(lat)*sin(d) + cos(lat)*cos(d)*sin(S)) # Mj/m^2/d
    Ra = Ra/2.43 # mm/day

    PET = 0.0135*KT*Ra*sqrt(TD)*(Tmean_monthly + 17.8) # mm/day
    PET = PET*numpy.asarray(days, dtype=float) # mm/month

    return PET

## Modules
##

import arcpy
import numpy
import scipy.io
import matplotlib.pyplot as plt
import mcmtools
from mcmtools import *
import datetime
## Setup
##
# Required folders
req_paths  = ['K:/GIS/MODEL/output',
              'K:/GIS/MODEL/output/temp',
              'K:/GIS/MODEL/output/rasters',
              'K:/GIS/MODEL/output/tables']
for path in req_paths:
    make_sure_path_exists(path)

# Workspace
base_workspace = 'K:/GIS/MODEL/input/'
arcpy.env.workspace = base_workspace
arcpy.env.scratchWorkspace = 'K:/GIS/MODEL/output'
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('Spatial')

# Input and output shapefiles
in_shp = 'K:/GIS/MODEL/input/sites/matlab_out.shp'
out_shp = 'K:/GIS/MODEL/output/python_out_precip.shp'
out_csv = 'K:/GIS/MODEL/output/precip_results.csv'

# Raster processing parameters (used in calcPrecip module)
basins_shp = 'K:/GIS/MODEL/input/basins/04172016.shp'
res_cellsz = '0.000833333' # deg = 3 arc-second ~ 90m
resampling_method = 'CUBIC' # NEAREST, BILINEAR, CUBIC, MAJORITY

# Copy input sites shapefile to new shapefile
arcpy.CopyFeatures_management(in_shp, out_shp)

headers = {'id':'id',
            'BasinID':'BasinID',
            'lat':'lat',
            'P_monthly':'P_monthly',
            'Tmean_monthly':'Tmean_monthly',
            'Tmax_monthly':'Tmax_monthly',
            'Tmin_monthly':'Tmin_monthly',
            'PET_monthly':'PET_monthly',
            'R_monthly':'R_monthly',
            'f_monthly':'f_monthly',
            'alpha_monthly':'alpha_monthly'}
write_dict2csv(headers, out_csv) # Create/overwrite csv and write headers
count = arcpy.GetCount_management(in_shp)
nrows = int(count.getOutput(0))
ncols = len(headers)

with arcpy.da.SearchCursor(in_shp, ['BasinID', 'id', 'date_y1', 'lat', 'long','date_y3']) as cur:
    dct_site = headers.fromkeys(headers, 'NaN')

    # Numpy object array to hold results. First row is headers
    results = numpy.zeros((nrows+1,ncols), dtype=numpy.object)
    results[0,:] = numpy.array([
        'id', 'BasinID', 'lat', 'P_monthly',
        'Tmean_monthly', 'Tmax_monthly', 'Tmin_monthly',
        'PET_monthly', 'R_monthly', 'f_monthly','alpha_monthly'],
        dtype=numpy.object)


    for i, row in enumerate(cur, start=1):

        # Site must have 5 attributes, otherwise something won't work
        try:
            sdate, sdate12 = mcmtools.dates_from_site(row[2]) # to 12 months after starting date
            edate, edate12 = mcmtools.dates_from_site(row[5]) # edate is exclusive

            basin_id = int(row[0])
            site_id = row[1]
            lat = row[3]
            lon = row[4]
        except:
            continue

        print '\n---------------------'
        print '  Processing ' + site_id
        print '-----------------------'

        # Return a numpy vector of the site results
        P_monthly, Tmean_monthly, Tmax_monthly, Tmin_monthly, f_monthly, alpha_monthly =\
            calc_PRISM_monthly(sdate,edate,basin_id,site_id, lat, lon)
##
##        # Calculate Potential Evapotranspiration
##        PET_monthly = calc_PET(Tmean_monthly, Tmax_monthly,
##                               Tmin_monthly, lat, sdate,edate)
##
##        R_monthly = P_monthly - PET_monthly

        # Write results to dictionary for CSV
        # Write results to numpy object array

        dct_site['id'] = site_id
        results[i,0] = site_id

        dct_site['BasinID'] = basin_id
        results[i,1] = basin_id

##        dct_site['lat'] = lat
##        results[i,2] = lat
##
        dct_site ['P_monthly'] = P_monthly
        results[i,3] = P_monthly
##
##        dct_site['Tmean_monthly'] = Tmean_monthly
##        results[i,4] = Tmean_monthly
##
##        dct_site['Tmax_monthly'] = Tmax_monthly
##        results[i,5] = Tmax_monthly
##
##        dct_site['Tmin_monthly'] = Tmin_monthly
##        results[i,6] = Tmin_monthly
##
##        dct_site['PET_monthly'] = PET_monthly
##        results[i,7] = PET_monthly
##
##        dct_site['R_monthly'] = R_monthly
##        results[i,8] = R_monthly

        dct_site['f_monthly'] = f_monthly
        results[i,9] = f_monthly
##
##        dct_site['alpha_monthly'] = alpha_monthly
##        results[i,10] = alpha_monthly

        # Append calculation results to the csv as backup to mat file
        append_dict2csv(dct_site,'K:/GIS/MODEL/output/results.csv')

        if i == 0:
            pass # in case you wanna check on some stuff


#--------------------------------------------------------------
# Change this if you don't want to overwrite an older mat file!
#
output_type = 'matlab' # <-------------------------------------
#--------------------------------------------------------------
if output_type == 'matlab':
    scipy.io.savemat('K:/GIS/MODEL/output/precip_final.mat', {'results':results})



