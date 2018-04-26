from dateutil.relativedelta import relativedelta
from osgeo import gdal,ogr
import struct

def point_value(ras, lat, lon):
    src_ds=gdal.Open(arcpy.env.workspace+'/'+ras)
    gt=src_ds.GetGeoTransform()
    rb=src_ds.GetRasterBand(1)
    px = int((lon - gt[0]) / gt[1]) #x pixel
    py = int((lat - gt[3]) / gt[5]) #y pixel
    pval = rb.ReadAsArray(px,py,1,1) #np array
    return float(pval[0])

def calc_PRISM_daily( sdate, edate, site_id, lat, lon): # called for each site
    """
    Calculate average monthly precip over a given area
    """
    print 'Starting date: '
    print sdate
    arcpy.env.workspace = base_workspace + 'PRISM/precip/daily'
    # Select the drainage basin matching BasinID from the site
    P_daily = []
    dates = []
    date=sdate
    j=0
    while date < edate:
        wildcard = '%04d%02d%02d' % (date.year, date.month, date.day)
        if date.day==1:
            print wildcard
        rlist = arcpy.ListRasters('*_'+wildcard+'_*')
        ras = rlist[0]
        P_day = point_value(ras, lat, lon)
        P_daily.append(P_day)
        dates.append(int(wildcard))
        date = date + relativedelta(days=+1)
        j+=1

    print '\n'
    arcpy.env.workspace=base_workspace
    return (numpy.asarray(P_daily),numpy.asarray(dates)) # NP array

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
              'K:/GIS/MODEL/output/tables',
              'K:/GIS/MODEL/output/daily_rainfall']
for path in req_paths:
    make_sure_path_exists(path)

# Workspace
base_workspace = 'K:/GIS/MODEL/input/'
arcpy.env.workspace = base_workspace
arcpy.env.scratchWorkspace = 'K:/GIS/MODEL/output'
arcpy.env.overwriteOutput = True
#arcpy.CheckOutExtension('Spatial')

# Input and output
in_shp = 'K:/GIS/MODEL/input/sites/matlab_out.shp'
out_csv = 'K:/GIS/MODEL/output/daily_rainfall.csv'

headers = {'id':'id',
            'BasinID':'BasinID',
            'P_daily':'P_daily'}
write_dict2csv(headers, out_csv) # Create/overwrite csv and write headers
count = arcpy.GetCount_management(in_shp)
nrows = int(count.getOutput(0))
ncols = len(headers)

with arcpy.da.SearchCursor(in_shp, ['BasinID', 'id', 'date_y1', 'lat', 'long','date_y3']) as cur:
    dct_site = headers.fromkeys(headers, 'NaN')

    # Numpy object array to hold results. First row is headers
    results = numpy.zeros((nrows+1,ncols), dtype=numpy.object)
    results[0,:] = numpy.array([
        'id', 'BasinID', 'P_daily'],
        dtype=numpy.object)


    for i, row in enumerate(cur, start=1):

        # Site must have 5 attributes, otherwise something won't work
        try:
            sdate = mcmtools.translate_date(row[2]) # to 12 months after starting date
            sdate = datetime.date(2013,01,01)
            edate = mcmtools.translate_date(row[5]) # edate is exclusive
            site_id = row[1]
            lat = row[3]
            lon = row[4]
            basin_id = int(row[0])

        except:
            print 'abc'
            continue

        print '\n---------------------'
        print '  Processing ' + site_id
        print '-----------------------'

        # Return a numpy vector of the site results
        P_daily,dates = calc_PRISM_daily(sdate,edate,site_id, lat, lon)
        fpath = 'K:/GIS/MODEL/output/daily_rainfall/'+site_id+'_'+row[2]+'_'+row[5]+'.txt'
        print fpath
        numpy.savetxt(fpath,numpy.c_[P_daily,dates],fmt=['%.2f','%08d'])

        dct_site['id'] = site_id
        results[i,0] = site_id

        dct_site['BasinID'] = basin_id
        results[i,1] = basin_id

        dct_site ['P_daily'] = P_daily
        results[i,2] = P_daily

#--------------------------------------------------------------
# Change this if you don't want to overwrite an older mat file!
#
output_type = 'matlab' # <-------------------------------------
#--------------------------------------------------------------
if output_type == 'matlab':
    scipy.io.savemat('K:/GIS/MODEL/output/daily_rainfall.mat', {'results':results})



