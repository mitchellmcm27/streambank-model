import os
import arcpy
import numpy
import csv
import json
import mcmtools
import datetime
import utm
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import ogr
import gdal
import gdal_zonal_stats

def zonal_stats_NLDAS_timeseries(sdate,edate,feat,band):
    """
    Given a month and location, get values from NLDAS raster
    """
    date = sdate
    ts = []
    dates = []

    while date < edate:
        datecode = ('%04d%02d' % (date.year, date.month))
        rpath = nldas_path + '/NLDAS_VIC0125_M.A'+datecode+'.002.grb'
        print rpath
        # this allows GDAL to throw Python Exceptions
        gdal.UseExceptions()

        try:
            src_ds = gdal.Open(rpath)
        except RuntimeError, e:
            print 'Unable to open INPUT.tif'
            print e
            sys.exit(1)

        # Save band as temporary geotiff
        #Open output format driver, see gdal_translate --formats for list
        fmt = "GTiff"
        driver = gdal.GetDriverByName( fmt )

        #Output to new format
        dst_filename = 'K:/GIS/MODEL/output/temp/bandcopy.tif'
        os.system("gdal_translate -of GTiff -b %d " %band + rpath + " " +  dst_filename)
        #Properly close the datasets to flush to disk
        src_ds = None

        lowres = arcpy.Raster(dst_filename)
        print lowres.extent
        break
    return
##        ts.append(avg)
##        dates.append(wildcard)
##        date = date + relativedelta(months=+1)
##    return numpy.asarray(ts), dates

## Setup
##
nldas_path = 'K:/GIS/MODEL/input/NLDAS/NLDAS_VIC0125_M.002'

sites_shp_path = 'K:/GIS/MODEL/input/sites/matlab_out_17-Apr-2016.shp'
csv_path = 'K:/GIS/MODEL/output/NLDAS_VIC.csv'
basins_shp_path = 'K:/GIS/Model/input/basins/04172016.shp'

## Loop through points in 'sites' shapefile
##
with open(csv_path,'w+b') as f:
        writer = csv.writer(f)
        writer.writerow(['id','dates','runoff','baseflow','soil_moisture','lai'])

basins_shp = ogr.Open(basins_shp_path)
basins_lyr = basins_shp.GetLayer()
basins_featList = range(basins_lyr.GetFeatureCount())

sites_shp = ogr.Open(sites_shp_path)
sites_lyr = sites_shp.GetLayer()
sites_featList = range(sites_lyr.GetFeatureCount())

for sFID in sites_featList:
    site = sites_lyr.GetFeature(sFID)
    basinid = int(site.GetField('BasinID'))
    site_id = site.GetField('id')
    d1 = site.GetField('date_y1')
    d2 = site.GetField('date_y2')

    if d1=='NaT':
        continue

    sdate, edate = mcmtools.dates_from_site(d1, d2) # one month before to 12 months after starting date

    # find feature from Basin shapefile
    for bFID in basins_featList:
        basin = basins_lyr.GetFeature(bFID)
        thisid = int(basin.GetField('BasinID'))
        if thisid is basinid:
            break
    print basinid
    [runoff_timeseries, dates] = zonal_stats_NLDAS_timeseries(sdate,edate,basin,12) # kg/m2

    break
##            [baseflow_timeseries, dates] = sample_NLDAS_timeseries(sdate,edate,basin,13) # kg/m2
##            [soil_moisture_timeseries, dates] = sample_NLDAS_timeseries(sdate,edate,basin,25) # kg/m2
##            [lai_timeseries, dates] = sample_NLDAS_timeseries(sdate,edate,basin,50) # 0-9
##
##            writer.writerow([site_id,list(dates),list(runoff_timeseries),list(baseflow_timeseries),list(soil_moisture_timeseries),list(lai_timeseries)])
##
##            i += 1
##            if i == 3:
##                pass # in case you wanna check on some stuff