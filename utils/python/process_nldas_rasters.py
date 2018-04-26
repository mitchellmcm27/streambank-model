import os
import arcpy
import numpy
import csv
import json
import mcmtools
import datetime
import utm
from dateutil.relativedelta import relativedelta

def sample_NLDAS_timeseries(sdate,edate,lonlat_tuple,band):
    """
    Given a month and location, get values from NLDAS raster
    """
    xy_tuple = utm.from_latlon(lonlat_tuple[1],lonlat_tuple[0])
    date = sdate
    ts = []
    while date < edate:
        wildcard = ('*%04d%02d*' % (date.year, date.month))
        date = date + relativedelta(months=+1)
        raster_list = arcpy.ListRasters(wildcard)
        raster = arcpy.Raster(raster_list[0]) # should only be one in the list

        result = arcpy.GetCellValue_management(raster,str(lonlat_tuple[0])+" "+str(lonlat_tuple[1]),str(band))
        val = float(result.getOutput(0))
        ts.append(val)
    return numpy.asarray(ts)

## Setup
##
arcpy.env.workspace = 'K:/GIS/r/NLDAS/NLDAS_NOAH0125_M.002'
arcpy.env.scratchWorkspace = 'K:/GIS/MODEL/output'
arcpy.env.overwriteOutput = True
in_shp = 'K:/GIS/MODEL/input/sites/matlab_out_17-Apr-2016.shp'
csv_path = 'K:/GIS/MODEL/output/NLDAS_NOAH.csv'
arcpy.CheckOutExtension('Spatial')
INT_NAN =  -99999

## Copy input sites shapefile to new shapefile, which we will use throughout
##
out_shp = 'K:/GIS/MODEL/input/sites/with_NLDAS.shp'
arcpy.CopyFeatures_management(in_shp, out_shp)

## Loop through points in 'sites' shapefile
##
with open(csv_path,'w+b') as f:
        writer = csv.writer(f)
        writer.writerow(['id','runoff','baseflow','soil_moisture','lai'])
with open(csv_path,'a') as f:
    writer = csv.writer(f)
    with arcpy.da.SearchCursor(out_shp, ['id','date_y1','date_y2', 'SHAPE@XY']) as cur:
        i = 0
        for row in cur:
            site_id = row[0]
            try:
                # See if the site has 2 valid dates
                sdate, edate = mcmtools.dates_from_site(row[1], row[2]) # one month before to 12 months after starting date
            except:
                print "No date or dates not valid"
                writer.writerow([site_id,'[]','[]','[]','[]'])
                continue
            print 'Processing site ' + site_id

            # Let's crunch the numbers
            print sdate
            print edate
            runoff_timeseries = sample_NLDAS_timeseries(sdate,edate,row[3],12) # kg/m2
            baseflow_timeseries = sample_NLDAS_timeseries(sdate,edate,row[3],13) # kg/m2
            soil_moisture_timeseries = sample_NLDAS_timeseries(sdate,edate,row[3],25) # kg/m2
            lai_timeseries = sample_NLDAS_timeseries(sdate,edate,row[3],50) # 0-9

            writer.writerow([site_id,list(runoff_timeseries),list(baseflow_timeseries),list(soil_moisture_timeseries),list(lai_timeseries)])

            i += 1
            if i == 3:
                pass # in case you wanna check on some stuff