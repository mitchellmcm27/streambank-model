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
import gdal
import ogr
import sys


def zonal_NLDAS_timeseries(sdate,edate,site_id,band,a_basin):
    """
    Given a month and location, get values from NLDAS raster
    """
    date = sdate
    ts = []
    dates = []
    oldws = arcpy.env.workspace
    arcpy.env.workspace = 'K:/GIS/MODEL/output/temp/'
    while date < edate:
        datecode = ('%04d%02d' % (date.year, date.month))
        if datecode =='201605':
            datecode = '201505' # use last year's data until May 2015 is released
        # Output file setup
        base_date_band_code = base[6] + '' + datecode + '_' + ('%02d' % (band))
        clipped = base_date_band_code + '_clip.tif'
        resampled = base_date_band_code + '_res.tif'

        try:
            # if I already have the needed data, save a little time
            resampled = arcpy.Raster(resampled)
        except:
            # this is the first time I'm doing these calculations
            rpath = nldas_path + base + datecode + '.002.grb'
            print rpath

            #Output to new format
            dst_filename = arcpy.env.workspace + datecode + '.tif'
            os.system("gdal_translate -of GTiff -b %d " %band + rpath + " " +  dst_filename)

            lowres = arcpy.Raster(dst_filename)
            print 'Clipping...'
            arcpy.Clip_management(lowres,clip_bbox,clipped, '#', '#', 'NONE')
            print 'Resampling...'
            arcpy.Resample_management(clipped, resampled, res_cellsz, resampling_method)

        print 'Zonal statistics...'
        table = 'in_memory/'+site_id
        arcpy.sa.ZonalStatisticsAsTable(a_basin, "BasinID", resampled, table, "NODATA", "MEAN")

        # Get value from table
        with arcpy.da.SearchCursor(table,['MEAN']) as tcur:
            for val in tcur:
                ts.append(float(val[0]))
                dates.append(int(datecode))
        date = date + relativedelta(months=+1)
        table = None
    print ts
    print dates
    arcpy.env.workspace=oldws
    return numpy.asarray(ts), numpy.asarray(dates)

## Setup
##

############## Choose Model ####################
model = 'NOAH'

nldas_path = 'K:/GIS/MODEL/input/NLDAS/NLDAS_'+model+'0125_M.002/'
arcpy.env.workspace=nldas_path
base = 'NLDAS_'+model+'0125_M.A'
arcpy.env.scratchWorkspace = 'K:/GIS/MODEL/output/'
arcpy.env.overwriteOutput = True
in_shp = 'K:/GIS/MODEL/input/sites/matlab_out_17-May-2016.shp'
csv_path = 'K:/GIS/MODEL/output/NLDAS_'+model+'.csv'
basins_shp = 'K:/GIS/Model/input/basins/04172016.shp'
arcpy.CheckOutExtension('Spatial')
INT_NAN =  -99999
clip_bbox = '-88 29 -84 32'
res_cellsz = '0.002778' # .0002778 deg. = 1 arc-second ~ 30m
resampling_method = 'CUBIC' # NEAREST, BILINEAR, CUBIC, MAJORITY

nldas_lut = {('VIC','baseflow'):13,
             ('VIC','runoff'):12,
             ('VIC','soil moisture'):24,
             ('VIC','LAI'):40,
             ('NOAH','baseflow'):13,
             ('NOAH','runoff'):12,
             ('NOAH','soil moisture'):23,
             ('NOAH','LAI'):50}


## Copy input sites shapefile to new shapefile, which we will use throughout
##


## Loop through points in 'sites' shapefile
##
with open(csv_path,'w+b') as f:
        writer = csv.writer(f)
        writer.writerow(['id','BasinID','dates','runoff','baseflow','soil_moisture','lai'])
with open(csv_path,'a') as f:
    writer = csv.writer(f)
    with arcpy.da.SearchCursor(in_shp, ['id','date_y1','date_y2', 'date_y3','BasinID']) as cur:
        i = 0
        for row in cur:
            site_id = row[0]
            try:
                # See if the site has 2 valid dates
                sdate, sdate12 = mcmtools.dates_from_site(row[1]) # to 12 months after starting date
                if row[3]=='NaT':
                    sdate = datetime.date(2016,5,1)
                    edate = datetime.date(2016,6,1)
                else:
                    continue
                    edate, edate12 = mcmtools.dates_from_site(row[3])
            except:
                print "No date or dates not valid"
                writer.writerow([site_id,int(row[4]),'[]','[]','[]','[]'])
                continue
            print 'Processing site ' + site_id

            # Let's crunch the numbers
            print sdate
            print edate
            # Select the drainage basin matching BasinID from the site
            with arcpy.da.SearchCursor(basins_shp, ['BasinID','SHAPE@']) as bcur:
                for basin in bcur:
                    if float(basin[0])==float(row[4]):
                        print basin[0]
                        print row[4]
                        a_basin = 'a_basin'
                        where = "{0} = '{1}'".format(arcpy.AddFieldDelimiters(basins_shp, 'BasinID'), basin[0])
                        arcpy.MakeFeatureLayer_management(basins_shp, a_basin, where)
                        break

            # sample_xxx_timeseries(starting date, ending date, unique id, band, basin)
            site_id = row[0]
            basin_id = int(row[4])
            runoff_timeseries, dates = zonal_NLDAS_timeseries(sdate,edate,site_id,nldas_lut[(model,'runoff')],a_basin) # kg/m2
            baseflow_timeseries, dates = zonal_NLDAS_timeseries(sdate,edate,site_id,nldas_lut[(model,'baseflow')],a_basin) # kg/m2
            soil_moisture_timeseries, dates = zonal_NLDAS_timeseries(sdate,edate,site_id,nldas_lut[(model,'soil moisture')],a_basin) # kg/m2
            lai_timeseries, dates = zonal_NLDAS_timeseries(sdate,edate,site_id,nldas_lut[(model,'LAI')],a_basin) # 0-9

            writer.writerow([site_id,basin_id,list(dates),list(runoff_timeseries),list(baseflow_timeseries),list(soil_moisture_timeseries),list(lai_timeseries)])