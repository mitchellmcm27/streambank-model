import os
import arcpy
import csv
import json
import mcmtools
import datetime
from dateutil.relativedelta import relativedelta

def calc_precip(start_date, end_date, basin_id_site, site_id): # called once for each site
    """
    Calculate precipitation statistics over a given area throughout given dates.

    Arguments:
    -input: rasters_sorted - a list created by the above function (all rasters)
    -input: start_date, end_date - dates created by the above function
    -input: basin_id_site - BasinID from study site shapefile, will be matched
    -input: site_id - unique identifier of study site, e.g. PR01Y1

    -output: dct_site - a python dictionary with the calculations for input site id
    -output: fields - python list of fields created to store calculations
    """
    print 'Starting date: '
    print start_date
    print 'Ending date:'
    print end_date
    print 'BasinID from site:'
    print basin_id_site

    rasters_month = []
    dct_site = {}
    # Select the drainage basin matching BasinID from the site
    cur_basins = arcpy.da.SearchCursor(basins_shp, 'BasinID')
    basins = [row[0] for row in cur_basins]

    # create a feature layer of current basin only
    layer = 'a_basin'
    where = "{0} = '{1}'".format(arcpy.AddFieldDelimiters(basins_shp, 'BasinID'), basin_id_site)
    arcpy.MakeFeatureLayer_management(basins_shp, layer, where)


    arcpy.env.workspace = 'K:/GIS/MODEL/input/PRISM/precip/monthly'
    thedate = start_date
    while thedate != end_date:
        y = thedate.year
        m = thedate.month
        raster = mcmtools.get_monthly_prism_raster(y, m)
        print thedate
        print 'Processing '+ raster
        monthly = arcpy.Raster(raster)

        ## Do the calculations for the month
        ##
        # Output file setup
        date_code = format(y,'04') + format(m,'02')
        clipped = 'K:/GIS/MODEL/output/temp/' + date_code + '_clip.tif'
        resampled = 'K:/GIS/MODEL/output/temp/' + date_code + '_res.tif'
        table = 'K:/GIS/MODEL/output/tables/' + site_id

        ## Monthly stats
        try:
            # if I already have the needed data, save a little time
            hires = arcpy.Raster(resampled)
        except:
            # this is the first time I'm doing these calculations
            print 'Clipping...'
            arcpy.Clip_management(monthly,clip_bbox,clipped, '#', '#', 'NONE')
            print 'Resampling...'
            arcpy.Resample_management(clipped, resampled, res_cellsz, resampling_method)
            hires = arcpy.Raster(resampled)


        # Create an arcpy table for this site,
        # and add the results to a dict that will eventually
        # have all sites' calculations
        arcpy.sa.ZonalStatisticsAsTable(layer, 'BasinID', hires, table, 'NODATA', 'ALL')
        with arcpy.da.SearchCursor(table, lst_flds) as cur:
            for row in cur:
                k = 0
                for fld in lst_flds:
                    field = fld + date_code # e.g. 'SUM201401' for january 2014
                    if field not in fields:
                        fields.append(field)
                    key = (site_id, field) # e.g. ('PR01Y1', 'SUM201401')
                    dct_site[key] = row[k]
                    print dct_site[key]
                    k += 1
        del rasters_month[:]
        thedate = thedate + relativedelta(months=+1)
    del cur_basins
    return (dct_site, fields) # dictionary, list

"""
Script to match up study sites, drainage basins, and PRISM precipitation data,
calculating statistics on the precip data, and outputting the results to
a shapefile, rasters, and a csv.

As an example, the total and mean precipitation during a given 12-month period
in a given area (drainage basin) can be determined.
"""

## Setup
##
# Required folders
paths_req  = ['K:/GIS/MODEL/output',
                  'K:/GIS/MODEL/output/temp',
                  'K:/GIS/MODEL/output/rasters',
                  'K:/GIS/MODEL/output/tables']
for path in paths_req:
    mcmtools.make_sure_path_exists(path)

# Workspace
arcpy.env.workspace = 'K:/GIS/MODEL/input/PRISM/precip/daily/'
arcpy.env.scratchWorkspace = 'K:/GIS/MODEL/output'
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('Spatial')

# Shapefiles
in_shp = 'K:/GIS/MODEL/input/sites/matlab_out_17-Apr-2016.shp'
out_shp = 'K:/GIS/MODEL/output/python_out.shp'

# Raster processing parameters (used in calcPrecip module)
basins_shp = 'K:/GIS/MODEL/input/basins/04172016.shp'
res_cellsz = '0.0002778' # .0002778 deg. = 1 arc-second ~ 30m
lst_flds = ['MEAN']
resampling_method = 'CUBIC' # NEAREST, BILINEAR, CUBIC, MAJORITY
clip_bbox = '-88 29 -84 32'
dct_results = {} # will hold all final results
fields = [] # will contain a list of fields eventually created by Zonal Statistics

# Copy input sites shapefile to new shapefile, which we will use throughout
arcpy.CopyFeatures_management(in_shp, out_shp)

# Let's get those rasters...
raster_list = arcpy.ListRasters('PRISM*')
# ... and output some useful properties
raster_example = arcpy.Raster(raster_list[0])
px_type = raster_example.pixelType
sref = raster_example.spatialReference
fmt = raster_example.format
print 'Pixel type: ' + px_type
print 'Spatial Reference: ' + sref.name
print 'Raster format: ' + fmt

# Sort the rasters so we can use 'em

## Loop through points in 'sites' shapefile
##
with arcpy.da.SearchCursor(out_shp, ['BasinID', 'date_y1','date_y2', 'id']) as cur:
    i = 0
    for row in cur:
        try:
            # See if the site has 2 valid dates
            sdate, edate = mcmtools.dates_from_site(row[1], row[2])
            basin_id = int(row[0])
            site_id = row[3]
        except:
            # No dates? No bueno. Go to next site.
            continue
        print 'Processing site ' + site_id

        # Let's crunch the numbers
        # Return little dict of site's results
        dct_site, fields = calc_precip(sdate, edate, basin_id, site_id)
        dct_results.update(dct_site) # Append calculation results to big dict
        i += 1
        if i == 3:
            pass # in case you wanna check on some stuff

## Write the returned dict to the site shapefile table and a CSV
##
for field in fields:
    arcpy.AddField_management(out_shp, field, 'DOUBLE')
fields.insert(0, 'id') # <--out_shp already had this, ^no need to add it again^
with arcpy.da.UpdateCursor(out_shp, fields) as cur:
    for row in cur:
        i = 0
        for field in fields:
            key = (row[0], field)
            if key in dct_results:
                row[i] = dct_results[key]
            i += 1
        cur.updateRow(row)

## output shapefile data to a CSV
##
mcmtools.dict2csv(dct_results, 'K:/GIS/MODEL/output/results.csv')


## Merge yearly rasters into 2 raster layers: SUM and MEAN
##
##arcpy.env.workspace = 'K:/GIS/MODEL/output/temp/'
##arcpy.env.compression = 'NONE'
##sum_rasters = arcpy.ListRasters('*year_SUM*')
##print 'There are ' + str(len(sum_rasters)) + ' rasters of SUM'
##avg_rasters = arcpy.ListRasters('*year_AVG*')
##print 'There are ' + str(len(avg_rasters)) + ' rasters of MEAN'
##example = arcpy.Raster(sum_rasters[0])
##sref = example.spatialReference
##cellsz = str(example.meanCellHeight)
##sum_list = ";".join(sum_rasters)
##avg_list = ";".join(avg_rasters)
##try:
##    arcpy.MosaicToNewRaster_management(sum_list, 'K:/GIS/MODEL/output/rasters/',\
##        'SUM.tif', sref, '32_BIT_FLOAT', cellsz, '1')
##    arcpy.MosaicToNewRaster_management(avg_list, 'K:/GIS/MODEL/output/rasters/',\
##        'MEAN.tif',sref, '32_BIT_FLOAT', cellsz, '1')
##except:
##    print arcpy.GetMessages()
##arcpy.CheckInExtension('Spatial')
print 'Finished!'





