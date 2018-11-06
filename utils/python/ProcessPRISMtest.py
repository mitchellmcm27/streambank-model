## Modules
##
import os
import arcpy
import csv
import json
import mcmtools

## Main Program
"""
Module to match up study sites, drainage basins, and PRISM precipitation data,
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
in_shp = 'K:/GIS/MODEL/input/sites/matlab_out.shp'
out_shp = 'K:/GIS/MODEL/output/python_out.shp'

# Raster processing parameters (used in calcPrecip module)
basins_shp = 'K:/GIS/MODEL/input/basins/Aug6.shp'
res_cellsz = '0.0002778' # .0002778 deg. = 1 arc-second ~ 30m
lst_flds = ['MEAN', 'SUM']
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
rasters_sorted = mcmtools.sort_prism_rasters(raster_list)

## Loop through points in 'sites' shapefile
##
with arcpy.da.SearchCursor(out_shp, ['BasinID', 'date_y1','date_y2', 'id']) as cur:
    i = 0
    for row in cur:
        try:
            # See if the site has 2 valid dates
            sdate, edate = dates_from_site(row[1], row[2])
            basin_id = int(row[0])
            site_id = row[3]
        except:
            # No dates? No bueno. Go to next site.
            continue
        print 'Processing site ' + site_id

        # Let's crunch the numbers
        # Return little dict of site's results
        dct_site, fields = calc_precip(rasters_sorted, sdate,
                                       edate, basin_id, site_id)
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
mcmtools.shp2csv(out_shp, fields, 'K:/GIS/MODEL/output/results.csv')

## Merge yearly rasters into 2 raster layers: SUM and MEAN
##
arcpy.env.workspace = 'K:/GIS/MODEL/output/temp/'
arcpy.env.compression = 'NONE'
sum_rasters = arcpy.ListRasters('*year_SUM*')
print 'There are ' + str(len(sum_rasters)) + ' rasters of SUM'
avg_rasters = arcpy.ListRasters('*year_AVG*')
print 'There are ' + str(len(avg_rasters)) + ' rasters of MEAN'
example = arcpy.Raster(sum_rasters[0])
sref = example.spatialReference
cellsz = str(example.meanCellHeight)
sum_list = ";".join(sum_rasters)
avg_list = ";".join(avg_rasters)
try:
    arcpy.MosaicToNewRaster_management(sum_list, 'K:/GIS/MODEL/output/rasters/',\
        'SUM.tif', sref, '32_BIT_FLOAT', cellsz, '1')
    arcpy.MosaicToNewRaster_management(avg_list, 'K:/GIS/MODEL/output/rasters/',\
        'MEAN.tif',sref, '32_BIT_FLOAT', cellsz, '1')
except:
    print arcpy.GetMessages()
arcpy.CheckInExtension('Spatial')
print 'Finished!'

def dates_from_site(date1, date2):
    """
    Function to convert 2 formatted datestrings into lists of ints
    -input format must match 'dd MMM yyyy', e.g. '01 Jan 2014'
    -output format list = [yy][mm][dd]
    """
    if ( date1 == 'NaT' ):
        return
    starting_day = int(date1[:2])
    starting_month = mcmtools.convert_month(date1[3:6])
    starting_year = int(date1[7:])
    ending_day = int(date2[:2])
    ending_month = mcmtools.convert_month(date2[3:6])
    ending_year = int(date2[7:])

    # Force ending date to be 1 year after start
    if starting_day > 15: # if more than halfway through the month,
        starting_month += 1 # start at the first of the next month
    if starting_month > 12:
        starting_month = 1 # go to January of next year if needed
        starting_year += 1

    starting_day = 1
    ending_day = 31
    ending_month = starting_month - 1 # Exactly 12 months

    if ending_month == 0: # Special case of starting on January 1st
        ending_month = 12
        ending_year = starting_year
    else:
        ending_year = starting_year + 1


    starting_date = [starting_year, starting_month, starting_day]
    ending_date = [ending_year, ending_month, ending_day]
    return (starting_date, ending_date)
def calc_precip( rasters_sorted, start_date, end_date, basin_id_site, site_id): # called for each site
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

    starting_year = start_date[0] - 2012
    starting_month = start_date[1]
    starting_day = start_date[2]
    ending_year = end_date[0] - 2012
    ending_month = end_date[1]
    ending_day = end_date[2]
    rasters_month = []
    dct_site = {}
    # Select the drainage basin matching BasinID from the site
    cur_basins = arcpy.da.SearchCursor(basins_shp, 'BasinID')
    basins = [row[0] for row in cur_basins]
    for basin in basins:
        basin_id_basins = int(basin)
        if basin_id_basins == basin_id_site:
            print 'Matched basin ' + str(basin_id_basins) + ' to ' + str(site_id)

            # create a feature layer of current basin only
            layer = 'a_basin'
            where = "{0} = '{1}'".format(arcpy.AddFieldDelimiters(basins_shp, 'BasinID'), basin_id_basins)
            arcpy.MakeFeatureLayer_management(basins_shp, layer, where)

            # Select the appropriate rasters
            j = 0
            for y in range(starting_year, ending_year+1):
                if y==starting_year:
                    m1 = starting_month
                else:
                    m1 = 1
                if y==ending_year:
                    m2 = ending_month
                else:
                    m2 = 12
                for m in range(m1, m2+1):

                    # Most of the calculations take place within this month loop
                    if y==starting_year and m==starting_month:
                        d1 = starting_day
                    else:
                        d1 = 1
                    if y==ending_year and m==ending_month:
                        d2 = ending_day
                    else:
                        d2 = 31

                    # Loop through the valid days in the month,
                    # getting all the needed rasters into a list
                    if d1 == 1 and d2 == 31:
                        # Use monthly rasters. They are found in...
                        arcpy.env.workspace =\
                            'K:/GIS/MODEL/input/PRISM/precip/monthly'
                        raster = mcmtools.get_monthly_prism_raster(y+2012, m)
                        sum_month = arcpy.Raster(raster)
                        days = 31
                    else:
                        # Use the daily rasters. They are found in...
                        arcpy.env.workspace = (
                            'K:/GIS/MODEL/input/PRISM/precip/daily')
                        days = 0
                        for d in range(d1, d2 + 1):
                            if rasters_sorted[y][m][d] != '':
                                days += 1
                                rasters_month.append(rasters_sorted[y][m][d])

                        # Still in month loop...
                        # Loop through rasters between dates adding them to a
                        # monthly raster
                        rasters_month = filter(None, rasters_month)
                        i = 0
                        for raster in rasters_month:
                            ras = arcpy.Raster(raster)
                            if i==0:
                                sum_month = ras
                            else:
                                sum_month += ras
                            i += 1
                    if j==0: # let's also get all the rasters from a 12 mo. period
                        sum_all = sum_month
                        print sum_all
                    elif j<12:
                        sum_all += sum_month
                    else:
                        raise Exception('More than 12 months given in a year!')
                    j += 1

                    ## Do the calculations for the month
                    ##
                    # Output file setup
                    date_code = str(y+2012) + '-' + str(m) + '_' + str(days) + 'days'
                    out_sum = 'K:/GIS/MODEL/output/rasters/' \
                                + site_id + '_' + str(basin_id_site) + '_' + date_code + '_SUM.tif'
                    out_avg = 'K:/GIS/MODEL/output/rasters/' \
                                + site_id + '_' + str(basin_id_site) + '_' + date_code + '_MEAN.tif'
                    clipped = 'K:/GIS/MODEL/output/temp/' + date_code + '_clip.tif'
                    resampled = 'K:/GIS/MODEL/output/temp/' + date_code + '_res.tif'
                    table = 'K:/GIS/MODEL/output/tables/' + site_id

                    ## Monthly stats
                    try:
                        # if I already have the needed data, save a little time
                        inraster = arcpy.Raster(resampled)
                    except:
                        # this is the first time I'm doing these calculations
                        print 'Clipping...'
                        arcpy.Clip_management(sum_month,clip_bbox,clipped, '#', '#', 'NONE')
                        print 'Resampling...'
                        arcpy.Resample_management(clipped, resampled, res_cellsz, resampling_method)
                        inraster = arcpy.Raster(resampled)


                    # Create an arcpy table for this site,
                    # and add the results to a dict that will eventually
                    # have all sites' calculations
                    arcpy.sa.ZonalStatisticsAsTable( layer, 'BasinID', inraster, table, 'NODATA', 'ALL')
                    with arcpy.da.SearchCursor(table, lst_flds) as cur:
                        for row in cur:
                            k = 0
                            for fld in lst_flds:
                                field = fld + '_%02d' % m # e.g. 'SUM_01' for january
                                if field not in fields:
                                    fields.append(field)
                                key = (site_id, field) # e.g. ('PR01Y1', 'SUM_01')
                                dct_site[key] = row[k]
                                k += 1
                    del rasters_month[:]

            ## 12-month stats (back in the year loop)
            ##
            resampled_all = 'K:/GIS/MODEL/output/temp/' + site_id + '_1year_res.tif'
            clipped_all = 'K:/GIS/MODEL/output/temp/' + site_id + '_1year_clip.tif'
            try:
                # if I already have the needed data, save a little time
                inraster = arcpy.Raster(resampled_all)
            except:
                # this is the first time I'm doing these calculations
                arcpy.Clip_management(sum_all, clip_bbox, clipped_all, '#', '#', 'NONE')
                arcpy.Resample_management(clipped_all, resampled_all, res_cellsz, resampling_method)
                inraster = arcpy.Raster(resampled_all)
            # Output file setup
            sum_basin_raster = arcpy.sa.ZonalStatistics( layer, 'BasinID', inraster, 'SUM', 'NODATA')
            avg_basin_raster = arcpy.sa.ZonalStatistics( layer, 'BasinID', inraster, 'MEAN', 'NODATA')
            out_sum = 'K:/GIS/MODEL/output/temp/' + site_id + '_year_SUM.tif'
            out_avg = 'K:/GIS/MODEL/output/temp/' + site_id + '_year_AVG.tif'
            sum_basin_raster.save(out_sum)
            avg_basin_raster.save(out_avg)
            print '12 month SUM saved to ' + out_sum
            print '12 month MEAN saved to ' + out_avg
            break # no need to check the rest of the sites after finding a match
    del cur_basins
    return (dct_site, fields) # dictionary, list





