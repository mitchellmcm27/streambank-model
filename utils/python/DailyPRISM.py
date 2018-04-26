def make_sure_path_exists( path ):
    """
    Function to create a directory and file if needed.
    Does nothing if it already exists.
    """
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

def sortPRISMRasters( RasterList ):
    """
    Function to arrange a  bunch of unsorted PRISM rasters, while
        first checking to see if a txt file of arranged rasters exists.
    -input: RasterList, a list of strings (filenames)
    -zero-indexes are left empty for convenience.
    -1st dimension: year (2013 = 1 and so on).
    -2nd dimension: month (1-12)
    -3rd dimension: day (1-31)
    -precip for 3/15/2014 = rasters[2][3][15]

    """
    file = 'K:/GIS/MODEL/output/temp/rasterlist.txt'
    try:
        # I'll try to get the arranged rasters from the text file
        infile = open(file, 'r')
    except:
        # Boo, we have to do the work ourselves
        pass
    else:
        # Yay, I'll just load those rasters for you...
        Rasters = json.load(infile)
        infile.close()
        print 'Read raster list from file: ' + file
        return Rasters

    ## Create 3-dimensional python list for storing rasters
    ##
    Rasters = [];
    for i in range(0,4):
        Rasters.append([])
        for j in range (0,13):
            Rasters[i].append([])
            for k in range(0,32):
                Rasters[i][j].append('')
    ## Loop through unsorted rasters
    ##
    for raster in RasterList:
        # Get dates from filename and place raster into the correct slot
        name = raster
        date = name[-16:-8:1]
        year = int( date[:4] )
        year = year - 2012;
        month = int( date[4:6] )
        day = int(date[6:])
        Rasters[year][month][day] = raster

    ## Dump the arranged list to a json for future loading (faster!)
    ##
    with open(file, 'w+') as outfile:
        json.dump(Rasters, outfile)
    print 'Dumped raster list to file ' + file
    return Rasters

def getMonthlyRaster(year,month):
    """
    Function that, given a year and a month, returns the corresponding
        monthly PRISM raster (filename as a string)
    """
    wildcard = '*%04d%02d*' % (year, month) # matches default PRISM naming
    RasterList = arcpy.ListRasters(wildcard)
    raster = RasterList[0] # but there should only be one raster in the list
    return raster

def datesFromSite(date1, date2):
    """
    Function to convert 2 formatted datestrings into lists of ints
    -input format must match 'dd MMM yyyy', e.g. '01 Jan 2014'
    -output format list = [yy][mm][dd]
    """
    if ( date1 == 'NaT' ):
        return
    starting_day = int(date1[:2])
    starting_month = convertMonth(date1[3:6])
    starting_year = int(date1[7:])
    ending_day = int(date2[:2])
    ending_month = convertMonth(date2[3:6])
    ending_year = int(date2[7:])

    StartingDate = [starting_year, starting_month, starting_day]
    EndingDate = [ending_year, ending_month, ending_day]
    return (StartingDate, EndingDate)
def calcPrecip( RastersSorted, StartingDate, EndingDate, basin_id_site, site_id): # called for each site
    """
    Function to calculate precipitation statistics over a given area throughout
        given dates.
    -input: RastersSorted - a list created by the above function (all rasters)
    -input: StartingDate, EndingDate - dates created by the above function
    -input: basin_id_site - BasinID from study site shapefile, will be matched
    -input: site_id - unique identifier of study site, e.g. PR01Y1

    -output: dct_site - a python dictionary with the calculations for input site id
    -output: Fields - python list of fields created to store calculations
    """
    print 'Starting date: '
    print StartingDate
    print 'Ending date:'
    print EndingDate
    print 'BasinID from site:'
    print basin_id_site

    starting_year = StartingDate[0] - 2012
    starting_month = StartingDate[1]
    starting_day = StartingDate[2]
    ending_year = EndingDate[0] - 2012
    ending_month = EndingDate[1]
    ending_day = EndingDate[2]
    RastersAll = []
    dct_site = {}
    # Select the drainage basin matching BasinID from the site
    cBasins = arcpy.da.SearchCursor(basins_shp, 'BasinID')
    Basins = [row[0] for row in cBasins]
    day = 0
    for basin in Basins:

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
                    # Use the daily rasters. They are found in...
                    arcpy.env.workspace = \
                        'K:/GIS/MODEL/input/PRISM/precip/daily'
                    for d in range(d1, d2 + 1):
                        if RastersSorted[y][m][d] != '':
                            day += 1
                            raster = RastersSorted[y][m][d]

                        ## Do the calculations
                        ##
                        # Output file setup
                        date_code = site_id + str(day)
                        clipped = 'K:/GIS/MODEL/output/temp/' + date_code + '_clip.tif'
                        resampled = 'K:/GIS/MODEL/output/temp/' + date_code + '_res.tif'
                        table = 'K:/GIS/MODEL/output/tables/' + site_id
                        try:
                            # if I already have the needed data, save a little time
                            inraster = arcpy.Raster(resampled)
                        except:
                            # this is the first time I'm doing these calculations
                            print 'Clipping...'
                            arcpy.Clip_management(raster,clip_bbox,clipped, '#', '#', 'NONE')
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
                                    if field not in Fields:
                                        Fields.append(field)
                                    key = (site_id, field) # e.g. ('PR01Y1', 'SUM_01')
                                    dct_site[key] = row[k]
                                    k += 1
                        del RastersMonth[:]

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
    del cBasins
    return (dct_site, Fields) # dictionary, list

def convertMonth( inpt ): # switch between 2-int or 3-char month notation
    """
    Converts between int and string representations of months
    """
    MonthStr = {'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4,
                'May':5, 'Jun':6, 'Jul':7, 'Aug':8,
                'Sep':9, 'Oct':10, 'Nov':11, 'Dec':12}
    MonthNum = {1:'Jan', 2:'Feb', 3:'Mar', 4:'Apr',
                5:'May', 6:'Jun', 7:'Jul', 8:'Aug',
                9:'Sep', 10:'Oct', 11:'Nov', 12:'Dec'}
    try:
        ans = MonthStr[inpt]
    except:
        ans = MonthNum[int(inpt)]
    return ans

def shp2CSV( shp, Fields_shp, csv_path ):
    """
    Function to write the given fields of a shapefile to a .csv file
    -input: string shp = 'path\filename of input.shp'
    -input: list Fields_shp = [list of fields to write to csv]
    -input: string csv_path = 'path\filename of output.csv'
    """
    with open(csv_path,'w+b') as f:
        writer = csv.writer(f)
        writer.writerow(Fields_shp) # write headers
        with arcpy.da.SearchCursor(shp, Fields_shp) as cur:
            for row in cur:
                writer.writerow(row[:])

## Main Program
"""
Module to match up study sites, drainage basins, and PRISM precipitation data,
calculating statistics on the precip data, and outputting the results to
a shapefile, rasters, and a csv.

As an example, the total and mean precipitation during a given 12-month period
in a given area (drainage basin) can be determined.
"""
## Modules
##
import os
import arcpy
import csv
import json

## Setup
##
# Required folders
RequiredPaths  = ['K:/GIS/MODEL/output',
                  'K:/GIS/MODEL/output/temp',
                  'K:/GIS/MODEL/output/rasters',
                  'K:/GIS/MODEL/output/tables']
for path in RequiredPaths:
    make_sure_path_exists(path)

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
Fields = [] # will contain a list of fields eventually created by Zonal Statistics

# Copy input sites shapefile to new shapefile, which we will use throughout
arcpy.CopyFeatures_management(in_shp, out_shp)

# Let's get those rasters...
RasterList = arcpy.ListRasters('PRISM*')
# ... and output some useful properties
exampleRas = arcpy.Raster(RasterList[0])
pixelType = exampleRas.pixelType
sr = exampleRas.spatialReference
fmt = exampleRas.format
print 'Pixel type: ' + pixelType
print 'Spatial Reference: ' + sr.name
print 'Raster format: ' + fmt

# Sort the rasters so we can use 'em
RastersSorted = sortPRISMRasters(RasterList)

## Loop through points in 'sites' shapefile
##
with arcpy.da.SearchCursor(out_shp, ['BasinID', 'date_y1','date_y2', 'id']) as cur:
    i = 0
    for row in cur:
        try:
            # See if the site has 2 valid dates
            StartingDate, EndingDate = datesFromSite(row[1], row[2])
            basin_id = int(row[0])
            site_id = row[3]
        except:
            # No dates? No bueno. Go to next site.
            continue
        print 'Processing site ' + site_id

        # Let's crunch the numbers
        # Return little dict of site's results
        dct_site, Fields = calcPrecip(RastersSorted, StartingDate,
                                       EndingDate, basin_id, site_id)
        dct_results.update(dct_site) # Append calculation results to big dict
        i += 1
        if i == 3:
            pass # in case you wanna check on some stuff
## output shapefile data to a CSV
##
shp2CSV(out_shp, Fields, 'K:/GIS/MODEL/output/results.csv')

arcpy.CheckInExtension('Spatial')


