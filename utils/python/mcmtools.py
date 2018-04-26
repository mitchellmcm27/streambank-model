import csv
import arcpy
import json
import os
import datetime
from dateutil.relativedelta import relativedelta

def make_sure_path_exists(path):
    """Create a directory only if needed."""
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

def shp2csv(shp, fields_shp, csv_path):
    """
    Write the given fields of a shapefile to a .csv file

    Args:
        shp: string shp = 'path\filename of input.shp'
        fields_shp: list fields_shp = [list of fields to write to csv]
        csv_path: string csv_path = 'path\filename of output.csv'
    """
    with open(csv_path,'w+b') as f:
        writer = csv.writer(f)
        writer.writerow(fields_shp) # write headers
        with arcpy.da.SearchCursor(shp, fields_shp) as cur:
            for row in cur:
                writer.writerow(row[:])

def append_dict2csv(dic, csv_path):
    """
    Append the given keys of a dict to a .csv file.
    """
    with open(csv_path,'a') as f:
        writer = csv.DictWriter(f, sorted(dic.keys()), restval='NaN')
        writer.writerow(dic)

def write_dict2csv(dic, csv_path):
    """
    Write the given keys of a dict to a .csv file, creating it if needed
    """
    with open(csv_path,'w+b') as f:
        writer = csv.DictWriter(f, sorted(dic.keys()), restval='NaN')
        writer.writerow(dic)

def dates_from_site(date1):
    """
    Convert 2 formatted datestrings into lists of ints
    -input format must match 'dd MMM yyyy', e.g. '01 Jan 2014'
    -output format list = [yy][mm][dd]

    Spans one month before to 12 months after given date (1 year with the month preceeding)
    """
    if ( date1 == 'NaT' ):
        return
    starting_day = int(date1[:2])
    starting_month = convert_month(date1[3:6])
    starting_year = int(date1[7:])

    # Force ending date to be 1 year after start
    if starting_day > 15: # if more than halfway through the month,
        starting_month += 1 # start at the first of the next month
    if starting_month > 12:
        starting_month = 1 # go to January of next year if needed
        starting_year += 1
    starting_day = 1

    starting_date = datetime.date(starting_year, starting_month, starting_day)
#    starting_date = starting_date + relativedelta(months=-1)
    ending_date = starting_date + relativedelta(months=+12)
    return (starting_date, ending_date)

def translate_date(date_str):
    """
    Convert 2 formatted datestrings into lists of ints
    -input format must match 'dd MMM yyyy', e.g. '01 Jan 2014'
    -output format list = [yy][mm][dd]

    Spans one month before to 12 months after given date (1 year with the month preceeding)
    """
    if ( date_str == 'NaT' ):
        return
    dy = int(date_str[:2])
    mn = convert_month(date_str[3:6])
    yr = int(date_str[7:])

    dt_date = datetime.date(yr, mn, dy)
    return dt_date

def convert_month(inpt): # switch between 2-int or 3-char month notation
    """
    Convert between int and string representations of months ('Jan' <-> 01).
    """
    month_str = {'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4,
                'May':5, 'Jun':6, 'Jul':7, 'Aug':8,
                'Sep':9, 'Oct':10, 'Nov':11, 'Dec':12}
    month_num = {1:'Jan', 2:'Feb', 3:'Mar', 4:'Apr',
                5:'May', 6:'Jun', 7:'Jul', 8:'Aug',
                9:'Sep', 10:'Oct', 11:'Nov', 12:'Dec'}
    try:
        ans = month_str[inpt]
    except:
        ans = month_num[int(inpt)]
    return ans

def get_monthly_prism_raster(year,month):
    """
    Find monthly PRISM precip raster corresponding to given year and month.

    Arguments:
        year: 4 digit int
        month: 2 digit int

    Returns string of matching raster file.
    """
    if year==0:
        wildcard = '*%02d*' % (month)
    else:
        wildcard = '*%04d%02d*' % (year, month) # matches default PRISM naming
    raster_list = arcpy.ListRasters(wildcard)
    raster = raster_list[0] # but there should only be one raster in the list
    return raster

def sort_prism_rasters( raster_list_in ):
    """
   Arrange a bunch of unsorted PRISM rasters more logically.

    Arguments:
    -raster_list_in: a list of strings (raster filenames)

    Returns a 3-dimensional list of strings, where
        zero-indexes are left empty for convenience.
        1st dimension: year (2013 = 1 and so on).
        2nd dimension: month (1-12)
        3rd dimension: day (1-31)
    Example: precip for 3/15/2014 = rasters[2][3][15]

    """

    print 'Creating raster list...'
    ## Create 3-dimensional python list for storing rasters
    ##
    rasters = [];
    for i in range(0,4):
        rasters.append([])
        for j in range (0,13):
            rasters[i].append([])
            for k in range(0,32):
                rasters[i][j].append('')
    ## Loop through unsorted rasters
    ##
    for raster in raster_list_in:
        # Get dates from filename and place raster into the correct slot
        name = raster
        date = name[-16:-8:1]
        year = int( date[:4] )
        year = year - 2012;
        month = int( date[4:6] )
        day = int(date[6:])
        rasters[year][month][day] = raster
    return rasters