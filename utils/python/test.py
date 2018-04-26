def sortPRISMRasters( RasterList ):
    # Return a 3 dimemnsional list of rasters
    # zero-indexes are left empty for convenience
    # 1st dimension: year (2013 = 1 and so on).
    # 2nd dimension: month (1-12)
    # 3rd dimension: day (1-31)
    # precip for 3/15/2014 = rasters[2][3][15]
    file = 'K:/GIS/MODEL/output/rasterlist.txt'
    try:
        infile = open(file, 'r')
    except:
        pass
    else:
        Rasters = json.load(infile)
        infile.close()
        print 'Read raster list from file: ' + file
        return Rasters

    Rasters = [];
    for i in range(0,4):
        Rasters.append([])
        for j in range (0,13):
            Rasters[i].append([])
            for k in range(0,32):
                Rasters[i][j].append('')

    for raster in RasterList:
        name = raster
        date = name[-16:-8:1]
        year = int( date[:4] )
        year = year - 2012;
        month = int( date[4:6] )
        day = int(date[6:])
        Rasters[year][month][day] = raster

    with open(file, 'w+') as outfile:
        json.dump(Rasters, outfile)
    print 'Dumped raster list to file ' + file
    return Rasters

import os
import arcpy
import csv
import json
## Setup
##

arcpy.env.workspace = 'K:/GIS/MODEL/input/PRISM/precip/daily/'
arcpy.env.scratchWorkspace = 'K:/GIS/MODEL/output'
arcpy.env.overwriteOutput = True
INT_NAN =  -99999

## Copy input sites shapefile to new shapefile which we will use throughout
##
RasterList = arcpy.ListRasters('PRISM*')
RastersSorted = sortPRISMRasters(RasterList)

exampleRas = arcpy.Raster(RastersSorted[1][1][1])
name = exampleRas.name
pixelType = exampleRas.pixelType
sr = exampleRas.spatialReference
fmt = exampleRas.format
print 'Name: ' + name
print 'Pixel type: ' + pixelType
print 'Spatial Reference: ' + sr.name
print 'Raster format: ' + fmt
