# Import system modules
import arcpy
from arcpy import env
from arcpy.sa import *

# Set environment settings
env.workspace = "K:/GIS/r/elevation/"
arcpy.CheckOutExtension("Spatial")
# Set local variables
for raster in arcpy.ListRasters():
    outSlope = Slope(raster, "PERCENT_RISE", 1)
    outpath = env.workspace + "/slope/" + arcpy.Describe(raster).baseName + "_slope.tif"
    print(outpath)
    outSlope.save(outpath)
print("finished")
