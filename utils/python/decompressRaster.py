def main():
    import arcpy
    from arcpy import env
    from arcpy.sa import *
    # Set environment settings
    env.workspace = "K:/GIS/r/elevation/"
    env.compression_type = "NONE"
if __name__ == '__main__':
    main()

for raster in arcpy.ListRasters():
    try:
        out_raster = env.workspace + "/exported/" + arcpy.Describe(raster).baseName + ".tif"
        arcpy.management.CopyRaster(raster, out_raster)
    except:
        print "Failed to copy " + arcpy.Describe(raster).baseName
        print arcpy.GetMessages()
print("finished")
