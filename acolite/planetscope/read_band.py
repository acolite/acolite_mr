## def read_band
## simple image reading for PlanetScope files
## sub keyword (xoff, yoff, xcount, ycount)
##
## written by Quinten Vanhellemont, RBINS 
## 2018*-03-12
## modifications: 

def read_band(file, band=None, sub=None):
    import os, sys, fnmatch
    if not os.path.isfile(file):
        print('File '+file+' not found.')
        sys.exit()
        
    if fnmatch.fnmatch(file,'*.tif'):
        from osgeo import gdal
        gdal.UseExceptions()
        bandf = gdal.Open(file)
        nrows=bandf.RasterYSize
        ncols=bandf.RasterXSize
        if band is not None:
            bandsel = bandf.GetRasterBand(band)
            if sub is None:
                data = bandsel.ReadAsArray()
            else:
                data = bandsel.ReadAsArray(sub[0],sub[1],sub[2],sub[3])
        else:
            if sub is None:
                data = bandf.ReadAsArray()
            else:
                data = bandf.ReadAsArray(sub[0],sub[1],sub[2],sub[3])

    return(data)
