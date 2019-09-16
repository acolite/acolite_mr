## def read_band
## WV2 image reading
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-27
## modifications:
##

def read_band(file, idx, sub=None):
    import os, sys, fnmatch
    if not os.path.isfile(file):
        print('File '+file+' not found.')
        sys.exit()
        
    if fnmatch.fnmatch(file,'*.TIF'):
        from osgeo import gdal
        gdal.UseExceptions()
        ds = gdal.Open(file)
        band = ds.GetRasterBand(idx)
        nrows=band.YSize
        ncols=band.XSize
        if sub is None:
             data = band.ReadAsArray()
        else:
             data = band.ReadAsArray(sub[0],sub[1],sub[2],sub[3])
        ds = None
                
    return data
