## def read_band
## simple image reading
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-06-29
## modifications: 2016-06-30 QV added JP2 support via the glymur package
##                2016-07-06 QV changed indexing for glymur
##                2016-10-21 QV glymur is broken, switched to gdal for JP2 reading
##                2016-11-23 QV added sub keyword (xoff, yoff, xcount, ycount)

def read_band(file, idx=None, sub=None):
    import os, sys, fnmatch
    if not os.path.isfile(file):
        print('File '+file+' not found.')
        sys.exit()
        
    if fnmatch.fnmatch(file,'*.TIF'):
        from osgeo import gdal
        gdal.UseExceptions()
        ds = gdal.Open(file)
        if idx != None: band = ds.GetRasterBand(idx)
        else: band = ds.GetRasterBand()
        nrows=band.YSize
        ncols=band.XSize
        if sub is None:
             data = band.ReadAsArray()
        else:
             data = band.ReadAsArray(sub[0],sub[1],sub[2],sub[3])
        ds = None
        
    #if fnmatch.fnmatch(file,'*.JP2'):
    #    import glymur
    #    ds = glymur.Jp2k(file)
    #    if idx != None: data =  ds[:,:,idx-1]
    #    else: data = ds[:]
    #    ds = None
    
    if fnmatch.fnmatch(file,'*.JP2'):
        from osgeo import gdal
        gdal.UseExceptions()
        ds = gdal.Open(file)
        if idx != None: band = ds.GetRasterBand(idx)
        else: band = ds.GetRasterBand()
        nrows=band.YSize
        ncols=band.XSize
        if sub is None:
             data = band.ReadAsArray()
        else:
             data = band.ReadAsArray(sub[0],sub[1],sub[2],sub[3])
        ds = None
        
    return data
