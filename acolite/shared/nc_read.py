# read dataset and global attributes from netcdf
def nc_read(file, dataset):
    from netCDF4 import Dataset
    nc = Dataset(file)
    gatts = {attr : getattr(nc,attr) for attr in nc.ncattrs()}
    out_array = nc.variables[dataset][:]
    nc.close()
    return (out_array, gatts)

# read dataset from netcdf
# Last updates: 2016-12-19 (QV) added crop (x0,x1,y0,y1)
##              2017-03-16 (QV) added sub keyword (xoff, yoff, xcount, ycount)
def nc_data(file, dataset, crop=False, sub=None):
    from netCDF4 import Dataset
    nc = Dataset(file)
    if sub is None:
        if crop is False:
            data = nc.variables[dataset][:]
        else:
            if len(crop) is 4: data = nc.variables[dataset][crop[2]:crop[3]:1,crop[0]:crop[1]:1]
            else: data = nc.variables[dataset][:]
    else:
        if len(sub) is 4: data = nc.variables[dataset][sub[1]:sub[1]+sub[3]:1,sub[0]:sub[0]+sub[2]:1]
        else: data = nc.variables[dataset][:]
    nc.close()
    return data

# read dataset and global attributes from netcdf
def nc_gatts(file):
    from netCDF4 import Dataset
    nc = Dataset(file)
    gatts = {attr : getattr(nc,attr) for attr in nc.ncattrs()}
    nc.close()
    return gatts

# read datasets in netcdf
def nc_datasets(file):
    from netCDF4 import Dataset
    nc = Dataset(file)
    ds = list(nc.variables.keys())
    nc.close()
    return ds
