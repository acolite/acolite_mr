## def get_ll
## gets lon and lat for given PlanetScope metadata with optional limit
## written by Quinten Vanhellemont, RBINS
## 2018-03-13
## modifications: 2020-09-15 (QV) full tile y coordinates do not need to be flipped
##

def get_ll(metadata, limit=None, xy=False, extend_limit=False):
    from numpy import linspace, tile, flipud
    from acolite.planetscope.geo import get_sub
    from acolite.planetscope.geo import get_projection

    pixelsize = metadata["resolution"]

    if limit is not None:
        sub, p, (xrange, yrange, grid_region) = get_sub(metadata, limit)
        if extend_limit is False:
            dims = [sub[2],sub[3]]
        else:
            dims = grid_region['dims']
            xrange, yrange = grid_region['xrange'], grid_region['yrange']
    else:
        p, (xrange,yrange) = get_projection(metadata)
        dims = metadata["dims"]

    xdim =  linspace(xrange[0],xrange[1],dims[0]).reshape(1,dims[0])
    ydim =  linspace(yrange[0],yrange[1],dims[1]).reshape(dims[1],1)

    xdim = tile(xdim, (dims[1],1))
    if limit is not None:
        ydim = flipud(tile(ydim, (1,dims[0])))
    else:
        ydim = tile(ydim, (1,dims[0]))

    lon,lat = p(xdim,ydim,inverse=True)

    if xy:
        return(lon,lat,xdim,ydim)
    else:
        xdim = None
        ydim = None
        return(lon,lat)
