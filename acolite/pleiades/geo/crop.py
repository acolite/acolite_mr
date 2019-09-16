## def crop
## finds crop position (x0, y0, ns, nl) for Pleiades image
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-17
## modifications:
##

def crop(metadata, limit):
    from math import ceil, floor
    from acolite.pleiades import geo
    
    ## get interpolator for row / col
    _, (zcol, zrow) = geo.init(metadata)
    
    south = limit[0]
    east = limit[1]
    north = limit[2]
    west = limit[3]

    se_x = floor(zcol(east,south)[0])
    se_y = ceil(zrow(east,south)[0])

    nw_x = ceil(zcol(west,north)[0])
    nw_y = floor(zrow(west,north)[0])

    ns = nw_x - se_x
    nl = se_y - nw_y

    sub = [se_x, nw_y, ns, nl]
    return sub
