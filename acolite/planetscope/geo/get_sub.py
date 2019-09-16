## def get_sub
## gets sub in PlanetScope image for given metadata and limit
## written by Quinten Vanhellemont, RBINS
## 2018-03-13
## modifications:  2018-03-15 (QV) fixed bugs

def get_sub(metadata, limit):
    from acolite.planetscope.geo import get_projection

    p, (xscene,yscene) = get_projection(metadata)
    
    dims = metadata["dims"]
    pixelsize = metadata["pixelsize"]

    ## compute x and y limits, round to pixel sizes
    xrange_raw, yrange_raw = p([limit[1],limit[3]],[limit[0],limit[2]])
    xrange = [xrange_raw[0] - (xrange_raw[0] % pixelsize[0]), xrange_raw[1]+pixelsize[0]-(xrange_raw[1] % pixelsize[0])]
    yrange = [yrange_raw[0] - (yrange_raw[0] % pixelsize[1]), yrange_raw[1]+pixelsize[1]-(yrange_raw[1] % pixelsize[1])]
    
    ## crop size
    x_size = int((xrange[1]-xrange[0])/pixelsize[0])
    y_size = int((yrange[0]-yrange[1])/pixelsize[1])

    ## reduce one pixel
    xrange[1] -= pixelsize[0]
    yrange[1] -= pixelsize[1]
        
    grid_region = {'dims':(x_size,y_size), 'xrange':xrange, 'yrange':yrange}

    if (xrange[1] < min(xscene)) or (xrange[0] > max(xscene)):
        print('Limits out of scene longitude')
        return(1)
    elif (yrange[1] < min(yscene)) or (yrange[0] > max(yscene)):
        print('Limits out of scene latitude')
        return(1)
    else:
        xoff = [(i - min(xscene))/pixelsize[0] for i in xrange]
        xoff_region = [(i - min(xscene))/pixelsize[0] for i in xrange]
        
        ## flip the y subset
        yoff = [(i - max(yscene))/pixelsize[1] for i in (yrange[1], yrange[0])]
        yoff_region = [(i - max(yscene))/pixelsize[1] for i in (yrange[1], yrange[0])]

        ## check if the x/y off are within the scene
        if xoff[0] < 0: xoff[0] = 0
        if yoff[0] < 0: yoff[0] = 0
        if xoff[1] >= dims[0]: xoff[1] = dims[0]-1
        if yoff[1] >= dims[1]: yoff[1] = dims[1]-1

        sub = [xoff[0], yoff[0], xoff[1]-xoff[0]-1, yoff[1]-yoff[0]-1]
        sub = [int(s) for s in sub]

        sub_region = [xoff_region[0], yoff_region[0], xoff_region[1]-xoff_region[0]-1, yoff_region[1]-yoff_region[0]-1]
        sub_region = [int(s) for s in sub_region]
        
        off = [sub[0]-sub_region[0], sub[1]-sub_region[1]]
        grid_region['off'] = off
        grid_region['sub'] = sub_region

        if (sub[2] < 0) or (sub[3] < 0):
            print('Problem with computing subset: negative dimensions.')
            return(1)

        return sub, p, (xrange,yrange,grid_region)

