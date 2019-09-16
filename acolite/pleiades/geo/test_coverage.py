## def test_coverage
## tests if ROI given by limit [S,W,N,E] is covered by given image
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-17
## modifications: QV 2017-03-23 added quick and dirty subset to returned file lists
##

def test_coverage(file, limit, verbose=False):
    import acolite as pp
    
    #scname = os.path.basename(file)

    ## check if we are dealing with Pléiades image bundle
    ifile,mfile = pp.pleiades.image_test(file)

    ## 2017-03-23 added subset to the list
    ifile = ifile[0]
    mfile = mfile[0]
    print('test coverage subset the pleiades image test output lists - this is a very clear error message')

    ## parse metadata
    metadata = pp.pleiades.parse_metadata(mfile)

    ## check for coverage
    sub = pp.pleiades.geo.crop(metadata, limit)
    
    ncols = float(metadata['NCOLS'])
    nrows = float(metadata['NROWS'])

    yu=0
    yl=0
    xl=0
    xr=0
    
    ### check longitudes
    if sub[0] >= ncols:
        if verbose: print('Western longitude out of eastern bound of scene.')
        xl=1
    
    if (sub[0]+sub[2]) >= ncols:
        if verbose: print('Eastern longitude out of eastern bound of scene.')
        xr=1

    if sub[0] <= 1:
        if verbose: print('Eastern longitude out of western bound of scene.')
        xl=-1
    
    if (sub[0]+sub[2]) <= 1:
        if verbose: print('Western longitude out of western bound of scene.')
        xr=-1

    ## check latitudes
    if sub[1] <= 1:
        if verbose: print('Northern longitude out of northern bound of scene.')
        yu=1
        
    if (sub[1]+sub[3]) <= 1:
        if verbose: print('Southern longitude out of northern bound of scene.')
        yl=1
        
    if sub[1] >= nrows:
        if verbose: print('Northern longitude out of southern bound of scene.')
        yu=-1
        
    if (sub[1]+sub[3]) >= nrows:
        if verbose: print('Southern longitude out of southern bound of scene.')
        yl=-1

    ## check dimensions
    if (sub[2] == 1) or (sub[3] == 1):
        if verbose: print('Crop dimensions of 1: region probably out of scene.')
        yl=-1
        yu=-1
        xl=-1
        xr=-1    

    if (yl == 0) & (yu == 0) & (xl == 0) & (xr == 0):
        if verbose: print('Region fully in scene.')
        return True

    if (xl == 1):
        if verbose: print('Region fully west of scene')
        return False

    if (xr == -1):
        if verbose: print('Region fully east of scene')
        return False

    if (yl == 1):
        if verbose: print('Region fully north of scene')
        return False

    if (yu == -1):
        if verbose: print('Region fully south of scene')
        return False

    if verbose: print('Region partially within bounds of scene')
    #print(sub)

    return True
