## def cover_test
## checks if ROI is within Pleiades image
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-18
## modifications: QV 2017-03-23 added quick and dirty subset to returned file lists
##

def cover_test(file, limit, region='', verbose=False):
    import sys,os
    import acolite as pp

    ret=pp.pleiades.geo.test_coverage(file, limit, verbose=verbose)
    scname = os.path.basename(file)
    if ret is True:
        ## check if we are dealing with Pléiades image bundle
        ifile,mfile = pp.pleiades.image_test(file)
      
        ## 2017-03-23 added subset to the list
        ifile = ifile[0]
        mfile = mfile[0]
        print('test coverage subset the pleiades image test output lists - this is a very clear error message')

        ## parse metadata
        metadata = pp.pleiades.parse_metadata(mfile)

        print('Region {} in scene {}:'.format(region, scname))
        print('Date/time: {} {}'.format(metadata['IMAGING_DATE'], metadata['IMAGING_TIME']))

        print('\tlon: {}; {}'.format(limit[1],limit[3]))
        print('\tlat: {}; {}'.format(limit[0],limit[2]))



        print('\tImage dimensions: {}x{}'.format(metadata['NCOLS'], metadata['NROWS']))

        sub = pp.pleiades.geo.crop(metadata, limit)
        print('\tCrop position: {}; {}'.format(sub[0],sub[1]))
        print('\tCrop dimensions: {}x{}'.format(sub[2],sub[3]))
        return True
    return False
