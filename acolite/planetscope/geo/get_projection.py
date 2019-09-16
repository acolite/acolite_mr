## def get_projection
## get projection from PlanetScope metadata
##
## written by Quinten Vanhellemont, RBINS 
## 2018-03-13
## modifications:  2018-03-15 (QV) added retrieval from GeoTiff image (bugs in other method)

def get_projection(metadata):
    from pyproj import Proj

    if 'image_file' in metadata:
        from osgeo import gdal,osr
        ds   = gdal.Open(metadata['image_file'])
        transform = ds.GetGeoTransform()
        projection = ds.GetProjection()
        src = osr.SpatialReference()
        src.ImportFromWkt(projection)
        proj4_string = src.ExportToProj4()
        p = Proj(proj4_string)
        
        x0 = transform[0]
        dx = transform[1]
        y0 = transform[3]
        dy = transform[5]
        pixelsize = (dx, dy)
        dimx, dimy = metadata['dims']
        metadata['pixelsize'] = pixelsize
        xrange = (x0,x0+dimx*dx)
        yrange = (y0,y0+dimy*dy)
    else:
        print('Warning: projection not derived from GeoTiff - has bugs')
        pixelsize = (metadata["resolution"][0],metadata["resolution"][1]*-1)
        proj = metadata['projection'].split()[3]
        datum = metadata['datum'].replace('_19','')

        if (proj == 'UTM') & (datum == 'WGS84'): is_utm = True
        else: is_utm = False

        if (not is_utm): print('Projection not implemented')

        if is_utm:
            zone = int(metadata['zone'])
            zonestr = metadata['projection'].split()[5]
            proj4_list = ['+proj=utm',
                          '+zone={}'.format(zone),
                          '+datum={}'.format(datum),
                          '+units=m',
                          '+no_defs ']
            if 'S' in zonestr: proj4_list+=['+south']

        proj4_string = ' '.join(proj4_list)
    
        p = Proj(proj4_string)

        ## check corners of image
        x,y = [],[]
        lons,lats = [],[]
        for corner in ['LL','UL','UR','LR']:
                lon = metadata['{}_LONGITUDE'.format(corner)]
                lat = metadata['{}_LATITUDE'.format(corner)]
                lons.append(lon)
                lats.append(lat)
                xc, yc = p(lon,lat)
                x.append(xc)
                y.append(yc)
        xrange = [min(x),max(x)]
        yrange = [min(y),max(y)]
        
    return(p, (xrange,yrange))

