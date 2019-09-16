## def get_rtoa
## reads PlanetScope DN and converts to RTOA
## written by Quinten Vanhellemont, RBINS
## 2018-03-12
## modifications: 2018-03-15 (QV) updated to new metadata
##                2018-03-21 (QV) updated to new metadata for Radiance option

def get_rtoa(file, band, name, metadata, sub=None, radiance=False):
    from acolite.planetscope import read_band
    from numpy import nan, uint16, pi, cos
    
    data = read_band(file, band, sub=sub)
    nodata = data == uint16(0)

    print(metadata['SATELLITE_SENSOR'])
    if (metadata['SATELLITE_SENSOR'] == 'RapidEye') | ('RapidEye' in metadata['SATELLITE_SENSOR']) | (radiance):
        data = data.astype(float) * float(metadata['{}-{}'.format(name,'to_radiance')])
        #print(metadata)
        if not radiance:
            cossza = cos(metadata['THS']*(pi/180.))
            d = metadata['SE_DISTANCE']
            #f0 = metadata['bands'][name]['f0']*10.
            f0 = metadata['{}-f0'.format(name)]*10.
            data *= (pi * d * d) / (f0 * cossza)
    else:
        data = data.astype(float) * float(metadata['{}-{}'.format(name,'to_reflectance')])
    
    data[nodata] = nan
    return(data)
