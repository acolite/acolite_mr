## def get_rtoa
## reads Pléiades DN and converts to RTOA - needs metadata and sun zenith
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-06-29
## modifications: 2017-03-16 (QV) added nc option (Pléiades only at the moment)
##                2017-05-03 (QV) added new reflectance option for Pléiades
##                2017-11-21 (QV) added ths and se_distance from metadata
##                2017-12-06 (QV) added gains option

def get_rtoa(file, idx, name, metadata, sun_zenith=0.0, se_distance=1.0, sub=None, gains=None):
    from acolite.pleiades import read_band
    from acolite import dn_to_rtoa, distance_se
    from numpy import nan, uint16, pi,cos,nanmedian, asarray

    if '.nc' in file:
        from ponder_processor.shared.nc_read import nc_data
        ## add check here for different sensors
        dsets = {'Blue':'rhot_495', 'Green':'rhot_558','Red':'rhot_655','NIR':'rhot_842',
                 'B0':'rhot_495', 'B1':'rhot_558','B2':'rhot_655','B3':'rhot_842'}
        data = nc_data(file, dsets[name], sub=sub)
    else:
        if 'THS' in metadata.keys():
            sun_zenith = metadata['THS']
        if 'SE_DISTANCE' in metadata.keys():
            se_distance = metadata['SE_DISTANCE']

        #sun_zenith = 90. - metadata['GEOMETRY'][1]['SUN_ELEVATION']
        #se_distance = distance_se(metadata['DOY'])
        #se_distance=1.0
        data = read_band(file, idx=idx, sub=sub)
        #saturated = data == uint16(metadata['SATURATED'])
        nodata = data == uint16(metadata['NODATA'])

        if (metadata['RADIOMETRIC_PROCESSING'] == 'RADIANCE') | (metadata['RADIOMETRIC_PROCESSING'] == 'BASIC'):
            data = dn_to_rtoa(data, metadata['BAND_INFO'][name]['F0'], sun_zenith, 
                        slope=1./metadata['BAND_INFO'][name]['radiance_gain'], offset=metadata['BAND_INFO'][name]['radiance_bias'], d=se_distance)

            #test=dn_to_rtoa(asarray((0,1,2,3,4,5,6,7,8,9,10)), metadata['BAND_INFO'][name]['F0'], sun_zenith, 
            #            slope=1./metadata['BAND_INFO'][name]['radiance_gain'], offset=metadata['BAND_INFO'][name]['radiance_bias'], d=se_distance)

        elif (metadata['RADIOMETRIC_PROCESSING'] == 'LINEAR_STRETCH'):
            data = data.astype('float')
            data = dn_to_rtoa(data, metadata['BAND_INFO'][name]['F0'], sun_zenith, 
                        slope=1./metadata['BAND_INFO'][name]['radiance_gain'], offset=metadata['BAND_INFO'][name]['radiance_bias'], d=se_distance)
            print('Warning linear stretch data')
        else:
            print("{} RADIOMETRIC_PROCESSING not recognised".format(metadata['RADIOMETRIC_PROCESSING']))

        if False:
            print(name)
            print(idx)
            print(metadata['BAND_INFO'][name])
            print('F0',metadata['BAND_INFO'][name]['F0'])
            print('sun_zenith',sun_zenith)
            print('se_distance',se_distance)
            print('slope',1./metadata['BAND_INFO'][name]['radiance_gain'])
            print('offset',metadata['BAND_INFO'][name]['radiance_bias'])
            print('digitization', test[1]-test[0])
            print(test)
            print(test[1:9]-test[0:8])

        if metadata['RADIOMETRIC_PROCESSING'] == 'REFLECTANCE':
            data=data.astype('float')
            #RHO=RHO'/GAIN+BIAS
            data /= metadata['BAND_INFO'][name]['reflectance_gain']
            data += metadata['BAND_INFO'][name]['reflectance_bias']
            dtor = pi/180.
            data /= cos(sun_zenith*dtor)
            #L=RHO/GAIN+BIAS
            #data /= metadata['BAND_INFO'][name]['radiance_gain']
            #data += metadata['BAND_INFO'][name]['radiance_bias']
            #dtor = pi/180.
            #print(nanmedian(data))
            #print(sun_zenith)
            #print(metadata['BAND_INFO'][name]['F0'])
            #data *= (pi * se_distance * se_distance) / (metadata['BAND_INFO'][name]['F0'] * cos(sun_zenith*dtor))

        data[nodata] = nan

    ## can be cleaner
    if gains is not None:
        bidx= [i for i, j in enumerate(metadata["BANDS"]) if j == name]
        if len(bidx) == 0:
            bidx = [i for i, j in enumerate(metadata["BAND_NAMES"]) if j == name]
        if len(bidx) == 1: bname = metadata["BAND_NAMES"][bidx[0]]
        if bname in gains:
            data*=gains[bname]
            print('Gain {} applied for band "{}" at TOA.'.format(gains[bname], bname)) 
    return data
