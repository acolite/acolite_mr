## def fixed_model
## 
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-30
## modifications: QV 2017-01-18 converted to function
##                QV 2017-01-27 added vis bestfit option
##                2017-11-28 (QV) moved PP data directory
##                2018-03-14 (QV) added pressure

def fixed_model(metadata, tau550, rsr_file=None, lutdir=None, lut='PONDER-LUT-201607-MOD2-1013mb', pressure=None):
    import acolite as pp
    import os
    
    ## get scene geometry and default bands from metadata
    try:
        if 'SE_DISTANCE' in metadata.keys():
            se_distance = metadata['SE_DISTANCE']
        else:
            se_distance = pp.distance_se(metadata['DOY'])
            
        ths = metadata['THS']
        thv = metadata['THV']
        azi = metadata['AZI']
        
        if 'SATELLITE_SENSOR' in metadata.keys():
            sensor = metadata['SATELLITE_SENSOR']
        else:
            sensor = metadata['SENSOR']

        #rednir_bands = metadata['BANDS_REDNIR']
        #vis_bands = metadata['BANDS_VIS']
        #nir_bands = metadata['BANDS_NIR']

        #bestfit_bands_defaults = metadata['BANDS_BESTFIT']
        
        bands_sorted = metadata['BANDS_ALL']
    except:
        print('Could not get appropriate metadata for model selection for satellite {}'.format(metadata['SATELLITE']))
        print(metadata.keys())
        return(1)

#    ## get scene geometry
#    if metadata['SATELLITE']=='Pléiades':
#        se_distance = pp.distance_se(metadata['DOY'])
#        ths = 90. - metadata['GEOMETRY'][1]['SUN_ELEVATION']
#        thv = metadata['GEOMETRY'][1]['VIEWING_ANGLE']
#        azi = abs(metadata['GEOMETRY'][1]['SUN_AZIMUTH'] - metadata['GEOMETRY'][1]['AZIMUTH_ANGLE'])
#        if azi > 180: azi -= 180
#        sensor = metadata['INSTRUMENT']+metadata['MISSION_INDEX']
#        rednir_bands = ['Red','NIR']
#        vis_bands = ['Blue','Green','Red']
#    if metadata['SATELLITE']=='WorldView2':
#        se_distance = pp.distance_se(metadata['DOY'])
#        ths = metadata['THS']
#        thv = metadata['THV']
#        azi = metadata['AZI']
#        sensor = 'WorldView2'
#        rednir_bands = ['RED','NIR1','NIR2']
#        vis_bands = ['COASTAL','BLUE','GREEN','YELLOW','RED']

    ## set LUT dir and rsr_file
    #pypath='/storage/Python/' ## to improve!
    #pp_path = os.path.dirname(pp.__file__)
    #if lutdir is None:
    #    lutdir=pp_path+'/data/LUT/'
    #if rsr_file is None:
    #    rsr_file = pp_path+'/data/RSR/'+sensor+'.txt'

    ## set LUT dir and rsr_file
    from acolite import config
    pp_path = config['pp_data_dir']
    if lutdir is None:
        lutdir=pp_path+'/LUT/'
    if rsr_file is None:
        rsr_file = pp_path+'/RSR/'+sensor+'.txt'

    rsr, rsr_bands = pp.rsr_read(file=rsr_file)

    ## get sensor LUT
    #lut_sensor, meta_sensor = pp.aerlut.get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut, override=0)
    #(ratm, rorayl, dtotr, utotr, dtott, utott, astot) = pp.aerlut.lut_get_ac_parameters_fixed_tau_sensor(lut_sensor,meta_sensor,azi,thv,ths,tau550)

    ## get sensor LUT
    lut_sensor, meta_sensor = pp.aerlut.get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut, override=0)

    ## read luts at other pressures if needed
    if pressure is not None:
         lut_data_dict = {}
         lut_data_dict[lut] = {'lut':lut_sensor, 'meta':meta_sensor}
         lut_split = lut.split('-')
         lut0 = '-'.join(lut_split[0:-1]+['0500mb'])
         lut_sensor, meta_sensor = pp.aerlut.get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut0, override=0)
         lut_data_dict[lut0] = {'lut':lut_sensor, 'meta':meta_sensor}

         lut1 = '-'.join(lut_split[0:-1]+['1100mb'])
         lut_sensor, meta_sensor = pp.aerlut.get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut1, override=0)
         lut_data_dict[lut1] = {'lut':lut_sensor, 'meta':meta_sensor}
     
         lut_sensor, meta_sensor = pp.aerlut.aerlut_pressure(lut, lutdir, pressure, sensor, rsr_file, lut_data_dict=lut_data_dict)
    
    (ratm, rorayl, dtotr, utotr, dtott, utott, astot) = pp.aerlut.lut_get_ac_parameters_fixed_tau_sensor(lut_sensor,meta_sensor,azi,thv,ths,tau550)

    return (ratm, rorayl, dtotr, utotr, dtott, utott, astot, tau550), lut_sensor, meta_sensor
