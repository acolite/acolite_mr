## def acolite_mr_ac
## 
## written by Quinten Vanhellemont, RBINS for the PONDER project
## Sep 2019
## modifications: QV 2019-09-16 initial version

def acolite_mr_ac(bundle, 
                  output=None, limit=None,
                  ancillary_data=True,
                  luts=['PONDER-LUT-201704-MOD1-1013mb', 'PONDER-LUT-201704-MOD2-1013mb'], 
                  uoz=0.3,
                  uwv=1.5,
                  map_rgb=True, map_rgb_rhos=True,
                  pan_sharpen_rgb=False,
                  sky_correction=True,
                  dem_pressure=False, 
                  force_band=None,
                  fixed_aot550=None, fixed_lut=None,
                  dark_spectrum_full_scene=True):

    import acolite as ac

    ## can we identify the bundle
    ## does not work with zipped bundles
    data_type = None
    for dt in ['pleiades', 'planet']:
        if data_type == None:
            try: 
                if dt == 'pleiades':
                    tmp = ac.pleiades.image_test(bundle, listpan=True)
                    data_type = str(dt)
                if dt == 'planet':
                    tmp = ac.planetscope.bundle_test(bundle)
                    data_type = str(dt)
            except:
                pass
        else: continue

    ## exit if bundle not recognised
    if data_type is None:
        print('Input bundle {} not recognised'.format(bundle))
        return()
    else:
        print('Processing {} file'.format(data_type))

    
    ## Pléiades/SPOT
    if data_type == 'pleiades':
        ac.pleiades.pleiades_ac(bundle, output=output, limit=limit, 
                            ancillary_data=ancillary_data,
                            luts=luts, write_netcdf_geo=False,
                            uoz=uoz,
                            uwv=uwv,
                            map_rgb=map_rgb, map_rgb_rhos=map_rgb_rhos,
                            pan_sharpen_rgb=pan_sharpen_rgb,
                            sky_correction=sky_correction,
                            dem_pressure=dem_pressure, 
                            force_band=force_band,
                            fixed_aot550=fixed_aot550, fixed_lut=fixed_lut,
                            dark_spectrum_full_scene=dark_spectrum_full_scene)
        return()

    ## PlanetScope/RapidEye
    if data_type == 'planet':
        ac.planetscope.planetscope_ac(bundle, output=output, limit=limit,
                            ancillary_data=ancillary_data,
                            luts=luts, 
                            uoz_default=uoz,
                            uwv_default=uwv,
                            map_rgb=map_rgb, map_rgb_rhos=map_rgb,
                            sky_correction=sky_correction)
        return()



    
