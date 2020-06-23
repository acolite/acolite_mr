## def acolite_mr_ac
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## Sep 2019
## modifications: QV 2019-09-16 initial version
##                2020-02-25 (QV) added ignore_sr_image keyword
##                2020-03-16 (QV) added worldview
##                2020-05-19 (QV) added elevation/pressure option

def acolite_mr_ac(bundle,
                  output=None, limit=None,
                  ancillary_data=True,
                  luts=['PONDER-LUT-201704-MOD1-1013mb', 'PONDER-LUT-201704-MOD2-1013mb'],
                  uoz=0.3,
                  uwv=1.5,
                  ignore_sr_image = True,
                  map_rgb=True, map_rgb_rhos=True,
                  pan_sharpen_rgb=False,
                  sky_correction=True,
                  dem_pressure=False, elevation=None, pressure=None,
                  force_band=None,
                  fixed_aot550=None, fixed_lut=None,
                  dark_spectrum_full_scene=True):

    import acolite as ac
    import glob

    ## can we identify the bundle
    ## does not work with zipped bundles
    data_type = None
    for dt in ['pleiades', 'planet', 'worldview']:
        if data_type == None:
            try:
                if dt == 'pleiades':
                    tmp = ac.pleiades.image_test(bundle, listpan=True)
                    data_type = str(dt)
                if dt == 'planet':
                    tmp = ac.planetscope.bundle_test(bundle)
                    if len(tmp) > 0:
                        data_type = str(dt)
                if dt == 'worldview':
                    metafile = glob.glob('{}/{}'.format(bundle,'*.XML'))[0]
                    tmp = ac.worldview.parse_metadata(metafile)
                    if len(tmp) > 0:
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

    ## compute pressure from elevation if given
    if (pressure is None) & (elevation is not None):
        elevation = float(elevation)
        pressure = ac.ac.pressure_elevation(elevation)
        print('Computed {:.2f} hPa pressure for {:.2f} m elevation'.format(pressure, elevation))

    ## Pléiades/SPOT
    if data_type == 'pleiades':
        ac.pleiades.pleiades_ac(bundle, output=output, limit=limit,
                            ancillary_data=ancillary_data,
                            luts=luts,
                            uoz=uoz,
                            uwv=uwv,
                            map_rgb=map_rgb, map_rgb_rhos=map_rgb_rhos,
                            pan_sharpen_rgb=pan_sharpen_rgb,
                            sky_correction=sky_correction,
                            dem_pressure=dem_pressure, pressure=pressure,
                            force_band=force_band,
                            fixed_aot550=fixed_aot550, fixed_lut=fixed_lut,
                            dark_spectrum_full_scene=dark_spectrum_full_scene)
        return()

    ## PlanetScope/RapidEye
    if data_type == 'planet':
        luts_ = []
        for l in luts:
            luts_.append('-'.join(l.split('-')[0:-1]))

        ac.planetscope.planetscope_ac(bundle, output=output, limit=limit,
                            ancillary_data=ancillary_data,
                            luts=luts_, pressure=pressure,
                            uoz_default=uoz,
                            uwv_default=uwv, ignore_sr_image=ignore_sr_image,
                            map_rgb=map_rgb, map_rgb_rhos=map_rgb,
                            sky_correction=sky_correction)
        return()



    ## Worldview (2)
    if data_type == 'worldview':
        luts_ = []
        for l in luts:
            luts_.append('-'.join(l.split('-')[0:-1]))

        ac.worldview.worldview_ac(bundle, output=output, limit=limit,
                            ancillary_data=ancillary_data,
                            luts=luts_, pressure=pressure,
                            uoz_default=uoz, uwv_default=uwv,
                            map_rgb=map_rgb, map_rgb_rhos=map_rgb,
                            sky_correction=sky_correction)
        return()
