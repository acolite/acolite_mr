## part of ACOLITE_MR - does DSF atmospheric correction for Worldview
## tested on a couple of WV2 scenes
## QV 2020-03-14
## modifications: QV 2020-03-16 converted to function, added ancillary data, sky correction
##                              added full tile geolocation
##                QV 2020-07-26 added WV3 and swir bundle
##                QV 2020-09-08 added SWIR reprojection (for single tile images?)
##                QV 2020-10-26 added option for separate SWIR processing
##                QV 2020-11-24 added fixes for WorldView2 "L1B" data

def worldview_ac(bundle,
                 swir_bundle=None,
                 output=None, limit=None,
                            ancillary_data=False,
                            #luts=['ACOLITE-LUT-202010-MOD1', 'ACOLITE-LUT-202010-MOD2'],
                            luts=['PONDER-LUT-201704-MOD1', 'PONDER-LUT-201704-MOD2'],
                            uoz_default=0.3, uwv_default=1.5, pressure = None,
                            map_rgb=True, map_rgb_rhos=True,
                            gas_transmittance = True, sky_correction=True):
    import os, glob, datetime
    import numpy as np
    import acolite as ac
    from scipy.interpolate import interp2d


    ## parse the metadata
    metafile = glob.glob('{}/{}'.format(bundle,'*.XML'))
    if len(metafile)>0:
        metafile = metafile[0]
        metadata = ac.worldview.parse_metadata(metafile)
    else:
        print('No metadata found for {}'.format(bundle))
        return()

    ## parse the metadata
    if swir_bundle is not None:
        swir_metafile = glob.glob('{}/{}'.format(swir_bundle,'*.XML'))
        if len(swir_metafile)>0:
            swir_metafile = swir_metafile[0]
            swir_metadata = ac.worldview.parse_metadata(swir_metafile)
        else:
            print('No metadata found for {}'.format(swir_bundle))
            return()

    print('{} - Started processing {}'.format(datetime.datetime.now().isoformat()[0:19], bundle))

    ## import luts
    lutd=ac.aerlut.import_luts(base_luts=luts)
    luts = list(lutd.keys())
    print('{} - Imported LUTs'.format(datetime.datetime.now().isoformat()[0:19]))

    ## make output directory
    if output is None:
        odir = os.path.dirname(bundle)
    else:
        odir = '{}'.format(output)
    if not os.path.exists(odir): os.makedirs(odir)

    ## no subsetting at the moment
    if limit is not None:
        print('Limit subsetting not yet supported for Worldview')
        limit = None

    ## get output filename
    scn = os.path.splitext(metadata['TILE_INFO'][0]['FILENAME'])[0][-4:]
    oname = '{}_{}_{}_L2R.nc'.format(metadata['SATELLITE'], metadata['TIME'].strftime('%Y_%m_%d_%H_%M_%S'), scn)

    ofile = '{}/{}'.format(odir, oname)

    ## identify sensor and bands
    sensor = metadata['SATELLITE']
    rsr_file = '{}/RSR/{}.txt'.format(ac.config['pp_data_dir'], sensor)
    rsr, rsr_bands = ac.shared.rsr_read(file=rsr_file)

    ## get band wavelengths
    #band_tags=["BAND_C","BAND_B","BAND_G","BAND_Y","BAND_R","BAND_RE","BAND_N", "BAND_N2"]
    #band_names = [b for ib, b in enumerate(rsr_bands) if b != 'PAN']
    #band_indices = [ib+1 for ib, b in enumerate(rsr_bands) if b != 'PAN']

    band_tags = list(metadata['BAND_INFO'].keys())
    band_names = [metadata['BAND_INFO'][b]['name'] for b in band_tags]
    band_indices = [metadata['BAND_INFO'][b]['index'] for b in band_tags]

    #band_tags = list(metadata['BAND_INFO'].keys())
    #if swir_bundle is None:
    #    band_names = [b for ib, b in enumerate(rsr_bands) if (b != 'PAN') & ('SWIR' not in b)]
    #    band_indices = [ib+1 for ib, b in enumerate(rsr_bands) if (b != 'PAN') & ('SWIR' not in b)]
    #else:
    #    band_names = [b for ib, b in enumerate(rsr_bands) if (b != 'PAN')]
    #    band_indices = [ib+1 for ib, b in enumerate(rsr_bands) if (b != 'PAN')]

    if swir_bundle is not None:
        band_tags += list(swir_metadata['BAND_INFO'].keys())
        band_names +=  [swir_metadata['BAND_INFO'][b]['name'] for b in list(swir_metadata['BAND_INFO'].keys())]
        band_indices += [swir_metadata['BAND_INFO'][b]['index'] for b in list(swir_metadata['BAND_INFO'].keys())]

    print(band_names)
    print(band_indices)

    wave_range=[0.2,2.4]
    wave_step=0.001
    wave_hyper = np.linspace(wave_range[0],wave_range[1],int(((wave_range[1]-wave_range[0])/wave_step)+2))
    waved = ac.shared.rsr_convolute_dict(wave_hyper, wave_hyper, rsr)
    wavelengths = [waved[b]*1000 for b in band_names]
    band_wavelengths = ['{0:.0f}'.format(i) for i in wavelengths]

    ## get observation geometry
    raa = abs(float(metadata['MEANSUNAZ']) - float(metadata['MEANSATAZ']))
    while raa >= 180.: raa -= 180.
    sza = 90. - float(metadata['MEANSUNEL'])
    vza = 90. - float(metadata['MEANSATEL'])

    ## LUT limits, LUT will otherwise give NaN values for higher zenith angles
    ## VZA will be updated in a new version of the LUT (running simulations in 202010)
    if True:
        sza = np.min((sza, 80))
        vza = np.min((vza, 56))

    metadata['se_distance'] = ac.shared.distance_se(metadata['DOY'])
    se_distance = metadata['se_distance']

    ## get pressure
    press = 1013 if pressure is None else pressure

    ## get ancillary data
    if ancillary_data:
        btag = list(metadata['BAND_INFO'].keys())[0]
        scene_lons = [metadata['BAND_INFO'][btag][k] for k in metadata['BAND_INFO'][btag] if 'LON' in k]
        scene_lats = [metadata['BAND_INFO'][btag][k] for k in metadata['BAND_INFO'][btag] if 'LAT' in k]
        pc_lon = np.nanmean(scene_lons)
        pc_lat = np.nanmean(scene_lats)
        pc_date = metadata['TIME'].strftime('%Y-%m-%d')
        pc_time=metadata['TIME'].hour + metadata['TIME'].minute/60. + metadata['TIME'].second/3600.
        pc_anc = ac.ac.ancillary.ancillary_get(pc_date, pc_lon, pc_lat, ftime=pc_time, kind='nearest')
        ## get pressure from ancillary data if not determined by user or by DEM
        if pressure is None: press = pc_anc['press']['interp']

    ## get gas transmittances
    if gas_transmittance:
        uwv = uwv_default * 1.0
        uoz = uoz_default * 1.0

        if ancillary_data:
            uoz=pc_anc['ozone']['interp']
            uwv=pc_anc['p_water']['interp']/10.

        ttg = ac.ac.gas_transmittance(sza, vza, uoz=uoz, uwv=uwv)
        ttg_wave = [w/1000. for w in ttg['wave']]
        ttg_res = {par:ac.shared.rsr_convolute_dict(ttg_wave, ttg[par], rsr) for par in ttg}

    ## add sky correction here
    if sky_correction:
        ## get rayleigh optical depth
        ip = [i for i,value in enumerate(lutd[luts[0]]['meta']['par']) if value == 'tray']
        ret = lutd[luts[0]]['rgi']((press, ip, lutd[luts[0]]['meta']['wave'], raa, vza, sza, 0.001))
        tray = ac.shared.rsr_convolute_dict(lutd[luts[0]]['meta']['wave'], ret, rsr)
        ## estimated sky reflectance
        rsky = {btag:ac.ac.rayleigh.ray_refl_onlysky(waved[btag],
                                                           sza * np.pi / 180., vza * np.pi / 180.,
                                                           0, raa * np.pi / 180.,
                                                           Patm=press, tau_ray=tray[btag]) for btag in tray}

    ## get list of tiles in this bundle
    tiles=None
    if tiles is None:
        tiles=[]
        for tile_mdata in metadata['TILE_INFO']:
            tiles.append(tile_mdata['FILENAME'].split('_')[1].split('-')[0])

    ## get dimensions and dark spectrum from each tile
    tiles_dims = {}
    tiles_dims_swir = {}
    rhod={}
    for ti, tile in enumerate(tiles):
        for tile_mdata in metadata['TILE_INFO']:
            if tile in tile_mdata['FILENAME']:
                file = '{}/{}'.format(bundle,tile_mdata['FILENAME'])

                ## check if the files were named .TIF instead of .TIFF
                if not os.path.exists(file):
                    file = file.replace('.TIFF', '.TIF')

                if not os.path.exists(file):
                    tiles_dims[tile] = (0,0)
                    continue

                ## get tile from wv3 swir bundle if provided
                swir_file=None
                if swir_bundle is not None:
                    for tile_mdata_swir in swir_metadata['TILE_INFO']:
                        if tile in tile_mdata_swir['FILENAME']:
                            swir_file = '{}/{}'.format(swir_bundle,tile_mdata_swir['FILENAME'])
                    if not os.path.exists(swir_file):
                        continue
                ## end swir bundle

                for b,band in enumerate(band_names):
                    if 'SWIR' not in band:
                        d=ac.worldview.get_rtoa(file, band_indices[b], band_tags[b], metadata,
                                                sun_zenith=sza, se_distance=se_distance)
                        if tile not in tiles_dims: tiles_dims[tile] = d.shape
                    else:
                        if swir_file is None:
                            swir_file='{}'.format(file)
                            swir_metadata = metadata.copy()

                        d=ac.worldview.get_rtoa(swir_file, band_indices[b], band_tags[b], swir_metadata,
                                                sun_zenith=sza, se_distance=se_distance)
                        if tile not in tiles_dims_swir: tiles_dims_swir[tile] = d.shape

                    print(band, d.shape)

                    npix=1000
                    tmp = d.ravel()
                    tmp.sort()
                    x = np.where(np.isfinite(tmp[0:npix]))[0]
                    y = [tmp[0:npix][i] for i in x]
                    dp = y

                    ## add to rho dark dict
                    if band_names[b] not in rhod:
                        rhod[band_names[b]] = {'rhod':dp, 'raa':raa,'vza':vza, 'sza':sza,
                                               'wave': waved[band_names[b]]*1000.,
                                               'rsky': 0., 'tt_gas': 1.}
                        if gas_transmittance:
                            rhod[band_names[b]]['tt_gas'] = ttg_res['tt_gas'][band_names[b]]
                        if sky_correction:
                            rhod[band_names[b]]['rsky'] = rsky[band_names[b]]
                    else:
                        rhod[band_names[b]]['rhod']+=dp
                    d = None


    ## sort the rhod in each band
    ## fit linear regression and get offset as estimated rhod
    for b,band in enumerate(band_names):
        rhod[band_names[b]]['rhod'] = np.sort(rhod[band_names[b]]['rhod'])
        tmp = rhod[band_names[b]]['rhod']
        x = np.where(np.isfinite(tmp[0:npix]))[0]
        y = [tmp[0:npix][i] for i in x]

        if gas_transmittance: y/= ttg_res['tt_gas'][band_names[b]]
        if sky_correction: y-= rsky[band_names[b]]

        my, by, ry, smy, sby = ac.shared.regression.lsqfity(x, y)
        rhod[band_names[b]]['rhod_tmp'] = tmp
        rhod[band_names[b]]['rhod'] = by

    ## fit model
    if pressure is None:
        pressure = 1013.25
    res = ac.ac.select_model2(rhod, sensor, pressure=pressure, rhod_tgas_cutoff=0.9, lutd=lutd)
    sel_model = res['lutid']
    sel_model_band_pair = [band_names[s] for s in (res['sel_idx'],res['sel_idx2'])]
    mod = {'1':'C', '2':'M', '3':'U'}[sel_model[-1]]
    print('{} - Fitted model {}, band pair {}, taua 550: {:.3f}'.format(datetime.datetime.now().isoformat()[0:19],
                                                                        mod, ':'.join(sel_model_band_pair), res['taua']), end='\n')

    ## using new sky reflectance correction
    if False:
        res2 = ac.ac.rhod_fit_model(raa, vza, sza, pressure = pressure,
                        rhod = {b:rhod[b]['rhod_tmp'][0] for b in rhod}, sensor=sensor, resample=True,luts=luts)

        sel_mod = res2['mod_sel']
        tau550 = res2[sel_mod]['tau_fit']
        if np.isnan(tau550): tau550 = 5

        mod = {'1':'C', '2':'M', '3':'U'}[sel_mod[-1]]

        rsky_new = {b: res2[sel_mod]['rhot_fit_rs'][b]-res2[sel_mod]['romix_fit_rs'][b] for b in res2[sel_mod]['rhot_fit_rs']}
        print('model:{} band:{} aot={:.3f}'.format(sel_mod,res2[sel_mod]['tau_band'],res2[sel_mod]['tau_fit']))

    ## for debugging
    if False:
        return(rhod, res)

    ## get A/C parameters per band
    pars = ['romix','dtott','utott','astot']
    rgi = res['rgi']
    ac_pars = {}
    for ip, par in enumerate(pars):
        ip = [i for i,value in enumerate(res['lut_meta']['par']) if value == par]
        if len(ip) == 1: ip = ip[0]
        else: continue
        ret = rgi((press, ip, res['lut_meta']['wave'], raa, vza, sza, res['taua']))
        ac_pars[par] = ac.shared.rsr_convolute_dict(res['lut_meta']['wave'], ret, rsr)
        ret = None

    ## write results to output file
    new=True
    global_dims = [int(metadata['NUMROWS']),int(metadata['NUMCOLUMNS'])]

    ## write lat/lon
    print('{} - Writing lat/lon'.format(datetime.datetime.now().isoformat()[0:19]))
    pcol = [0, global_dims[1], global_dims[1], 0]
    prow = [0, 0, global_dims[0], global_dims[0]]
    plon = []
    plat = []
    band_tag = list(metadata['BAND_INFO'].keys())[0]
    for bk in ['UL', 'UR', 'LR', 'LL']:
            k = '{}{}'.format(bk, 'LON')
            plon.append(metadata['BAND_INFO'][band_tag][k])
            k = '{}{}'.format(bk, 'LAT')
            plat.append(metadata['BAND_INFO'][band_tag][k])

    ## set up interpolator
    zlon = interp2d(pcol, prow, plon, kind='linear')
    zlat = interp2d(pcol, prow, plat, kind='linear')
    x = np.arange(1, 1+global_dims[1], 1)
    y = np.arange(1, 1+global_dims[0], 1)
    ac.output.nc_write(ofile, 'lat', zlat(x, y), global_dims=global_dims, new=new)
    ac.output.nc_write(ofile, 'lon', zlon(x, y))
    new = False
    ## end write lat/lon

    ## get SWIR lat/lon
    if swir_bundle is not None:
        swir_dims = [int(swir_metadata['NUMROWS']),int(swir_metadata['NUMCOLUMNS'])]
        pcol = [0, swir_dims[1], swir_dims[1], 0]
        prow = [0, 0, swir_dims[0], swir_dims[0]]
        plon = []
        plat = []
        for bk in ['UL', 'UR', 'LR', 'LL']:
                k = '{}{}'.format(bk, 'LON')
                plon.append(swir_metadata['BAND_INFO']['BAND_S1'][k])
                k = '{}{}'.format(bk, 'LAT')
                plat.append(swir_metadata['BAND_INFO']['BAND_S1'][k])

        ## set up interpolator
        swir_zlon = interp2d(pcol, prow, plon, kind='linear')
        swir_zlat = interp2d(pcol, prow, plat, kind='linear')
        swir_x = np.arange(1, 1+swir_dims[1], 1)
        swir_y = np.arange(1, 1+swir_dims[0], 1)
        swir_lat = swir_zlat(swir_x,swir_y)
        swir_lon = swir_zlon(swir_x,swir_y)
        ## also read vnir lat lon for reprojection
        vnir_lat = ac.shared.nc_data(ofile, 'lat')
        vnir_lon = ac.shared.nc_data(ofile, 'lon')
    ## end get SWIR lat/lon

    for ti, tile in enumerate(tiles):
        for tile_mdata in metadata['TILE_INFO']:
            if tile in tile_mdata['FILENAME']:
                file = '{}/{}'.format(bundle,tile_mdata['FILENAME'])

                ## check if the files were named .TIF instead of .TIFF
                if not os.path.exists(file):
                    file = file.replace('.TIFF', '.TIF')

                if not os.path.exists(file): continue

                ## get tile from wv3 swir bundle if provided
                if swir_bundle is not None:
                    for tile_mdata_swir in swir_metadata['TILE_INFO']:
                        if tile in tile_mdata_swir['FILENAME']:
                            swir_file = '{}/{}'.format(swir_bundle,tile_mdata_swir['FILENAME'])
                    if not os.path.exists(swir_file):
                        continue
                ## end swir bundle

                ## get tile offset
                offset = [int(tile_mdata['ULCOLOFFSET']), int(tile_mdata['ULROWOFFSET'])]
                print('Processing tile', tile, offset)

                ## function to nn resize an array to a given shape
                def nn_resize(arr, size):
                    return(np.asarray([[ arr[int(arr.shape[0] * r / size[0])][int(arr.shape[1] * c / size[1])]
                             for c in range(size[1])] for r in range(size[0])]))

                ## run through spectral bands
                for b,band in enumerate(band_names):
                    ds_att = {'wavelength':wavelengths[b], 'tt_gas': 1., 'rsky': 0}
                    wave  = band_wavelengths[b]

                    if sky_correction:
                        ds_att['rsky'] = rsky[band_names[b]]
                    if gas_transmittance:
                        for key in ttg_res: ds_att[key]=ttg_res[key][band_names[b]]

                    for key in ac_pars: ds_att[key]=ac_pars[key][band_names[b]]
                    if not os.path.exists(file): continue

                    if 'SWIR' not in band:
                        d=ac.worldview.get_rtoa(file, band_indices[b], band_tags[b], metadata,
                                                sun_zenith=sza, se_distance=se_distance)
                    else:
                        d=ac.worldview.get_rtoa(swir_file, band_indices[b], band_tags[b], swir_metadata,
                                                sun_zenith=sza, se_distance=se_distance)
                        ### next line needs to be reprojection to VNIR data extent, not nn interpolation
                        #d = nn_resize(d, tiles_dims[tile])
                        ## reproject SWIR data to VNIR data extent
                        if swir_file != file:
                            d = ac.shared.reproject2(d, swir_lon, swir_lat, vnir_lon, vnir_lat, fill=np.nan)
                    print(band, d.shape)


                    ac.output.nc_write(ofile, 'rhot_{}'.format(wave), d, dataset_attributes=ds_att,
                                           offset=offset, global_dims=global_dims, new=new)
                    new = False

                    if gas_transmittance: d/= ds_att['tt_gas']
                    if sky_correction: d-=  ds_att['rsky']

                    rhos_data = d - ds_att['romix']
                    rhos_data = (rhos_data) / ((ds_att['utott']*ds_att['dtott']) + ds_att['astot'] * rhos_data)

                    ac.output.nc_write(ofile, 'rhos_{}'.format(wave), rhos_data, dataset_attributes=ds_att,
                                      offset=offset, global_dims=global_dims, new=new)
                    d = None
                    rhos_data = None

    print('{} - Finished writing {}'.format(datetime.datetime.now().isoformat()[0:19], ofile))
    return(ofile)
