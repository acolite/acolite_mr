## def pleiades_ac
## performs Pleiades atmospheric correction
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07
## modifications: QV 2017-01-18 converted to function, added pan sharpened RGB
##                QV 2017-01-23 changed data reading to dict
##                QV 2017-01-25 added full scene dark spectrum option and output additional attributes to netcdf
##                QV 2017-03-23 added support for multiple scenes in the file
##                QV 2017-05-03 updated default luts, added absolute_pixel_list option
##                QV 2017-10-18 added NCEP+TOAST ancillary data
##                              added TOA gas_transmittance option
##                              added pressure support to LUT
##                              added scene center DEM based pressure estimate
##                QV 2017-10-19 added force_band option
##                QV 2017-11-21 added model_selection and rdark_list_selection keywords
##                              added sky correction option
##                              uses new model selection function
##                QV 2017-12-06 added gains option
##                QV 2018-01-16 updated to new nc_write function without metadata
##                QV 2018-05-24 added return of tau and model for return_rho_s keyword
##                QV 2019-02-25 changed import name
##                QV 2019-03-13 removed scene name in directory and some old plotting stuff
##                QV 2020-10-18 added support for multi tile images

def pleiades_ac(file, output=None, sub=None, limit=None, region=None,
                write_netcdf=True,
                write_netcdf_geo=True,
                write_netcdf_rhot=True,
                write_netcdf_rhos=True,
                write_netcdf_rhorc=True,
                map_rgb = True,
                map_rgb_rhos = True,
                #plot_dark_spectrum=True,
                rgb_autorange = False,
                rgb_percentiles = [1,90],
                rgb_range_min = [0.0,0.0,0.0],
                rgb_range_max = [0.15,0.15,0.15],
                pan_sharpen_rgb = False,
                pan_sharpen_rgb_netcdf=False,
                dark_spectrum_option='dark_list',
                #dark_spectrum_smooth = False, ## not used
                dark_spectrum_full_scene=False,
                pixel_idx=10,
                perc_idx=1,
                percentiles = [0,0.1,1,5,10,25,50,75,90,95,99,99.9,100],
                ## options for choosing model
                model_selection='min_tau',
                ## if dark_spectrum_option is 'dark_list' the next one selects how this list is treated in select_model
                rdark_list_selection='intercept',
                ## default luts
                luts=['PONDER-LUT-201704-MOD1-1013mb', 'PONDER-LUT-201704-MOD2-1013mb'],
                fixed_aot550=None,
                fixed_lut='PONDER-LUT-201704-MOD2-1013mb',
                force_band = None,
                bestfit='bands',
                bestfit_bands=["Red","NIR"],
                pixel_range_min=0,
                pixel_range_max=1000,
                ## pressure
                lut_pressure = True,
                dem_pressure = False,
                dem_pressure_percentile = 25,
                pressure = None,
                ## apply gas corrections at TOA
                gas_transmittance = True,
                uoz = 0.3,
                uwv = 1.5,
                wvlut = '201710C',
                ## do sky glint correction
                sky_correction = False,
                sky_correction_option = 'all',
                ## use ancillary data for gas transmittances rather than defaults
                ancillary_data = True,
                ## only return rho dark and metadata
                return_rho_dark = False,
                ## return rho_s, lat, lon and metadata
                return_rho_s = False,
                gains=None
                ):

    import os, time
    import acolite as ac
    import numpy as np

    ## check if we are dealing with Pléiades image bundle
    ifiles,mfiles,pifiles,pmfiles = ac.pleiades.image_test(file, listpan=True)

    mfiles_set = set(mfiles)
    if len(mfiles_set)>1:
        print('Multiple metadata files found')
        return(1)

    ## read metadata only once
    for mfile in mfiles_set:
        metadata = ac.pleiades.parse_metadata(mfile)
        sensor = metadata['SENSOR']
        bands = metadata['BANDS']
        band_names = metadata['BAND_NAMES']
        band_indices = metadata['BAND_INDICES']
        wavelengths = metadata['WAVELENGTHS']
        band_wavelengths = metadata['BAND_WAVELENGTHS']
        rgb_bands = ['Red','Green','Blue','NIR']
        tile_dims = metadata['tiles_nrows'], metadata['tiles_ncols']
        data = {}

    ## make output directory
    scname = os.path.basename(file)
    if output is None:
        output='{}'.format(os.path.dirname(file))
    else:
        output='{}'.format(output)

    ## add the region to the output directory structure
    odir = '{}'.format(output)
    if (region is not None):
        odir = '{}/{}'.format(odir, region)
    if not os.path.exists(odir): os.makedirs(odir)

    ## make output filename base
    obase = '{}{}_{}'.format(metadata['INSTRUMENT'],metadata['INSTRUMENT_INDEX'],metadata['TIME'].strftime('%Y_%m_%d_%H_%M_%S'))
    metadata['OBASE'] = obase
    metadata['ODIR'] = odir
    metadata['OUTPUT_BASE'] = '{}/{}'.format(odir,obase)

    ## get pressure from dem (scene center)
    if (dem_pressure) & (pressure == None) & (lut_pressure):
            pc_lon=metadata['VERTICES']['C']['LON']
            pc_lat=metadata['VERTICES']['C']['LAT']
            ## interpolate DEM
            dem = ac.dem.hgt_lonlat(pc_lon, pc_lat)
            pressure = (ac.ac.pressure_elevation(dem, ratio=False))[0]

    ## get NCEP & TOAST ancillary data
    if ancillary_data:
        pc_lon=metadata['VERTICES']['C']['LON']
        pc_lat=metadata['VERTICES']['C']['LAT']
        pc_date = metadata['TIME'].strftime('%Y-%m-%d')
        pc_time=metadata['TIME'].hour + metadata['TIME'].minute/60. + metadata['TIME'].second/3600.
        pc_anc = ac.ac.ancillary.ancillary_get(pc_date, pc_lon, pc_lat, ftime=pc_time, kind='nearest')
        ## get pressure from ancillary data if not determined by user or by DEM
        if (pressure == None) & (lut_pressure): pressure = pc_anc['press']['interp']

    ## get gas transmittances
    if gas_transmittance:
        if ancillary_data:
            uoz=pc_anc['ozone']['interp']
            uwv=pc_anc['p_water']['interp']/10.
        tt_oz = ac.ac.o3_transmittance(sensor, metadata, uoz=uoz)
        tt_wv = ac.ac.wvlut_interp(metadata['THS'], metadata['THV'], uwv=uwv, sensor=sensor, config=wvlut)
        tt_gas = {btag: tt_oz[btag] * tt_wv[btag] for btag in tt_oz.keys()}

    ## check if we need to find crop position
    if limit is not None:
        if len(limit) is 4:
            ncols = int(metadata['NCOLS'])
            nrows = int(metadata['NROWS'])
            sub = ac.pleiades.geo.crop(metadata, limit)
            if (sub[2] <= 50) or (sub[3] <= 50):
                print("Region not in image.")
                print("Crop dimensions: {}x{}.".format(sub[2],sub[3]))
                #return 1
            if (sub[0] <= 0): sub[0] = 0
            if (sub[1] <= 0): sub[1] = 0
            if (sub[0]+sub[2] >= ncols):
                sub[2] -= ((sub[0]+sub[2])-ncols)
            if (sub[1]+sub[3] >= nrows):
                sub[3] -= ((sub[1]+sub[3])-nrows)

    if sub is not None: cropname = '_'.join([str(i) for i in sub])


    if sub is None:
        dims = int(metadata['NROWS']), int(metadata['NCOLS'])
    else:
        dims = sub[3], sub[2]


    ## run through image tiles
    t = time.process_time()
    tile_id = 0
    for r_tile in range(metadata['ntiles_R']):
        for c_tile in range(metadata['ntiles_C']):
            tile_id+=1
            tile_name = 'R{}C{}'.format(r_tile+1,c_tile+1)
            print(tile_id, tile_name)

            ## tile offsets
            tile_row_off = r_tile * tile_dims[0]
            tile_col_off = c_tile * tile_dims[1]

            ifile = None
            for it,tfile in enumerate(ifiles):
                if tile_name not in tfile: continue
                ifile=ifiles[it]
                mfile=mfiles[it]
                pifile=pifiles[it]
                pmfile=pmfiles[it]

            print('Reading: {}'.format(ifile))

            ## read in TOA reflectances
            for bidx,band in enumerate(band_names):
                data_b = ac.pleiades.get_rtoa(ifile, band_indices[bidx]+1, bands[bidx], metadata, sub=sub, gains=gains)
                cur_shape = data_b.shape
                print(band, cur_shape)
                if tile_id == 1:
                    data[band] = np.zeros(dims)
                data[band][tile_row_off:tile_row_off+cur_shape[0],
                           tile_col_off:tile_col_off+cur_shape[1]] = data_b

            ## read pan data
            if (pifile is not None) & (pan_sharpen_rgb is True):
                if sub is None:
                    print('Pan processing currently only with subsetting.')
                    return(1)

                if tile_id == 1:
                    pan_data = np.zeros((dims[0]*4, dims[1]*4))


                pan_sub = [int(i)*4 for i in sub]
                data_pan = ac.pleiades.get_rtoa(pifile, 1, 'P', pan_metadata,sub=pan_sub)
                cur_shape = data_pan.shape

                pan_data[tile_row_off*4:tile_row_off*4+cur_shape[0],
                         tile_col_off*4:tile_col_off*4+cur_shape[1]] = data_pan
            elapsed = time.process_time()-t
    print('Time elapsed reading data: ',elapsed)
    ## finished reading data

    ## keep scene edges
    mask = data[band_names[0]] == np.nan

    ## get dark spectrum
    #if dark_spectrum_full_scene is False:
    #    rtoa_dict, perc_all = ac.ac.get_dark_spectrum(data, option=dark_spectrum_option,
    #                                                      percentiles=percentiles, perc_idx=perc_idx,
    #                                                      pixel_idx=pixel_idx,pixel_range_min=pixel_range_min,pixel_range_max=pixel_range_max)

    rtoa_dict, perc_all = ac.ac.get_dark_spectrum(data, option=dark_spectrum_option,
                                                          percentiles=percentiles, perc_idx=perc_idx,
                                                          pixel_idx=pixel_idx,pixel_range_min=pixel_range_min,pixel_range_max=pixel_range_max)


    ### TOA gas correction
    if gas_transmittance:
        rtoa_dict = {band:rtoa_dict[band]/tt_gas[band] for band in rtoa_dict}
        data = {band:data[band]/tt_gas[band] for band in data.keys()}

    ## Sky reflectance correction
    if sky_correction:
        rsky = ac.ac.toa_rsky(metadata, pressure=pressure)
        if sky_correction_option == 'all':
            rtoa_dict = {band:rtoa_dict[band]-rsky[band] for band in rtoa_dict.keys()}
            data = {band:data[band]-rsky[band] for band in data.keys()}

        if sky_correction_option == 'auto':
            print('Sky correction option {} not configured'.format(sky_correction_option))


    ## return only rho_dark and metadata, no further processing
    if return_rho_dark:
        if return_rho_s: return_rho_dark=False
        else: return(rtoa_dict,metadata)

    (ratm_s,rorayl_s,dtotr_s,utotr_s,dtott_s,utott_s,astot_s, tau550),\
    (bands_sorted, tau550_all_bands, dark_idx, sel_rmsd, rdark_sel, pixel_idx), \
    (sel_model_lut, sel_model_lut_meta) = ac.ac.select_model(metadata, rtoa_dict, luts=luts,
                                                           bestfit=bestfit, bestfit_bands=bestfit_bands,
                                                           model_selection=model_selection,
                                                           rdark_list_selection=rdark_list_selection,
                                                           pressure=pressure,force_band=force_band)

    if fixed_aot550 is not None:
        (ratm_s,rorayl_s,dtotr_s,utotr_s,dtott_s,utott_s,astot_s, tau550),\
        sel_model_lut, sel_model_lut_meta = ac.ac.fixed_model(metadata, fixed_aot550, lut=fixed_lut, pressure=pressure)

    print('tau {}, band {}, model {}, pixel {}'.format(tau550, dark_idx, sel_model_lut_meta['aermod'][0], pixel_idx))



    ## convert dicts to lists
    ratm = [ratm_s[band] for band in bands_sorted if band in rdark_sel]
    rorayl = [rorayl_s[band] for band in bands_sorted if band in rdark_sel]
    dtotr = [dtotr_s[band] for band in bands_sorted if band in rdark_sel]
    utotr = [utotr_s[band] for band in bands_sorted if band in rdark_sel]
    dtott = [dtott_s[band] for band in bands_sorted if band in rdark_sel]
    utott = [utott_s[band] for band in bands_sorted if band in rdark_sel]
    astot = [astot_s[band] for band in bands_sorted if band in rdark_sel]

    if 'aermod' in sel_model_lut_meta.keys():
        if sel_model_lut_meta['aermod'][0] == "1": model_char = 'C'
        if sel_model_lut_meta['aermod'][0] == "2": model_char = 'M'
        if sel_model_lut_meta['aermod'][0] == "3": model_char = 'U'
    else:
        model_char = '4C'
        model_char = '4C: {}/{}/{}/{}'.format(sel_model_lut_meta['mod1'],sel_model_lut_meta['mod2'],sel_model_lut_meta['mod3'],sel_model_lut_meta['mod4'])
    ac_tuple = (ratm,rorayl,dtotr,utotr,dtott,utott,astot)

    if (dark_spectrum_full_scene is False) and (sub != None): ds_origin = 'sub scene'
    else: ds_origin = 'full scene'

    ## set up attributes
    attributes={'inputfile':file,'ths':metadata['THS'], 'thv':metadata['THV'], 'azi':metadata['AZI'], 'se_distance':metadata['SE_DISTANCE']}
    attributes['ds_origin']=ds_origin
    attributes['ds_option']=dark_spectrum_option
    attributes['ds_pixel_idx']=pixel_idx
    attributes['ds_percentile']=percentiles[perc_idx]
    attributes['ds_bestfit']=bestfit

    attributes['ac_model']=sel_model_lut_meta['base']#[0]
    attributes['ac_model_char']=model_char
    attributes['ac_band']=dark_idx
    attributes['ac_aot550']=tau550#tau550_all_bands[dark_idx]
    attributes['ac_rmsd']=sel_rmsd

    if gas_transmittance is True:
        attributes["uoz"] = uoz
        attributes["uwv"] = uwv
        for band in tt_gas.keys():
            attributes['{}_t_gas'.format(band)] = tt_gas[band]

    ## add ancillary data information
    if pressure is not None: attributes["pressure"] = pressure
    if ancillary_data is True:
        for k in pc_anc.keys():
            if type(pc_anc[k]) is dict:
                if 'series' in pc_anc[k].keys(): attributes['anc_{}_series'.format(k)] = pc_anc[k]['series']
                if 'interp' in pc_anc[k].keys(): attributes['anc_{}'.format(k)] = pc_anc[k]['interp']
            else:
                attributes['anc_{}'.format(k)] = pc_anc[k]


    print('model:{}, band:{}, aot={:.3f}'.format(attributes['ac_model_char'],attributes['ac_band'],attributes['ac_aot550']))

    if True:
        ds_plot = '{}/{}_{}.{}'.format(odir,obase, 'dark_spectrum','png')
        #band_names = bands
        data_type = 'NetCDF'
        #waves = [metadata['bands'][b]['wave_name'] for b in bands]
        for ib,band in enumerate(bands):
            bk='{}-{}'.format(band, 'wave_name')
            metadata[bk] = metadata['BAND_WAVELENGTHS'][ib]

        waves = [metadata['{}-{}'.format(b,'wave_name')] for b in bands]
        dsf_spectrum_option = 'fixed'

        #metadata['SENSOR']='PlanetScope'
        metadata['SATELLITE_SENSOR'] = '{}_{}'.format(metadata['SATELLITE'],metadata['SENSOR'])
        print(metadata['SENSOR'], metadata['SATELLITE'])

        ac.plotting.plot_dark_spectrum(metadata, ds_plot, waves, ratm_s, rorayl_s, rdark_sel, dark_idx, tau550,sel_model_lut_meta, xlim=(450,900))

    for li,lut in enumerate(luts):
        attributes['ac_lut{}'.format(li)] = lut

    #mkeys = ['GEOMETRY', 'RESAMPLING_SPACING', 'MISSION_INDEX', 'NBANDS', 'VERTICES', 'INSTRUMENT_INDEX', 'MISSION', 'NROWS', 'IMAGING_DATE', 'ALPHA_CHANNEL', 'NCOLS', 'INSTRUMENT', 'BAND_MODE', 'DOY', 'GREEN_CHANNEL', 'EXTENT_TYPE', 'RED_CHANNEL', 'SATURATED', 'TIME', 'BLUE_CHANNEL', 'NODATA', 'IMAGING_TIME', 'BAND_INFO']
    mkeys = ['MISSION','MISSION_INDEX','NBANDS','NCOLS', 'NROWS',  'IMAGING_DATE','IMAGING_TIME', 'DOY','RESAMPLING_SPACING']
    for key in mkeys: attributes[key]=metadata[key]



    ## calculate surface reflectance
    rhos_data={}
    for bidx,band in enumerate(bands_sorted):
        if band == 'Pan': continue
        rhos_data[band] = ac.rtoa_to_rhos(data[band], ratm_s[band], utott_s[band], dtott_s[band], astot_s[band], tt_gas = 1.) ## note tt_gas has been applied to data dict already

    if return_rho_s:
        if sub is None:
             lon, lat = ac.pleiades.geo.ll(metadata)
        else:
             lon, lat = ac.pleiades.geo.ll(metadata, sub=sub)
        #return(rtoa_dict,metadata,rhos_data,lat,lon,tau550,sel_model_lut_meta)

    ## map rgb images
    ## keep image 3d matrix for further plotting (if needed)
    rgb_image = None
    if map_rgb:
            rgb_rhot_dir = '{}'.format(odir)
            if not os.path.exists(rgb_rhot_dir): os.makedirs(rgb_rhot_dir)
            rgb_rhot_file = '{}/{}_{}.{}'.format(rgb_rhot_dir,obase,'RGB_RHOT','png')
            rgb_image = ac.output.write_rgb(rgb_rhot_file, data[rgb_bands[0]], data[rgb_bands[1]], data[rgb_bands[2]],
                            rgb_autorange=rgb_autorange, rgb_percentiles=rgb_percentiles,
                            rgb_range_min=rgb_range_min,rgb_range_max=rgb_range_max, return_image=True)

    if map_rgb_rhos:
            rgb_rhos_dir = '{}'.format(odir)
            if not os.path.exists(rgb_rhos_dir): os.makedirs(rgb_rhos_dir)
            rgb_rhos_file = '{}/{}_{}.{}'.format(rgb_rhos_dir,obase,'RGB_RHOS','png')
            rgb_image = ac.output.write_rgb(rgb_rhos_file, rhos_data[rgb_bands[0]], rhos_data[rgb_bands[1]], rhos_data[rgb_bands[2]],
                            rgb_autorange=rgb_autorange, rgb_percentiles=rgb_percentiles,
                            rgb_range_min=rgb_range_min,rgb_range_max=rgb_range_max, return_image=True)

    ## make pan sharpened RGBs
    if pan_sharpen_rgb:
        import scipy.ndimage
        if map_rgb:
            red_ps = scipy.ndimage.zoom(data[rgb_bands[0]], 4, order=0)
            green_ps = scipy.ndimage.zoom(data[rgb_bands[1]], 4, order=0)
            blue_ps = scipy.ndimage.zoom(data[rgb_bands[2]], 4, order=0)
            nir_ps = scipy.ndimage.zoom(data[rgb_bands[3]], 4, order=0)
            pf = 1./ (((red_ps+green_ps+blue_ps+nir_ps)/4.)/pan)
            red_ps*=pf
            green_ps*=pf
            blue_ps*=pf
            rgb_rhot_dir = '{}'.format(odir)
            if not os.path.exists(rgb_rhot_dir): os.makedirs(rgb_rhot_dir)
            rgb_rhot_file = '{}/{}_{}.{}'.format(rgb_rhot_dir,obase,'RGB_RHOT_PAN','png')
            ac.output.write_rgb(rgb_rhot_file, red_ps, green_ps, blue_ps,
                            rgb_autorange=rgb_autorange, rgb_percentiles=rgb_percentiles,
                            rgb_range_min=rgb_range_min,rgb_range_max=rgb_range_max, return_image=False)
        if map_rgb_rhos:
            red_ps = scipy.ndimage.zoom(data[rgb_bands[0]], 4, order=0)
            green_ps = scipy.ndimage.zoom(data[rgb_bands[1]], 4, order=0)
            blue_ps = scipy.ndimage.zoom(data[rgb_bands[2]], 4, order=0)
            nir_ps = scipy.ndimage.zoom(data[rgb_bands[3]], 4, order=0)
            pf = 1./ (((red_ps+green_ps+blue_ps+nir_ps)/4.)/pan)

            red_ps = scipy.ndimage.zoom(rhos_data[rgb_bands[0]], 4, order=0)
            green_ps = scipy.ndimage.zoom(rhos_data[rgb_bands[1]], 4, order=0)
            blue_ps = scipy.ndimage.zoom(rhos_data[rgb_bands[2]], 4, order=0)
            red_ps*=pf
            green_ps*=pf
            blue_ps*=pf

            rgb_rhos_dir = '{}'.format(odir)
            if not os.path.exists(rgb_rhos_dir): os.makedirs(rgb_rhos_dir)
            rgb_rhos_file = '{}/{}_{}.{}'.format(rgb_rhos_dir,obase,'RGB_RHOS_PAN','png')
            ac.output.write_rgb(rgb_rhos_file, red_ps, green_ps, blue_ps,
                            rgb_autorange=rgb_autorange, rgb_percentiles=rgb_percentiles,
                            rgb_range_min=rgb_range_min,rgb_range_max=rgb_range_max, return_image=False)

    ## output netcdf
    if write_netcdf:
        #netcdf_dir = '{}/{}'.format(odir,'BestFit')
        netcdf_dir = '{}'.format(odir)
        if not os.path.exists(netcdf_dir): os.makedirs(netcdf_dir)
        netcdf_file = '{}/{}.{}'.format(netcdf_dir,obase,'nc')
        new = True

        netcdf_file_pan = '{}/{}_pan.{}'.format(netcdf_dir,obase,'nc')
        new_pan = True

        if write_netcdf_geo:
            if sub is None:
                    lon, lat = ac.pleiades.geo.ll(metadata)
            else:
                    lon, lat = ac.pleiades.geo.ll(metadata, sub=sub)
            ac.output.nc_write(netcdf_file, 'lon', lon, new=new, attributes=attributes)
            ac.output.nc_write(netcdf_file, 'lat', lat)
            new = False

        if pan_sharpen_rgb_netcdf:
            lon = scipy.ndimage.zoom(lon, 4)
            ac.output.nc_write(netcdf_file_pan, 'lon', lon, new=new_pan, attributes=attributes)
            lon = None

            lat = scipy.ndimage.zoom(lat, 4)
            ac.output.nc_write(netcdf_file_pan, 'lat', lat)
            lat = None


            ac.output.nc_write(netcdf_file_pan, 'Red', ac.shared.datascl(red_ps, dmin=rgb_range_min[0], dmax=rgb_range_max[0]))
            ac.output.nc_write(netcdf_file_pan, 'Green', ac.shared.datascl(green_ps, dmin=rgb_range_min[1], dmax=rgb_range_max[1]))
            ac.output.nc_write(netcdf_file_pan, 'Blue', ac.shared.datascl(blue_ps, dmin=rgb_range_min[2], dmax=rgb_range_max[2]))

        if write_netcdf_rhot:
            for bidx,band in enumerate(band_names):
                if band == 'Pan': continue
                ac.output.nc_write(netcdf_file, 'rhot_'+band_wavelengths[bidx], data[band], \
                                       wavelength=float(wavelengths[bidx])*1000.,new=new, attributes=attributes)
                new=False

        if write_netcdf_rhos:
            for bidx,band in enumerate(band_names):
                if band == 'Pan': continue
                ac.output.nc_write(netcdf_file, 'rhos_'+band_wavelengths[bidx], rhos_data[band], \
                                       wavelength=float(wavelengths[bidx])*1000.,new=new, attributes=attributes)
                new=False

        if write_netcdf_rhorc:
            for bidx,band in enumerate(band_names):
                if band == 'Pan': continue
                ac.output.nc_write(netcdf_file, 'rhorc_'+band_wavelengths[bidx], data[band]-rorayl[bidx], \
                                       wavelength=float(wavelengths[bidx])*1000.,new=new, attributes=attributes)
                new=False
    data = None

    print('Finished processing {}'.format(file))
