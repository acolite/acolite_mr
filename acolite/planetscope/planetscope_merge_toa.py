## def planetscope_merge_toa
## merges PlanetScope TOA tiles
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-15
## modifications: 2018-03-19 (QV) added unzipped folder delete, fixed bundle test and added replace_nan for nc writing
##                2018-04-23 (QV) added nc_compression
##                2018-05-08 (QV) changed Rapideye output name
##                2019-09-25 (QV) added sr if exists in separate file and include_planet_sr keyword
##                2019-10-01 (QV) skipping files with different UTM zone

def planetscope_merge_toa(files, output, limit=None, subname = 'merged', nc_compression=True, include_planet_sr=False):
    import os
    from ponder_processor import planetscope
    from ponder_processor.output import nc_write

    nc_l1r_new,nc_file_l1r = True, ''
    extend_limit = True

    for bundle in files:
        zipped=False
        if bundle[-4:] == '.zip':
             zipped=True
             import zipfile, shutil
             bundle_orig = '{}'.format(bundle)
             bundle,ext = os.path.splitext(bundle_orig)
             if os.path.exists(bundle) is False:
                 try:
                     print(bundle_orig)
                     zip_ref = zipfile.ZipFile(bundle_orig, 'r')
                     #zip_ref.extractall(bundle)
                     for z in zip_ref.infolist():
                         z.filename = os.path.basename(z.filename)
                         zip_ref.extract(z, bundle)
                     zip_ref.close()
                 except:
                     print("Error extracting {}".format(bundle_orig))
                     os.remove(bundle_orig)
                     return()

             if include_planet_sr:
                 bundle_sr = bundle_orig.replace('analytic', 'analytic_sr')
                 if os.path.exists(bundle_sr):
                     try:
                         print(bundle_sr)
                         zip_ref = zipfile.ZipFile(bundle_sr, 'r')
                         zip_ref.extractall(bundle)
                         zip_ref.close()
                     except:
                         print("Error extracting {}".format(bundle_sr))
                         os.remove(bundle_orig)
                         return()

        files = planetscope.bundle_test(bundle)

        sr_image_file = None
        metafile = None
        image_file = None
        process = True

        if 'metadata' in files: metafile = files['metadata']['path']
        if 'analytic' in files: image_file = files['analytic']['path']
        if 'sr' in files: sr_image_file = files['sr']['path']
        if (sr_image_file is None) & (include_planet_sr):
            print('SR image file not found')
            process = False

        if metafile == None:
            print('Could not find metadata for {}'.format(bundle))
            process = False

        if process:
            try:
                metadata = planetscope.parse_metadata(metafile)
            except:
                print('Could not parse metadata: {}'.format(metafile))
                return()

            metadata['image_file'] = image_file

            sub = None
            sensor = metadata['LUT_SENSOR']

            if sensor == 'PlanetScope_0d':
                print('Sensor {} not yet implemented'.format(sensor))
                process = False

        if process:
            if limit is not None:
                ret = planetscope.geo.get_sub(metadata, limit)
                if type(ret) == int:
                    print('Region not in scene.')
                    process = False
                else:
                    sub, p, (xrange,yrange,grid_region) = ret

        if process:
            if not os.path.exists(output): os.makedirs(output)

            if nc_l1r_new:
                zone = metadata['zone']

                bands = metadata['BANDS_ALL']
                ## add PS satellite id
                if metadata['SENSOR']=='PlanetScope':
                    obase = '{}_{}_{}_{}'.format(metadata['SENSOR'], metadata['TIME'].strftime('%Y_%m_%d_%H_%M'), 'PS{}'.format(metadata['SATELLITE_SENSOR'].split('_')[1]),subname)
                else:
                    #obase = '{}_{}_{}'.format(metadata['SATELLITE_SENSOR'], metadata['TIME'].strftime('%Y_%m_%d_%H_%M'),subname)
                    obase = '{}_{}_{}_{}'.format(metadata['SENSOR'], metadata['TIME'].strftime('%Y_%m_%d_%H_%M'), 'RE{}'.format(metadata['SATELLITE_SENSOR'].split('-')[1]),subname)
                nc_file_l1r = '{}/{}_L1R.nc'.format(output, obase)
                global_dims = (grid_region['dims'][1],grid_region['dims'][0])
                metadata['limit'] = limit
                metadata['obase'] = obase

                lon, lat = planetscope.get_ll(metadata,limit=limit,extend_limit=extend_limit)
                nc_write(nc_file_l1r, 'lon', lon, new=nc_l1r_new, attributes=metadata, nc_compression=nc_compression)
                nc_l1r_new=False
                lon = None
                nc_write(nc_file_l1r, 'lat', lat, new=nc_l1r_new, attributes=metadata, nc_compression=nc_compression)
                nc_l1r_new=False
                lat = None

            if metadata['zone'] != zone:
                print('UTM zone differs for file {}'.format(image_file))
                continue

            offset = grid_region['off']

            for bi, band in enumerate(bands):
                ds_att = planetscope.get_band_att(metadata, band)
                parname_t = 'rhot_{}'.format(ds_att['wave_name'])
                band_data = planetscope.get_rtoa(image_file, bi+1, band, metadata, sub=sub)

                ## write to L1R NetCDF
                nc_write(nc_file_l1r, parname_t, band_data, dataset_attributes=ds_att, replace_nan=True,
                                  new=nc_l1r_new, attributes=metadata, global_dims=global_dims, offset=offset, nc_compression=nc_compression)
                nc_l1r_new=False

                ## SR file
                if (sr_image_file is not None) & (include_planet_sr):
                        band_data = planetscope.get_rsur(sr_image_file, bi+1,  sub=sub)
                        parname_sr = 'rhosr_{}'.format(ds_att['wave_name'])
                        ## write to L1R NetCDF
                        nc_write(nc_file_l1r, parname_sr, band_data, dataset_attributes=ds_att, replace_nan=True,
                                  new=nc_l1r_new, attributes=metadata, global_dims=global_dims, offset=offset, nc_compression=nc_compression)
                        nc_l1r_new=False

        ## remove the extracted bundle
        if zipped:
             shutil.rmtree(bundle)
             bundle = '{}'.format(bundle_orig)

    return(nc_file_l1r)
