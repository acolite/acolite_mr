## def parse_metadata
## parses XML metadata from PlanetScope bundle images
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-12
## modifications: 2018-03-15 (QV) moved band stuff to separate keys
##                2018-03-27 (QV) added 0d = 0c
##                2019-06-12 (QV) added 11 = 0f

def parse_metadata(metafile):
    import os
    from acolite.shared import rsr_read, f0_band, rsr_convolute, distance_se
    from acolite import config

    from xml.dom import minidom
    import dateutil.parser

    xmldoc = minidom.parse(metafile)


    metadata = {}

    ## get platform info
    main_tag = 'eop:Platform'
    tags = ["eop:shortName", "eop:serialIdentifier","eop:orbitType"]   
    tags_out = ["platform", 'platform_id', 'orbit']
    for t in xmldoc.getElementsByTagName('eop:Platform'):
        for i,tag in enumerate(tags):
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    metadata[tags_out[i]] = node[0].firstChild.nodeValue

    if "RapidEye" in metadata['platform_id']:
        metadata['SATELLITE_ID'] = metadata['platform']
        metadata['SATELLITE_SENSOR'] = metadata['platform_id']
        metadata['LUT_SENSOR'] = "RapidEye"
        metadata['SENSOR']='RapidEye'
        metadata['SATELLITE_PREFIX'] = 're'
        bnames = {1:'Blue',2:'Green',3:'Red',4:'RedEdge', 5:'NIR'}
        metadata['BANDS_REDNIR'] = ['Red','RedEdge','NIR']
        metadata['BANDS_VIS'] = ['Blue','Green','Red']
        metadata['BANDS_NIR'] = ['RedEdge','NIR']
        metadata['BANDS_BESTFIT'] = ['Red','RedEdge','NIR']
        metadata['BANDS_ALL'] = ['Blue','Green','Red','RedEdge','NIR']


    if 'PlanetScope' in metadata['platform']:
        metadata['SATELLITE_ID'] = metadata['platform_id'][0:2]
        if metadata['SATELLITE_ID'] == '10': metadata['SATELLITE_ID'] = "0f"
        if metadata['SATELLITE_ID'] == '11': metadata['SATELLITE_ID'] = "0f"
        if metadata['SATELLITE_ID'] == '0d': metadata['SATELLITE_ID'] = "0c"

        metadata['SATELLITE_SENSOR'] = '{}_{}'.format(metadata['platform'], metadata['platform_id'])
        metadata['LUT_SENSOR'] = '{}_{}'.format(metadata['platform'], metadata['SATELLITE_ID'])
        metadata['SENSOR']='PlanetScope'
        metadata['SATELLITE_PREFIX'] = 'ps'
        bnames = {1:'Blue',2:'Green',3:'Red',4:'NIR'}
        metadata['BANDS_REDNIR'] = ['Red','NIR']
        metadata['BANDS_VIS'] = ['Blue','Green','Red']
        metadata['BANDS_NIR'] = ['NIR']
        metadata['BANDS_BESTFIT'] = ['Red','NIR']
        metadata['BANDS_ALL'] = ['Blue','Green','Red','NIR']

    ## get acquisition info
    main_tag = 'eop:acquisitionParameters'

    tags = ["eop:orbitDirection", "eop:incidenceAngle",
            "opt:illuminationAzimuthAngle", "opt:illuminationElevationAngle", 
            "{}:azimuthAngle".format(metadata['SATELLITE_PREFIX']),
            "{}:spaceCraftViewAngle".format(metadata['SATELLITE_PREFIX']),
            "{}:acquisitionDateTime".format(metadata['SATELLITE_PREFIX'])]
    tags_out = ["orbit",'ViewingIncidence', 
                'SunAzimuth', 'SunElevation', 
                'ViewingAzimuth', 'ViewZenith', 'isotime']

    for t in xmldoc.getElementsByTagName(main_tag):
        for i,tag in enumerate(tags):
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    if tag in ["eop:orbitDirection","{}:acquisitionDateTime".format(metadata['SATELLITE_PREFIX'])]: 
                        val=node[0].firstChild.nodeValue
                    else: 
                        val=float(node[0].firstChild.nodeValue)
                    metadata[tags_out[i]] = val

    metadata['THS'] = 90.-metadata['SunElevation']
    metadata['PHIS'] = metadata['SunAzimuth']
    metadata['THV'] = abs(metadata['ViewZenith'])
    metadata['PHIV'] = metadata['ViewingAzimuth']
    metadata['AZI'] = abs(metadata['PHIS'] - metadata['PHIV'])
    while (metadata['AZI']>180):
        metadata['AZI']=abs(180-metadata['AZI'])
    
    metadata["TIME"]=dateutil.parser.parse(metadata["isotime"])
    metadata["DOY"] = metadata["TIME"].strftime('%j')
    metadata["SE_DISTANCE"] = distance_se(metadata['DOY'])
    metadata["isodate"]=metadata["TIME"].strftime('%Y-%m-%dT%H:%M:%SZ')

    ## get band data
    main_tag='{}:bandSpecificMetadata'.format(metadata['SATELLITE_PREFIX'])
    bands = {}
    tags = ["{}:bandNumber".format(metadata['SATELLITE_PREFIX']),
            '{}:radiometricScaleFactor'.format(metadata['SATELLITE_PREFIX']),
            '{}:reflectanceCoefficient'.format(metadata['SATELLITE_PREFIX'])]
    tags_out = ["band_idx",'to_radiance', 'to_reflectance']

    ## import RSR to get F0
    #pp_path = os.path.dirname(pp.__file__)
    #if metadata['LUT_SENSOR'] == 'PlanetScope_0d':
    #    print('Sensor {} not yet supported'.format(metadata['LUT_SENSOR']))
    #    return(metadata)

    rsr_file="{}/RSR/{}.txt".format(config['pp_data_dir'],metadata['LUT_SENSOR'])
    rsr, rsr_bands = rsr_read(file=rsr_file)

    for t in xmldoc.getElementsByTagName(main_tag):
        band = {}
        for i,tag in enumerate(tags):
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    val = float(node[0].firstChild.nodeValue)
                    band[tags_out[i]] = val

        band['band_name']= bnames[band['band_idx']]
        
        band_rsr=rsr[band['band_name']]['response']
        band_wave=[float(i*1000) for i in rsr[band['band_name']]['wave']]
        
        f0=f0_band(band_wave, band_rsr)
        wave = rsr_convolute(rsr[band['band_name']]['wave'], 
                             rsr[band['band_name']]['wave'], 
                             band_rsr, rsr[band['band_name']]['wave'])
        
        band['f0'] = f0
        band['wave'] = wave
        band['wave_name'] = str(round(int(band['wave']*1000.),2))

        bands[band['band_name']] = band
        
        
    #metadata['bands']=bands
    for band in bands:
        for key in bands[band]:
            bk='{}-{}'.format(band, key)
            metadata[bk] = bands[band][key]

    ## get product info
    main_tag = '{}:spatialReferenceSystem'.format(metadata['SATELLITE_PREFIX'])
    tags = ["{}:epsgCode".format(metadata['SATELLITE_PREFIX']),
            "{}:geodeticDatum".format(metadata['SATELLITE_PREFIX']),
            "{}:projection".format(metadata['SATELLITE_PREFIX']),
            "{}:projectionZone".format(metadata['SATELLITE_PREFIX'])]
    tags_out = ["epsg",'datum', 'projection', 'zone']
    for t in xmldoc.getElementsByTagName(main_tag):
        for i,tag in enumerate(tags):
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    val=node[0].firstChild.nodeValue
                    metadata[tags_out[i]] = val
                    
    ## get resolution info
    #main_tag = 'ps:Sensor'
    #tags = ["eop:resolution"]
    #tags_out = ["resolution"]
    main_tag = '{}:ProductInformation'.format(metadata['SATELLITE_PREFIX'])
    tags = ["{}:numRows".format(metadata['SATELLITE_PREFIX']), 
            "{}:numColumns".format(metadata['SATELLITE_PREFIX']), 
            "{}:numBands".format(metadata['SATELLITE_PREFIX']), 
            "{}:rowGsd".format(metadata['SATELLITE_PREFIX']), 
            "{}:columnGsd".format(metadata['SATELLITE_PREFIX'])]
    tags_out = ["nrow","ncol","nband","resolution_row", "resolution_col"]
    for t in xmldoc.getElementsByTagName(main_tag):
        for i,tag in enumerate(tags):
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    val=node[0].firstChild.nodeValue
                    metadata[tags_out[i]] = val
    metadata['resolution'] = (float(metadata['resolution_row']),float(metadata['resolution_col']))
    metadata['dims'] = (int(metadata['ncol']),int(metadata['nrow']))

    ## get bounding box
    main_tag = '{}:geographicLocation'.format(metadata['SATELLITE_PREFIX'])
    tags = ["{}:topLeft".format(metadata['SATELLITE_PREFIX']),
            "{}:topRight".format(metadata['SATELLITE_PREFIX']),
            "{}:bottomRight".format(metadata['SATELLITE_PREFIX']),
            "{}:bottomLeft".format(metadata['SATELLITE_PREFIX'])]
    tags_out = ["UL",'UR', 'LR', 'LL']
    for t in xmldoc.getElementsByTagName(main_tag):
        for i,tag in enumerate(tags):
                node = t.getElementsByTagName(tag)
                for j,tag2 in enumerate(['{}:latitude'.format(metadata['SATELLITE_PREFIX']),'{}:longitude'.format(metadata['SATELLITE_PREFIX'])]):
                    node2 = node[0].getElementsByTagName(tag2)
                    if len(node2) > 0:
                        val=node2[0].firstChild.nodeValue
                        tout = '{}_{}'.format(tags_out[i], tag2.split(':')[1].upper())
                        metadata[tout] = float(val)


    return(metadata)

