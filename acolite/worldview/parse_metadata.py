## def parse_metadata
## parses XML metadata from WorldView2 bundle images
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-27
## modifications: QV 2019-05-06
##                2020-03-16 (QV) fixed tile metadata reads
##                2020-07-26 (QV) added WV3
##                2020-11-24 (QV) added fix for WorldView2 "L1B" data

def parse_metadata(metafile):
    from acolite.shared import rsr_read, f0_band, rsr_convolute
    import acolite as ac

    import os, sys, fnmatch, dateutil.parser
    from xml.dom import minidom

    if not os.path.isfile(metafile):
        print('Metadata file not found.')
        sys.exit()

    try:
        xmldoc = minidom.parse(metafile)
    except:
        print('Error opening metadata file.')
        sys.exit()

    metadata = {}
    #metadata['SATELLITE'] = 'WorldView2'
    #metadata['SENSOR'] = 'WorldView2'

    ## get image information
    metadata_tags = ['SATID', 'FIRSTLINETIME', 'NUMROWS','NUMCOLUMNS','PRODUCTLEVEL'
                    "MININTRACKVIEWANGLE", "MAXINTRACKVIEWANGLE", "MEANINTRACKVIEWANGLE",
                    "MINCROSSTRACKVIEWANGLE", "MAXCROSSTRACKVIEWANGLE", "MEANCROSSTRACKVIEWANGLE",
                    "MINOFFNADIRVIEWANGLE", "MAXOFFNADIRVIEWANGLE", "MEANOFFNADIRVIEWANGLE",
                    "MINSUNAZ","MAXSUNAZ","MEANSUNAZ",
                    "MINSUNEL","MAXSUNEL","MEANSUNEL",
                    "MINSATAZ","MAXSATAZ","MEANSATAZ",
                    "MINSATEL","MAXSATEL","MEANSATEL",
                    "EARLIESTACQTIME","LATESTACQTIME"]

    for tag in metadata_tags:
        node = xmldoc.getElementsByTagName(tag)
        if len(node) > 0: metadata[tag] = node[0].firstChild.nodeValue

    if 'SATID' in metadata:
        if metadata['SATID'] == 'WV03':
            metadata['SATELLITE'] = 'WorldView3'
            metadata['SENSOR'] = 'WorldView3'
            metadata["TIME"]=dateutil.parser.parse(metadata["FIRSTLINETIME"])
            band_names=['COASTAL','BLUE','GREEN','YELLOW','RED','REDEDGE','NIR1','NIR2',
                        'SWIR1','SWIR2','SWIR3','SWIR4','SWIR5','SWIR6','SWIR7','SWIR8']
            band_indices=[1,2,3,4,5,6,7,8,
                          9, 10, 11, 12, 13, 14, 15, 16]
            band_indices=[1,2,3,4,5,6,7,8,
                          1,2,3,4,5,6,7,8]

            band_tag_names = ["BAND_C","BAND_B","BAND_G","BAND_Y","BAND_R","BAND_RE","BAND_N", "BAND_N2",
                              "BAND_S1", "BAND_S2", "BAND_S3", "BAND_S4", "BAND_S5", "BAND_S6", "BAND_S7", "BAND_S8"]

        if metadata['SATID'] == 'WV02':
            metadata['SATELLITE'] = 'WorldView2'
            metadata['SENSOR'] = 'WorldView2'
            try:
                ## "L2A" data
                metadata["TIME"]=dateutil.parser.parse(metadata["EARLIESTACQTIME"])
            except:
                ## "L1B" data
                metadata["TIME"]=dateutil.parser.parse(metadata["FIRSTLINETIME"])
            band_names=['COASTAL','BLUE','GREEN','YELLOW','RED','REDEDGE','NIR1','NIR2']
            band_indices=[1,2,3,4,5,6,7,8]
            band_tag_names = ["BAND_C","BAND_B","BAND_G","BAND_Y","BAND_R","BAND_RE","BAND_N", "BAND_N2"]

    ## common stuff
    metadata["DOY"] = metadata["TIME"].strftime('%j')
    metadata["THS"] = 90.-float(metadata['MEANSUNEL'])
    metadata["THV"] = 90.-float(metadata['MEANSATEL'])
    metadata["AZI"] = abs(float(metadata['MEANSATAZ'])-float(metadata['MEANSUNAZ']))
    if metadata["AZI"] > 180.: metadata["AZI"]-=180.

    band_tags = ["ULLON","ULLAT","ULHAE",
                "URLON","URLAT","URHAE",
                "LRLON","LRLAT","LRHAE",
                "LLLON","LLLAT","LLHAE",
                "ABSCALFACTOR","EFFECTIVEBANDWIDTH","TDILEVEL"]

    ## import RSR to get F0
    sensor=metadata['SENSOR']
    pp_path = ac.config['pp_data_dir']
    rsr_file = pp_path+'/RSR/'+sensor+'.txt'
    rsr, rsr_bands = rsr_read(file=rsr_file)

    ## read band information of spatial extent
    band_values={}
    for b,band_tag in enumerate(band_tag_names):
        band_rsr=rsr[band_names[b]]['response']
        band_wave=[i*1000 for i in rsr[band_names[b]]['wave']]
        f0=f0_band(band_wave, band_rsr)
        wave = rsr_convolute(rsr[band_names[b]]['wave'], rsr[band_names[b]]['wave'], band_rsr, band_wave)

        band_data = {'F0':f0, 'wave':wave, 'name':band_names[b], 'index':band_indices[b]}

        band_data['wave_name'] = str(round(int(band_data['wave']*1000.),2))
        ## there are two tags in WV3 metadata
        for t in xmldoc.getElementsByTagName(band_tag):
            for tag in band_tags:
                    node = t.getElementsByTagName(tag)
                    if len(node) > 0:
                        band_data[tag]=float(node[0].firstChild.nodeValue)
            #band_data['band_index']=t
        if len(band_data)>5:
            band_values[band_tag]=band_data

    metadata['BAND_INFO'] = band_values


    ## get tile information
    tile_tags = ["FILENAME",
                "ULCOLOFFSET","ULROWOFFSET","URCOLOFFSET","URROWOFFSET",
                "LRCOLOFFSET","LRROWOFFSET","LLCOLOFFSET","LLROWOFFSET",
                "ULLON","ULLAT","URLON","URLAT",
                "LRLON","LRLAT","LLLON","LLLAT",
                "ULX","ULY","URX","URY",
                "LRX","LRY","LLX","LLY"]
    tile_values=[]
    for t in xmldoc.getElementsByTagName('TILE'):
        tile = {}
        for tag in tile_tags:
                node = t.getElementsByTagName(tag)
                if len(node) > 0:
                    if tag == "FILENAME": val=node[0].firstChild.nodeValue
                    else: val=float(node[0].firstChild.nodeValue)
                    tile[tag]=val
        tile_values.append(tile)
    metadata['TILE_INFO'] = tile_values

    return metadata
