## def get_band_att
## gets band specific attributes from metadata 
## written by Quinten Vanhellemont, RBINS
## 2018-03-15
## modifications:

def get_band_att(metadata, band):
    b_att = {}
    for mk in metadata:
        if '{}-'.format(band) not in mk: continue
        bk, bp = mk.split('-')
        b_att[bp] = metadata[mk]
    b_att['wavelength']=float(metadata['{}-{}'.format(band, 'wave_name')])
    return(b_att)

