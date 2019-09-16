## def get_rtoa
## reads WorldView2 DN and converts to RTOA - needs metadata
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-01-27
## modifications:

def get_rtoa(file, idx, name, metadata, sun_zenith=0.0, se_distance=1.0, sub=None):
    from acolite.worldview import read_band
    from acolite.shared import dn_to_rtoa
    from numpy import nan, uint16
    
    data = read_band(file, idx, sub=sub)
    nodata = data == uint16(0)

    data = (data.astype(float) * float(metadata['BAND_INFO'][name]['ABSCALFACTOR']))/float(metadata['BAND_INFO'][name]['EFFECTIVEBANDWIDTH'])
    
    data = dn_to_rtoa(data, metadata['BAND_INFO'][name]['F0']*10., sun_zenith, 
                slope=1.0, offset=0.0, d=se_distance)
    data[nodata] = nan
    return data
