## def get_rsur
## reads PlanetScope SR DN and converts to Surface Reflectance
## written by Quinten Vanhellemont, RBINS
## 2018-03-12
## modifications:

def get_rsur(file, band, sub=None):
    from numpy import nan, uint16
    
    data = read_band(file, band, sub=sub)
    nodata = data == uint16(0)

    data = data.astype(float) / float(10000.)
    
    data[nodata] = nan
    return(data)
