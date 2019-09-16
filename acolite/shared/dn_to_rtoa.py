## def dn_to_rtoa
## convert DN to TOA reflectance
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-06-29
## modifications:

def dn_to_rtoa(data, f0, sun_zenith, slope=1.0, offset=0.0, d=1.0):
    from numpy import pi
    from numpy import cos
    from numpy import ndarray
    
    dtor = pi/180.
    data = data.astype('float')
    if slope != 1.0: data *= slope
    if offset != 0.0: data += offset
    data *= (pi * d * d) / (f0 * cos(sun_zenith*dtor))
    return data
