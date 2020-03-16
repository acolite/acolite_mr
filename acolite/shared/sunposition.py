## def sunposition
## compute sun position
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-11-27
## modifications:

def sunposition(isotime, lon, lat):
    from numpy import floor, pi, sin, cos, arctan, arcsin, isfinite
    dtor = pi/180.
    import dateutil.parser
    time = dateutil.parser.parse(isotime)
    doy = float(time.strftime('%j'))

    # Get Julian date - 2400000
    ftime = time.hour + time.minute/60. + time.second / 3600.
    delta = time.year - 1949
    leap = floor(delta / 4.) ;# former leapyears
    jd = 32916.5 + delta * 365. + leap + doy + ftime / 24.

    # The input to the Atronomer's almanach is the difference between
    # the Julian date and JD 2451545.0 (noon, 1 January 2000)
    almanach_time = jd - 51545.

    # Ecliptic coordinates
    # Mean longitude
    mnlong = 280.460 + .9856474 * almanach_time
    mnlong = mnlong % 360.
    if mnlong<0:mnlong+=360.

    # Mean anomaly
    mnanom = 357.528 + .9856003 * almanach_time
    mnanom = mnanom % 360.
    if mnanom<0:mnanom+=360.
    mnanom = mnanom * dtor

    # Ecliptic longitude and obliquity of ecliptic
    eclong = mnlong + 1.915 * sin(mnanom) + 0.020 * sin(2 * mnanom)
    eclong = eclong % 360.
    if eclong<0:eclong+=360.
    oblqec = 23.429 - 0.0000004 * almanach_time
    oblqec = 23.439 - 0.0000004 * almanach_time
    eclong = eclong * dtor
    oblqec = oblqec * dtor

    # Celestial coordinates
    # Right ascension and declination
    num = cos(oblqec) * sin(eclong)
    den = cos(eclong)
    ra = arctan(num / den)
    if den < 0: ra+=pi
    if (den > 0) & (num > 0): ra+= pi * 2.
    dec = arcsin(sin(oblqec) * sin(eclong))

    # Local coordinates
    # Greenwich mean sidereal time
    gmst = 6.697375 + .0657098242 * almanach_time + ftime
    gmst = gmst % 24.
    if gmst < 0: gmst+=24.

    # Local mean sidereal time
    lmst = gmst + lon / 15.
    lmst = lmst % 24.
    if lmst < 0: lmst += 24.
    lmst = lmst * 15. * dtor

    # Hour angle
    ha = lmst - ra
    if ha < -pi: ha += pi * 2.
    if ha > pi: ha -= pi * 2.

    # Latitude to radians
    lat = lat * dtor

    # Azimuth and elevation
    el = arcsin(sin(dec) * sin(lat) + cos(dec) * cos(lat) * cos(ha))
    az = arcsin(-cos(dec) * sin(ha) / cos(el))

    cosAzPos = 0 < sin(dec) - sin(el) * sin(lat)
    sinAzNeg = sin(az) < 0
    if (cosAzPos) and (sinAzNeg): az += pi * 2
    if (cosAzPos == 0): az = pi - az
    if isfinite(az) is False:az=pi/2

    distance = 1.00014 - 0.01671 * cos(mnanom) - 0.00014 * cos(2.*mnanom)

    return({'elevation':el/dtor, 'azimuth':az/dtor, 'zenith':90.-el/dtor,
            'distance':distance, 'mnanom':mnanom, 'doy':doy})
