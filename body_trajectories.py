#!/usr/bin/env python
import ephem
from datetime import datetime, date, time, timedelta
import numpy as np
sun = ephem.Sun()
moon = ephem.Moon()

def radec2thetaphi(ra, dec):
    """
    Correspondance between RA/Dec and theta/phi coordinate systems.

    Parameters
    ----------
    ra : float or 1d array
        RA angle in radian.
    dec : float or 1d array
        Dec angle in radian.

    Returns
    ----------
    theta : float or 1d array
        Theta angle in radian
    phi : float or 1d array
        Phi angle in radian
    """
    theta = np.pi / 2 - dec
    phi = float(ra)
    return (theta, phi)

def thetaphi_of_body_oneday(body, date, lon, lat, tz_offset):
    """
    Compute the coordinate (theta, phi) of the Sun at a particular date

    Parameters
    ----------
    date : ephem.Date
        Floating point value used by ephem to represent a date.
        The value is the number of days since 1899 December 31 12:00 UT. When
        creating an instance you can pass in a Python datetime instance,
        timetuple, year-month-day triple, or a plain float.
        Run str() on this object to see the UTC date it represents.
        ...
        WTF?
    lon : string, optional
        Longitute (angle) of the telescope. String form: 0:00:00.0.
    lat : string, optional
        Latitude (angle) of the telescope. String form: 0:00:00.0.

    Returns
    ----------
    theta : float
        Theta angle in radian
    phi : float
        Phi angle in radian
    """
    ob = ephem.Observer()
    ob.lat, ob.lon = lat, lon
    ts = datetime.combine(date, time(12)) - timedelta(hours=tz_offset)
    ob.date = ts
    ALTAZ = alt_az(body, ts, ob)
    RADEC = radec_of(ob, ALTAZ)
    THETAPHI = radec2thetaphi(RADEC[0], RADEC[1])
    return THETAPHI[0], THETAPHI[1]

def alt_az(body, ts, ob):
    """
    Compute altitude (elevation) and azimuth of the Sun.
    """
    ob.date = ts
    body.compute(ob)
    return (body.alt, body.az)

def radec_of(ob, altaz):
    """
    Get RA/Dec from alt/az.
    """
    return ob.radec_of(altaz[1], altaz[0])

def thetaphi_of_body_oneyear(body, lon_observer, lat_observer, year):
    """
    Plot the course of the sun over one year.
    We assume each month is made of 28 days.
    """
    if body == 'sun':
        body = sun
    elif body == 'moon':
        body = moon
    coords = []
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    for month in range(1, 13):
        for day in range(1, 28):
            val = date(year, 1, 1).replace(month=month, day=day)
            t, p = thetaphi_of_body_oneday(
                body,
                val,
                lon_observer.znorm,
                lat_observer.znorm,
                tz_offset=0)
            coords.append((t, p))
    return coords, months
