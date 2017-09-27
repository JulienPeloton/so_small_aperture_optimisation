#!/usr/bin/python
"""
Script to predict the trajectory of a body (Sun, Moon, ...) on the sky
[(theta, phi) coordinates] over the course of a year.

Author: Julien Peloton, j.peloton@sussex.ac.uk
"""
from __future__ import division, absolute_import, print_function

import ephem
from datetime import datetime, date, time, timedelta
import numpy as np

## Define available objects
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

def radec_of(ob, altaz):
    """
    Get RA/Dec from alt/az.

    Parameters
    ----------
    ob : ephem.Observer instance
        ephem class with properties of the observer.
    altaz : tuple
        Tuple containing (alt, az) in radians.

    Returns
    ----------
    RA : float
        RA of the body [radian].
    Dec : float
        Dec of the body [radian].
    """
    return ob.radec_of(altaz[1], altaz[0])

def alt_az(body, ts, ob):
    """
    Compute altitude (elevation) and azimuth of the body.

    Parameters
    ----------
    body : ephem.<body> instance
        <body> can be Sun(), Moon(), and so on. See ephem.
    ts : Python datetime instance
        Floating point value used by ephem to represent a date.
        The value is the number of days since 1899 December 31 12:00 UT. When
        creating an instance you can pass in a Python datetime instance,
        timetuple, year-month-day triple, or a plain float.
        Run str() on this object to see the UTC date it represents.
        ...
        WTF?
    ob : ephem.Observer instance
        ephem class with properties of the observer.

    Returns
    ----------
    alt : float
        Altitude (elevation) of the body at the time of observation as
        seen by the observer [radian].
    az : float
        Azimuth of the body at the time of observation as seen
        by the observer [radian].
    """
    ob.date = ts
    body.compute(ob)
    return (body.alt, body.az)

def thetaphi_of_body_oneday(body, date, lon, lat, tz_offset):
    """
    Compute the coordinate (theta, phi) of the Sun at a particular date

    Parameters
    ----------
    body : ephem.<body> instance
        <body> can be Sun(), Moon(), and so on. See ephem.
    date : Python datetime instance
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
        Theta angle [radian].
    phi : float
        Phi angle [radian].
    """
    ob = ephem.Observer()
    ob.lat, ob.lon = lat, lon
    ts = datetime.combine(date, time(12)) - timedelta(hours=tz_offset)
    ob.date = ts
    ALTAZ = alt_az(body, ts, ob)
    RADEC = radec_of(ob, ALTAZ)
    THETAPHI = radec2thetaphi(RADEC[0], RADEC[1])
    return THETAPHI[0], THETAPHI[1]

def thetaphi_of_body_oneyear(body, lon_observer, lat_observer, year):
    """
    Plot the course of the sun over one year.
    We assume each month is made of 28 days.

    Parameters
    ----------
    body : ephem.<body> instance
        <body> can be Sun(), Moon(), and so on. See ephem.
    lon_observer : str
        Longitute (angle) of the telescope. String form: 0:00:00.0.
    lat_observer : str
        Latitude (angle) of the telescope. String form: 0:00:00.0.
    year : int
        Year of observation

    Returns
    ----------
    coords : 2D array of floats
        The (theta, phi) coordinate for all days of the year. We compute
        the coordinate for only 28 days per month.
    months : list of strings
        Name of the 12 months.

    Examples
    ----------
    >>> coords, months = thetaphi_of_body_oneyear('sun', '-67:46.816',
    ...     '-22:56.396', 2013)
    >>> print(round(coords[0][0], 2), round(coords[0][1], 2))
    1.97 4.92

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
                lon_observer,
                lat_observer,
                tz_offset=0)
            coords.append((t, p))
    return coords, months


if __name__ == "__main__":
    import doctest
    doctest.testmod()
