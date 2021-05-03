"""
Internal validation functions.

:copyright: (c) 2015 by Mark Richards.
:license: BSD 3-Clause, see LICENSE.txt for more details.
"""

# Internal constants
# Latitude
_MINLAT = -90.0
_MAXLAT = 90.0

# Longitude
_MINLON = -180.0
_MAXLON = 180.0

# Solar declination
_MINSOLDEC = -23.5
_MAXSOLDEC = 23.5

# Sunset hour angle
_MINSHA = 0.0
_MAXSHA = 180


def check_day_hours(hours, arg_name):
    """
    Check that *hours* is in the range 1 to 24.
    """
    if not 0 <= hours <= 24:
        raise ValueError(
            '{0} should be in range 0-24: {1!r}'.format(arg_name, hours))


def check_doy(doy):
    """
    Check day of the year is valid.
    """
    if not 1 <= doy <= 366:
        raise ValueError(
            'Day of the year (doy) must be in range 1-366: {0!r}'.format(doy))


def check_latitude_rad(latitude):
    if not _MINLAT <= latitude <= _MAXLAT:
        raise ValueError('latitude outside valid range {0!r} to {1!r} rad: {2!r}'.format(_MINLAT, _MAXLAT, latitude))


def check_longitude(longitude):
    if not _MINLON <= longitude <= _MAXLON:
        raise ValueError('longitude outside valid range {0!r} to {1!r} rad: {2!r}'.format(_MINLAT, _MAXLAT, longitude))


def check_sol_dec_rad(sd):
    """
    Solar declination can vary between -23.5 and +23.5 degrees.

    See http://mypages.iit.edu/~maslanka/SolarGeo.pdf
    """
    if not _MINSOLDEC <= sd <= _MAXSOLDEC:
        raise ValueError(
            'solar declination outside valid range {0!r} to {1!r} rad: {2!r}'.format(_MINSOLDEC, _MAXSOLDEC, sd))


def check_sunset_hour_angle_rad(sha):
    """
    Sunset hour angle has the range 0 to 180 degrees.

    See http://mypages.iit.edu/~maslanka/SolarGeo.pdf
    """
    if not _MINSHA <= sha <= _MAXSHA:
        raise ValueError(
            'sunset hour angle outside valid range {0!r} to {1!r} rad: {2!r}'.format(_MINSHA, _MAXSHA, sha))
