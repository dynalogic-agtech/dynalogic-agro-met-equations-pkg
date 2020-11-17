"""
Library of functions for meteorological equations that can be used for general calculations.

:author: Bruno Ducraux
:copyright: (c) 2019 by Dynalogic.
:license: BSD 3-Clause, see LICENSE.txt for more details.
"""

import math
from datetime import date


def frost_point(temperature_mean, _dew_point):
    """Compute the frost point in degrees Celsius
    :param temperature_mean: current ambient temperature in degrees Celsius
    :type temperature_mean: float
    :param _dew_point: current dew point in degrees Celsius
    :type _dew_point: float
    :return: the frost point in degrees Celsius
    :rtype: float
    """
    dew_point_k = 273.15 + _dew_point
    t_air_k = 273.15 + temperature_mean
    frost_point_k = dew_point_k - t_air_k + 2671.02 / ((2954.61 / t_air_k) + 2.193665 * math.log(t_air_k) - 13.3448)
    return frost_point_k - 273.15


def dew_point(temperature_mean, relative_humidity_mean):
    """Compute the dew point in degrees Celsius
    :param temperature_mean: current ambient temperature in degrees Celsius
    :type temperature_mean: float
    :param relative_humidity_mean: relative humidity in %
    :type relative_humidity_mean: float
    :return: the dew point in degrees Celsius
    :rtype: float
    """
    ref_a = 17.27
    ref_b = 237.3
    alpha = ((ref_a * temperature_mean) / (ref_b + temperature_mean)) + math.log(relative_humidity_mean / 100.0)
    return (ref_b * alpha) / (ref_a - alpha)


# todo create UnitTest for all degree day cases
def degree_day(temperature_mean, temperature_basal_inferior, **kwargs):
    # optional args for more accurate value
    temperature_basal_superior = kwargs.get('temperature_basal_superior')
    temperature_max = kwargs.get('temperature_max')
    temperature_min = kwargs.get('temperature_min')

    if all(key in kwargs for key in ('temperature_basal_superior', 'temperature_max', 'temperature_min')):
        # complete formula
        deg_day = None

        if temperature_basal_superior > temperature_max > temperature_min > temperature_basal_inferior:
            deg_day = ((temperature_max - temperature_min) / 2) + temperature_min + temperature_basal_inferior
        elif temperature_basal_superior > temperature_max > temperature_basal_inferior > temperature_min:
            tmp = math.pow((temperature_max - temperature_basal_inferior), 2)
            tmp1 = (2 * (temperature_max - temperature_min))
            deg_day = tmp / tmp1
        elif temperature_basal_superior > temperature_basal_inferior > temperature_max > temperature_min:
            deg_day = 0
        elif temperature_max > temperature_basal_superior > temperature_min > temperature_basal_inferior:
            tdif = temperature_max - temperature_min
            tmp = (2 * tdif) * (temperature_min - temperature_basal_inferior)
            tmp1 = math.pow(tdif, 2)
            tmp2 = math.pow((temperature_max - temperature_basal_superior), 2)
            tmp3 = 2 * tdif
            deg_day = (tmp + tmp1 - tmp2) / tmp3
        elif temperature_max > temperature_basal_superior > temperature_basal_inferior > temperature_min:
            tmp = math.pow((temperature_max - temperature_basal_inferior), 2)
            tmp1 = math.pow((temperature_max - temperature_basal_superior), 2)
            tmp2 = 2 * (temperature_max - temperature_min)
            deg_day = (tmp - tmp1) / tmp2
        return deg_day
    else:
        # simplified formula
        print("simplified formula")
        return temperature_mean - temperature_basal_inferior


def days_passed_on_current_year():
    """
    Calculate the days passed until today date
    :return: int
    """
    today = date.today()
    return today.timetuple().tm_yday


def sunrise_sunset_hour(_sunset_hour_angle):
    """
    :param _sunset_hour_angle:
    :return: tuple with sunset hour and sunrise hour
    """

    # Sunrise
    hour_to_sunrise = 12 + (- round(_sunset_hour_angle / 15, 2))

    split_hour = math.modf(hour_to_sunrise)

    sunrise_hour = split_hour[1]
    sunrise_minutes = 60 * split_hour[0]
    sunrise = "{}:{}".format(int(sunrise_hour), int(sunrise_minutes))

    # Sunset
    hour_to_sunset = 12 - (- round(_sunset_hour_angle / 15, 2))

    split_hour = math.modf(hour_to_sunset)

    sunset_hour = split_hour[1]
    sunset_minutes = 60 * split_hour[0]
    sunset = "{}:{}".format(int(sunset_hour), int(sunset_minutes))

    return sunrise, sunset


def vapour_pressure_deficit(es, ea):
    """
    VPD - Vapour Pressure Deficit
    :param es:
    :param ea: calculated by vp_from_rhmin_rhmax()
    :return:
    """
    return es - ea
