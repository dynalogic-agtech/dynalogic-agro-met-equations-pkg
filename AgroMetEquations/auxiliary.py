"""
Library of functions for meteorological equations that can be used for general calculations and conversions.

:author: Bruno Ducraux
:copyright: (c) 2019 by Dynalogic.
:license: BSD 3-Clause, see LICENSE.txt for more details.
"""

import math
from datetime import time, date
import calendar
from AgroMetEquations.config import (_MIN_LATITUDE,
                                     _MAX_LATITUDE,
                                     _MIN_LONGITUDE,
                                     _MAX_LONGITUDE,
                                     _MIN_SOLAR_DECLINATION,
                                     _MAX_SOLAR_DECLINATION,
                                     _MIN_SUNSET_HOUR_ANGLE,
                                     _MAX_SUNSET_HOUR_ANGLE,
                                     SOLAR_CONSTANT,
                                     REFERENTIAL_LONGITUDE,
                                     STEFAN_BOLTZMANN_CONSTANT,
                                     _MONTHDAYS,
                                     _LEAP_MONTHDAYS)


# ######################################################################################################################
# VALIDATIONS
# ######################################################################################################################

def check_day_hours(hours, arg_name):
    """
    Check that *hours* is in the range 1 to 24.
    """
    if not 0 <= hours <= 24:
        raise ValueError(
            '{0} should be in range 0-24: {1!r}'.format(arg_name, hours))


def check_day_of_year(day_of_year):
    """
    Check day of the year is valid.
    """
    if not 1 <= day_of_year <= 366:
        raise ValueError(
            'Day of the year (doy) must be in range 1-366: {0!r}'.format(day_of_year))


def check_latitude(latitude):
    if not _MIN_LATITUDE <= latitude <= _MAX_LATITUDE:
        raise ValueError('latitude outside valid range {0!r} to {1!r} rad: {2!r}'.format(_MIN_LATITUDE,
                                                                                         _MAX_LATITUDE,
                                                                                         latitude))


def check_longitude(longitude):
    if not _MIN_LONGITUDE <= longitude <= _MAX_LONGITUDE:
        raise ValueError('longitude outside valid range {0!r} to {1!r} rad: {2!r}'.format(_MIN_LATITUDE,
                                                                                          _MAX_LATITUDE,
                                                                                          longitude))


def check_solar_declination(sd):
    """
    Solar declination can vary between -23.5 and +23.5 degrees.
    See https://mypages.iit.edu/~maslanka/SolarGeo.pdf
    """
    if not _MIN_SOLAR_DECLINATION <= sd <= _MAX_SOLAR_DECLINATION:
        raise ValueError(
            'solar declination outside valid range {0!r} to {1!r} rad: {2!r}'.format(_MIN_SOLAR_DECLINATION,
                                                                                     _MAX_SOLAR_DECLINATION,
                                                                                     sd))


def check_sunset_hour_angle(sha):
    """
    Sunset hour angle has the range 0 to 180 degrees.
    See https://mypages.iit.edu/~maslanka/SolarGeo.pdf
    """
    if not _MIN_SUNSET_HOUR_ANGLE <= sha <= _MAX_SUNSET_HOUR_ANGLE:
        raise ValueError(
            'sunset hour angle outside valid range {0!r} to {1!r} rad: {2!r}'.format(_MIN_SUNSET_HOUR_ANGLE,
                                                                                     _MAX_SUNSET_HOUR_ANGLE,
                                                                                     sha))


# ######################################################################################################################
# CONVERSIONS
# ######################################################################################################################


def celsius2kelvin(celsius: float) -> float:
    """
    Convert temperature in degrees Celsius to degrees Kelvin.

    :param celsius: float -> Degrees Celsius
    :return: Degrees Kelvin
    :rtype: float
    """
    return celsius + 273.15


def kelvin2celsius(kelvin):
    """
    Convert temperature in degrees Kelvin to degrees Celsius.

    :param kelvin: float -> Degrees Kelvin
    :return: Degrees Celsius
    :rtype: float
    """
    return kelvin - 273.15


def deg2rad(degrees: float) -> float:
    """
    Convert angular degrees to radians

    :param degrees: float -> Value in degrees to be converted.
    :return: Value in radians
    :rtype: float
    """
    return degrees * (math.pi / 180.0)


def rad2deg(radians):
    """
    Convert radians to angular degrees

    :param radians: Value in radians to be converted.
    :return: Value in angular degrees
    :rtype: float
    """
    return radians * (180.0 / math.pi)


def watt2mj15min(watt):
    """
    Convert watt/m2 to MJ/m2/15min
    :param watt:
    :return: MJ/m2/15min
    :rtype: float
    """
    return watt * 0.0009


def daily_to_15min(daily_value):
    """
    Transform daily value to 15min value
    (daily_value / 24) / 4 = daily_value / 96
    :param daily_value:
    :return: 15min value
    """
    return daily_value / 96


def dms2dd(degrees, minutes, seconds, direction: str):
    """
    Convert degrees, minutes, seconds to decimal degrees.
    :param degrees: float
    :param minutes: float
    :param seconds: float
    :param direction: str -> 'N', 'S', 'E', 'W'
    :return:
    """
    dd = float(degrees) + float(minutes) / 60 + float(seconds) / (60 * 60)
    if direction == 'S' or direction == 'W':
        dd *= -1
    return dd


def dd2dms(dd) -> tuple:
    """
    Convert decimal degrees to degrees, minutes, seconds.
    :param dd:
    :return: tuple -> degrees, minutes, seconds
    """
    negative = dd < 0
    dd = abs(dd)
    minutes, seconds = divmod(dd * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    if negative:
        if degrees > 0:
            degrees = -degrees
        elif minutes > 0:
            minutes = -minutes
        else:
            seconds = -seconds

    return degrees, minutes, seconds


def energy2evap(energy, _latent_heat):
    """
    Convert energy (e.g. radiation energy) in MJ m-2 time to the equivalent
    evaporation, assuming a grass reference crop.

    Energy is converted to equivalent evaporation using a conversion
    factor equal to the inverse of the latent heat of vapourisation
    (1 / _latent_heat).

    :param energy: Energy e.g. radiation or heat flux [MJ m-2 time].
    :param _latent_heat Calculated by latent_heat(temperature_mean)
    :return: Equivalent evaporation [mm time-1].
    :rtype: float
    """
    return (1 / _latent_heat) * energy


# ######################################################################################################################
# GENERAL CALCULATIONS
# ######################################################################################################################


def get_frost_point(temperature_mean, _dew_point):
    """Compute the frost point in degrees Celsius
    :param temperature_mean: current ambient temperature in degrees Celsius
    :type temperature_mean: float
    :param _dew_point: current dew point in degrees Celsius
    :type _dew_point: float
    :return: the frost point in degrees Celsius
    :rtype: float
    """
    dew_point_k = 273.15 + _dew_point
    dew_point_k = 273.15 + _dew_point
    t_air_k = 273.15 + temperature_mean
    frost_point_k = dew_point_k - t_air_k + 2671.02 / ((2954.61 / t_air_k) + 2.193665 * math.log(t_air_k) - 13.3448)
    return frost_point_k - 273.15


def get_dew_point(temperature_mean, relative_humidity_mean):
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


def get_degree_day(temperature_mean, temperature_basal_inferior, **kwargs):
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
            temp_dif = temperature_max - temperature_min
            tmp = (2 * temp_dif) * (temperature_min - temperature_basal_inferior)
            tmp1 = math.pow(temp_dif, 2)
            tmp2 = math.pow((temperature_max - temperature_basal_superior), 2)
            tmp3 = 2 * temp_dif
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


def get_days_passed_on_current_year():
    """
    Calculate the days passed until today date
    :return: int
    """
    today = date.today()
    return today.timetuple().tm_yday


def get_sunrise_sunset_hour(_sunset_hour_angle):
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


def get_vapour_pressure_deficit(es, ea):
    """
    VPD - Vapour Pressure Deficit
    :param es:
    :param ea: calculated by vp_from_rhmin_rhmax()
    :return:
    """
    return es - ea


def get_atm_pressure(altitude) -> float:
    """
    Estimate atmospheric pressure from altitude.

    Calculated using a simplification of the ideal gas law, assuming 20 degrees
    Celsius for a standard atmosphere. Based on equation 7, page 62 in Allen
    et al. (1998).

    :param altitude: Elevation/altitude above sea level [m]
    :return: atmospheric pressure [kPa]
    :rtype: float
    """
    temp = (293.0 - (0.0065 * altitude)) / 293.0
    return math.pow(temp, 5.26) * 101.3


def get_avp_from_temperature_min(temperature_min):
    """
    Estimate actual vapour pressure (*ea*) from minimum temperature.

    This method is to be used where humidity data are lacking or are of
    questionable quality. The method assumes that the dewpoint temperature
    is approximately equal to the minimum temperature (*temperature_min*), i.e. the
    air is saturated with water vapour at *temperature_min*.

    **Note**: This assumption may not hold in arid/semi-arid areas.
    In these areas it may be better to subtract 2 deg C from the
    minimum temperature (see Annex 6 in FAO paper).

    Based on equation 48 in Allen et al. (1998).

    :param temperature_min: Daily minimum temperature [deg C]
    :return: Actual vapour pressure [kPa]
    :rtype: float
    """
    return 0.611 * math.exp((17.27 * temperature_min) / (temperature_min + 237.3))


def get_avp_from_rhmin_rhmax(svp_temperature_min, svp_temperature_max, relative_humidity_min, rh_max) -> float:
    """
    Estimate actual vapour pressure (*ea*) from saturation vapour pressure and
    relative humidity.

    Based on FAO equation 17 in Allen et al. (1998).

    :param svp_temperature_min: Saturation vapour pressure at daily minimum temperature
        [kPa]. Can be estimated using ``svp_from_t()``.
    :param svp_temperature_max: Saturation vapour pressure at daily maximum temperature
        [kPa]. Can be estimated using ``svp_from_t()``.
    :param relative_humidity_min: Minimum relative humidity [%]
    :param rh_max: Maximum relative humidity [%]
    :return: Actual vapour pressure [kPa]
    :rtype: float
    """
    tmp1 = svp_temperature_min * (rh_max / 100.0)
    tmp2 = svp_temperature_max * (relative_humidity_min / 100.0)
    return (tmp1 + tmp2) / 2.0


def get_avp_from_rhmax(svp_temperature_min, relative_humidity_max):
    """
    Estimate actual vapour pressure (*e*a) from saturation vapour pressure at
    daily minimum temperature and maximum relative humidity

    Based on FAO equation 18 in Allen et al (1998).

    :param svp_temperature_min: Saturation vapour pressure at daily minimum temperature
        [kPa]. Can be estimated using ``svp_from_t()``.
    :param relative_humidity_max: Maximum relative humidity [%]
    :return: Actual vapour pressure [kPa]
    :rtype: float
    """
    return svp_temperature_min * (relative_humidity_max / 100.0)


def get_avp_from_rhmean(svp_temperature_min, svp_temperature_max, relative_humidity_mean):
    """
    Estimate actual vapour pressure (*ea*) from saturation vapour pressure at
    daily minimum and maximum temperature, and mean relative humidity.

    Based on FAO equation 19 in Allen et al (1998).

    :param svp_temperature_min: Saturation vapour pressure at daily minimum temperature
        [kPa]. Can be estimated using ``svp_from_t()``.
    :param svp_temperature_max: Saturation vapour pressure at daily maximum temperature
        [kPa]. Can be estimated using ``svp_from_t()``.
    :param relative_humidity_mean: Mean relative humidity [%] (average of RH min and RH max).
    :return: Actual vapour pressure [kPa]
    :rtype: float
    """
    return (relative_humidity_mean / 100.0) * ((svp_temperature_max + svp_temperature_min) / 2.0)


def get_avp_from_tdew(tdew):
    """
    Estimate actual vapour pressure (*ea*) from dewpoint temperature.

    Based on equation 14 in Allen et al (1998). As the dewpoint temperature is
    the temperature to which air needs to be cooled to make it saturated, the
    actual vapour pressure is the saturation vapour pressure at the dewpoint
    temperature.

    This method is preferable to calculating vapour pressure from
    minimum temperature.

    :param tdew: Dewpoint temperature [deg C]
    :return: Actual vapour pressure [kPa]
    :rtype: float
    """
    return 0.6108 * math.exp((17.27 * tdew) / (tdew + 237.3))


def get_avp_from_twet_tdry(twet, tdry, svp_twet, psyc_const):
    """
    Estimate actual vapour pressure (*ea*) from wet and dry bulb temperature.

    Based on equation 15 in Allen et al (1998). As the dewpoint temperature
    is the temperature to which air needs to be cooled to make it saturated, the
    actual vapour pressure is the saturation vapour pressure at the dewpoint
    temperature.

    This method is preferable to calculating vapour pressure from
    minimum temperature.

    Values for the psychrometric constant of the psychrometer (*psy_const*)
    can be calculated using ``psyc_const_of_psychrometer()``.

    :param twet: Wet bulb temperature [deg C]
    :param tdry: Dry bulb temperature [deg C]
    :param svp_twet: Saturated vapour pressure at the wet bulb temperature
        [kPa]. Can be estimated using ``svp_from_t()``.
    :param psyc_const: Psychrometric constant of the pyschrometer [kPa deg C-1].
        Can be estimated using ``psy_const()`` or
        ``psy_const_of_psychrometer()``.
    :return: Actual vapour pressure [kPa]
    :rtype: float
    """
    return svp_twet - (psyc_const * (tdry - twet))


def get_clear_sky_radiation(altitude, ext_rad):
    """
    Estimate clear sky radiation from altitude and extraterrestrial radiation.

    Based on equation 37 in Allen et al (1998) which is recommended when
    calibrated Angstrom values are not available.

    :param altitude: Elevation above sea level [m]
    :param ext_rad: Extraterrestrial radiation [MJ m-2 time-1]. Can be
        estimated using ``et_rad()``.
    :return: Clear sky radiation [MJ m-2 time-1]
    :rtype: float
    """
    return (0.00002 * altitude + 0.75) * ext_rad


def get_daily_mean_temperature(temperature_min, temperature_max):
    """
    Estimate mean daily temperature from the daily minimum and maximum
    temperatures.

    :param temperature_min: Minimum daily temperature [deg C]
    :param temperature_max: Maximum daily temperature [deg C]
    :return: Mean daily temperature [deg C]
    :rtype: float
    """
    return (temperature_max + temperature_min) / 2.0


def get_daylight_hours(sunset_hour_angle):
    """
    Calculate daylight hours from sunset hour angle.

    Based on FAO equation 34 in Allen et al (1998).

    :param sunset_hour_angle: Sunset hour angle [rad]. Can be calculated using
        ``calculate_sunset_hour_angle()``.
    :return: Daylight hours.
    :rtype: float
    """
    check_sunset_hour_angle(sunset_hour_angle)
    return (24.0 / math.pi) * sunset_hour_angle


def get_delta_svp(t):
    """
    Estimate the slope of the saturation vapour pressure curve at a given
    temperature.

    Based on equation 13 in Allen et al. (1998). If using in the Penman-Monteith
    *t* should be the mean air temperature.

    :param t: Air temperature [deg C]. Use mean air temperature for use in
        Penman-Monteith.
    :return: Saturation vapour pressure [kPa degC-1]
    :rtype: float
    """
    tmp = 4098 * (0.6108 * math.exp((17.27 * t) / (t + 237.3)))
    return tmp / math.pow((t + 237.3), 2)


def get_latent_heat(temperature_mean):
    return (2500 - (2.37 * temperature_mean)) / 1000


def get_daily_extraterrestrial_radiation(latitude, solar_dec,
                                         sunset_hour_angle,
                                         inverse_relative_distance_earth_sun):
    """
    Estimate daily extraterrestrial radiation (*Ra*, 'top of the atmosphere
    radiation').

    Based on equation 21 in Allen et al. (1998). If monthly mean radiation is
    required make sure *sol_dec*. *sha* and *irl* have been calculated using
    the day of the year that corresponds to the middle of the month.

    **Note**: From Allen et al (1998): "For the winter months in latitudes
    greater than 55 degrees (N or S), the equations have limited validity.
    Reference should be made to the Smithsonian Tables to assess possible
    deviations."

    :param latitude: Latitude [radians]
    :param solar_dec: Solar declination [radians]. Can be calculated using
        ``calculate_solar_declination()``.
    :param sunset_hour_angle: Sunset hour angle [radians]. Can be calculated using
        ``calculate_sunset_hour_angle()``.
    :param inverse_relative_distance_earth_sun: Inverse relative distance earth-sun [dimensionless]. Can be
        calculated using ``inv_rel_dist_earth_sun()``.
    :return: Daily extraterrestrial radiation [MJ m-2 day-1]
    :rtype: float
    """
    check_latitude(latitude)
    check_solar_declination(solar_dec)
    check_sunset_hour_angle(sunset_hour_angle)

    tmp1 = (24.0 * 60.0) / math.pi
    tmp2 = sunset_hour_angle * math.sin(latitude) * math.sin(solar_dec)
    tmp3 = math.cos(latitude) * math.cos(solar_dec) * math.sin(sunset_hour_angle)
    return tmp1 * SOLAR_CONSTANT * inverse_relative_distance_earth_sun * (tmp2 + tmp3)


def get_15min_extraterrestrial_radiation(latitude, longitude, solar_dec, ird, day_of_year, register_hour):
    """
    Estimate 15 minutes extraterrestrial radiation (*Ra*, 'top of the atmosphere
    radiation').

    Based on equation 21 in Allen et al. (1998). If monthly mean radiation is
    required make sure *sol_dec*. *sha* and *irl* have been calculated using
    the day of the year that corresponds to the middle of the month.

    **Note**: From Allen et al (1998): "For the winter months in latitudes
    greater than 55 degrees (N or S), the equations have limited validity.
    Reference should be made to the Smithsonian Tables to assess possible
    deviations."

    :param latitude: Latitude [degree]
    :param longitude: longitude [degree]
    :param solar_dec: Solar declination [radians]. Can be calculated using
        ``sol_dec()``.
    :param ird: Inverse relative distance earth-sun [dimensionless]. Can be
        calculated using ``get_inverse_relative_distance_earth_sun()``.
    :param day_of_year: day of year of the register
    :param register_hour: time object
    :return: Daily extraterrestrial radiation [MJ m-2 day-1]
    :rtype: float
    """
    check_latitude(latitude)
    check_longitude(longitude)
    check_solar_declination(solar_dec)

    latitude = deg2rad(latitude)

    # b
    b = (2 * math.pi * (day_of_year - 81)) / 364

    # Seasonal Solar Time
    sc = 0.1645 * math.sin(2 * b) - 0.1255 * math.cos(b) - 0.025 * math.sin(b)

    # calculate t
    if not isinstance(register_hour, time):
        raise ValueError('Register hour is not a datetime object')

    t = register_hour.hour + (register_hour.minute / 60)

    # Solar Time Angle
    w = math.pi / 12 * ((t + 0.06667 * (REFERENTIAL_LONGITUDE - longitude) + sc) - 12)

    # w1
    w1 = w - ((math.pi * 0.25) / 24)

    # w2
    w2 = w + ((math.pi * 0.25) / 24)

    tmp1 = ((12.0 * 60.0) / math.pi) * SOLAR_CONSTANT * ird
    tmp2 = (w2 - w1) * math.sin(latitude) * math.sin(solar_dec)
    tmp3 = math.cos(latitude) * math.cos(solar_dec) * (math.sin(w2) - math.sin(w1))

    et_rad_value = tmp1 * (tmp2 + tmp3)

    if et_rad_value < 0:
        return 0

    return et_rad_value


def get_reference_evapotranspiration_over_grass(temperature_min, temperature_max, temperature_mean, et_radiation):
    """
    Estimate reference evapotranspiration over grass (ETo) using the Hargreaves
    equation.

    Generally, when solar radiation data, relative humidity data
    and/or wind speed data are missing, it is better to estimate them using
    the functions available in this module, and then calculate ETo
    the FAO Penman-Monteith equation. However, as an alternative, ETo can be
    estimated using the Hargreaves ETo equation.

    Based on equation 52 in Allen et al. (1998).

    :param temperature_min: Minimum daily temperature [deg C]
    :param temperature_max: Maximum daily temperature [deg C]
    :param temperature_mean: Mean daily temperature [deg C]. If emasurements not
        available it can be estimated as (*temperature_min* + *temperature_max*) / 2.
    :param et_radiation: Extraterrestrial radiation (Ra) [MJ m-2 day-1]. Can be
        estimated using ``et_rad()``.
    :return: Reference evapotranspiration over grass (ETo) [mm day-1]
    :rtype: float
    """
    # Note, multiplied by 0.408 to convert extraterrestrial radiation could
    # be given in MJ m-2 day-1 rather than as equivalent evaporation in
    # mm day-1
    return 0.0023 * (temperature_mean + 17.8) * (temperature_max - temperature_min) ** 0.5 * 0.408 * et_radiation


def get_inverse_relative_distance_earth_sun(day_of_year):
    """
    Calculate the inverse relative distance between earth and sun from
    day of the year.

    Based on FAO equation 23 in Allen et al (1998).

    :param day_of_year: Day of the year [1 to 366]
    :return: Inverse relative distance between earth and the sun
    :rtype: float
    """
    check_day_of_year(day_of_year)
    return 1 + (0.033 * math.cos((2.0 * math.pi / 365.0) * day_of_year))


def get_mean_saturation_vapor_pressure(temperature_min, temperature_max):
    """
    Estimate mean saturation vapour pressure, *es* [kPa] from minimum and
    maximum temperature.

    Based on equations 11 and 12 in Allen et al (1998).

    Mean saturation vapour pressure is calculated as the mean of the
    saturation vapour pressure at temperature_max (maximum temperature) and temperature_min
    (minimum temperature).

    :param temperature_min: Minimum temperature [deg C]
    :param temperature_max: Maximum temperature [deg C]
    :return: Mean saturation vapour pressure (*es*) [kPa]
    :rtype: float
    """
    return (get_svp_from_temp(temperature_min) + get_svp_from_temp(temperature_max)) / 2.0


# not used at the moment
def get_monthly_soil_heat_flux(t_month_prev, t_month_next):
    """
    Estimate monthly soil heat flux (Gmonth) from the mean air temperature of
    the previous and next month, assuming a grass crop.

    Based on equation 43 in Allen et al. (1998). If the air temperature of the
    next month is not known use ``monthly_soil_heat_flux2()`` instead. The
    resulting heat flux can be converted to equivalent evaporation [mm day-1]
    using ``energy2evap()``.

    :param t_month_prev: Mean air temperature of the previous month
        [deg Celsius]
    :param t_month_next: Mean air temperature of the next month [deg Celsius]
    :return: Monthly soil heat flux (Gmonth) [MJ m-2 day-1]
    :rtype: float
    """
    return 0.07 * (t_month_next - t_month_prev)


def get_monthly_soil_heat_flux2(previous_month_temperature, current_month_temperature):
    """
    Estimate monthly soil heat flux (Gmonth) [MJ m-2 day-1] from the mean
    air temperature of the previous and current month, assuming a grass crop.

    Based on equation 44 in Allen et al. (1998). If the air temperature of the
    next month is available, use ``monthly_soil_heat_flux()`` instead. The
    resulting heat flux can be converted to equivalent evaporation [mm day-1]
    using ``energy2evap()``.

    :param previous_month_temperature: Mean air temperature of the previous month
        [deg Celsius]
    :param current_month_temperature: Mean air temperature of the current month [deg Celsius]
    :return: Monthly soil heat flux (Gmonth) [MJ m-2 day-1]
    :rtype: float
    """
    return 0.14 * (current_month_temperature - previous_month_temperature)


def get_soil_heat_flux_by_night_or_day_period(net_radiation, is_day=True) -> float:
    """
    Calculate soil heat flux (G) [MJ m-2 15min] from net radiation by night or day period.

    :param net_radiation: Net radiation at crop surface
    :param is_day: variable to say if is day period or night period
    :return: Shorten period heat flux [MJ m-2 15min]
    :rtype: float
    """
    if is_day:
        soil_heat_flux = 0.1 * net_radiation
    else:
        soil_heat_flux = 0.5 * net_radiation

    return soil_heat_flux


def get_net_in_sol_rad(sol_rad, albedo=0.23):
    """
    Calculate net incoming solar (or shortwave) radiation from gross
    incoming solar radiation, assuming a grass reference crop.

    Net incoming solar radiation is the net shortwave radiation resulting
    from the balance between incoming and reflected solar radiation. The
    output can be converted to equivalent evaporation [mm day-1] using
    ``energy2evap()``.

    Based on FAO equation 38 in Allen et al. (1998).

    :param sol_rad: Gross incoming solar radiation [MJ m-2 day-1]. If
        necessary this can be estimated using functions whose name
        begins with 'sol_rad_from'.
    :param albedo: Albedo of the crop as the proportion of gross incoming solar
        radiation that is reflected by the surface. Default value is 0.23,
        which is the value used by the FAO for a short grass reference crop.
        Albedo can be as high as 0.95 for freshly fallen snow and as low as
        0.05 for wet bare soil. A green vegetation over has an albedo of
        about 0.20-0.25 (Allen et al, 1998).
    :return: Net incoming solar (or shortwave) radiation [MJ m-2 day-1].
    :rtype: float
    """
    return (1 - albedo) * sol_rad


def get_net_out_lw_rad(temperature_min, temperature_max, sol_rad, cs_radiation, avp, time_period="15min"):
    """
    Estimate net outgoing longwave radiation.

    This is the net longwave energy (net energy flux) leaving the
    earth's surface. It is proportional to the absolute temperature of
    the surface raised to the fourth power according to the Stefan-Boltzmann
    law. However, water vapour, clouds, carbon dioxide and dust are absorbers
    and emitters of long wave radiation. This function corrects the Stefan-Boltzmann
    law for humidity (using actual vapor pressure) and cloudiness
    (using solar radiation and clear sky radiation). The concentrations of all
    other absorbers are assumed to be constant.

    The output can be converted to equivalent evaporation [mm day-1] using
    ``energy2evap()``.

    Based on FAO equation 39 in Allen et al (1998).

    :param temperature_min: Absolute daily minimum temperature [degrees Kelvin]
    :param temperature_max: Absolute daily maximum temperature [degrees Kelvin]
    :param sol_rad: Solar radiation [MJ m-2 day-1]. If necessary this can be
        estimated using ``sol+rad()``.
    :param cs_radiation: Clear sky radiation [MJ m-2 day-1]. Can be estimated using
        ``cs_rad()``.
    :param avp: Actual vapour pressure [kPa]. Can be estimated using functions
        with names beginning with 'avp_from'.
    :param time_period The period of time that will be used to calculate the result
        ( Supported values: daily, hourly, half_hourly and 15min )
    :return: Net outgoing longwave radiation [MJ m-2 day-1]
    :rtype: float
    """
    if time_period == "15min":
        stefan_boltzmann = STEFAN_BOLTZMANN_CONSTANT / 96
    elif time_period == 'half_hourly':
        stefan_boltzmann = STEFAN_BOLTZMANN_CONSTANT / 48
    elif time_period == 'hourly':
        stefan_boltzmann = STEFAN_BOLTZMANN_CONSTANT / 24
    else:
        stefan_boltzmann = STEFAN_BOLTZMANN_CONSTANT

    if cs_radiation > 0 and sol_rad > 0:  # Fix to avoid division by Zero
        rad = sol_rad / cs_radiation
        if rad > 1:
            rad = 1.0
    else:
        rad = 0.5

    tmp1 = (stefan_boltzmann * ((math.pow(temperature_max, 4) + math.pow(temperature_min, 4)) / 2))
    tmp2 = (0.34 - (0.14 * math.sqrt(avp)))
    tmp3 = 1.35 * rad - 0.35

    # Apply constant if needed
    # constant = 1.35 * 0.5 - 0.35 = 0.325
    tmp3 = tmp3 if tmp3 < 0 else 0.325

    return tmp1 * tmp2 * tmp3


def get_net_rad(ni_sw_rad, no_lw_rad):
    """
    Calculate net radiation at the crop surface by time, assuming a grass
    reference crop.

    Net radiation is the difference between the incoming net shortwave (or
    solar) radiation and the outgoing net longwave radiation. Output can be
    converted to equivalent evaporation [mm time-1] using ``energy2evap()``.

    Based on equation 40 in Allen et al (1998).

    :param ni_sw_rad: Net incoming shortwave radiation [MJ m-2 time-1]. Can be
        estimated using ``net_in_sol_rad()``.
    :param no_lw_rad: Net outgoing longwave radiation [MJ m-2 time-1]. Can be
        estimated using ``net_out_lw_rad()``.
    :return: Daily net radiation [MJ m-2 time-1].
    :rtype: float
    """
    return ni_sw_rad - no_lw_rad


def get_psy_const(atmosphere_pressure, latent_ht):
    """
    Calculate the psychrometric constant.

    This method assumes that the air is saturated with water vapour at the
    minimum daily temperature. This assumption may not hold in arid areas.

    Based on equation 8, page 95 in Allen et al (1998).

    :param atmosphere_pressure: Atmospheric pressure [kPa]. Can be estimated using
        ``atm_pressure()``.
    :param latent_ht: Can be estimated using latent_heat()
    :return: Psychrometric constant [kPa degC-1].
    :rtype: float
    """
    cp = 0.001013
    ratio_molecular_weight_water_dryair = 0.622

    return (cp * atmosphere_pressure) / (ratio_molecular_weight_water_dryair * latent_ht)


def get_psy_const_of_psychrometer(psychrometer, atmos_pres):
    """
    Calculate the psychrometric constant for different types of
    psychrometer at a given atmospheric pressure.

    Based on FAO equation 16 in Allen et al (1998).

    :param psychrometer: Integer between 1 and 3 which denotes type of
        psychrometer:
        1. ventilated (Asmann or aspirated type) psychrometer with
           an air movement of approximately 5 m/s
        2. natural ventilated psychrometer with an air movement
           of approximately 1 m/s
        3. non ventilated psychrometer installed indoors
    :param atmos_pres: Atmospheric pressure [kPa]. Can be estimated using
        ``atm_pressure()``.
    :return: Psychometric constant [kPa degC-1].
    :rtype: float
    """
    # Select coefficient based on type of ventilation of the wet bulb
    if psychrometer == 1:
        psy_coeff = 0.000662
    elif psychrometer == 2:
        psy_coeff = 0.000800
    elif psychrometer == 3:
        psy_coeff = 0.001200
    else:
        raise ValueError(
            'psychrometer should be in range 1 to 3: {0!r}'.format(psychrometer))

    return psy_coeff * atmos_pres


def get_rh_from_avp_svp(avp, sat_vp):
    """
    Calculate relative humidity as the ratio of actual vapour pressure
    to saturation vapour pressure at the same temperature.

    See Allen et al (1998), page 67 for details.

    :param avp: Actual vapour pressure [units do not matter so long as they
        are the same as for *svp*]. Can be estimated using functions whose
        name begins with 'avp_from'.
    :param sat_vp: Saturated vapour pressure [units do not matter so long as they
        are the same as for *avp*]. Can be estimated using ``svp_from_t()``.
    :return: Relative humidity [%].
    :rtype: float
    """
    return 100.0 * avp / sat_vp


def get_solar_declination(day_of_year):
    """
    Calculate solar declination from day of the year.

    Based on FAO equation 24 in Allen et al (1998).

    :param day_of_year: Day of year integer between 1 and 365 or 366.
    :return: solar declination [radians]
    :rtype: float
    """
    check_day_of_year(day_of_year)
    return 0.409 * math.sin(((2.0 * math.pi / 365.0) * day_of_year - 1.39))


def get_sol_rad_from_sun_hours(dl_hours, sunshine_hours, et_radiation):
    """
    Calculate incoming solar (or shortwave) radiation, *Rs* (radiation hitting
    a horizontal plane after scattering by the atmosphere) from relative
    sunshine duration.

    If measured radiation data are not available this method is preferable
    to calculating solar radiation from temperature. If a monthly mean is
    required then divide the monthly number of sunshine hours by number of
    days in the month and ensure that *et_rad* and *daylight_hours* was
    calculated using the day of the year that corresponds to the middle of
    the month.

    Based on equations 34 and 35 in Allen et al (1998).

    :param dl_hours: Number of daylight hours [hours]. Can be calculated
        using ``daylight_hours()``.
    :param sunshine_hours: Sunshine duration [hours]. Can be calculated
        using ``sunshine_hours()``.
    :param et_radiation: Extraterrestrial radiation [MJ m-2 day-1]. Can be
        estimated using ``et_rad()``.
    :return: Incoming solar (or shortwave) radiation [MJ m-2 day-1]
    :rtype: float
    """
    check_day_hours(sunshine_hours, 'sun_hours')
    check_day_hours(dl_hours, 'daylight_hours')

    # 0.5 and 0.25 are default values of regression constants (Angstrom values)
    # recommended by FAO when calibrated values are unavailable.
    return (0.5 * sunshine_hours / dl_hours + 0.25) * et_radiation


def get_solar_radiation_from_temperature(et_radiation, cs_radiation, temperature_min, temperature_max, coastal):
    """
    Estimate incoming solar (or shortwave) radiation, *Rs*, (radiation hitting
    a horizontal plane after scattering by the atmosphere) from min and max
    temperature together with an empirical adjustment coefficient for
    'interior' and 'coastal' regions.

    The formula is based on equation 50 in Allen et al (1998) which is the
    Hargreaves radiation formula (Hargreaves and Samani, 1982, 1985). This
    method should be used only when solar radiation or sunshine hours data are
    not available. It is only recommended for locations where it is not
    possible to use radiation data from a regional station (either because
    climate conditions are heterogeneous or data are lacking).

    **NOTE**: this method is not suitable for island locations due to the
    moderating effects of the surrounding water.

    :param et_radiation: Extraterrestrial radiation [MJ m-2 day-1]. Can be
        estimated using ``et_rad()``.
    :param cs_radiation: Clear sky radiation [MJ m-2 day-1]. Can be estimated
        using ``cs_rad()``.
    :param temperature_min: Daily minimum temperature [deg C].
    :param temperature_max: Daily maximum temperature [deg C].
    :param coastal: ``True`` if site is a coastal location, situated on or
        adjacent to coast of a large land mass and where air masses are
        influenced by a nearby water body, ``False`` if interior location
        where land mass dominates and air masses are not strongly influenced
        by a large water body.
    :return: Incoming solar (or shortwave) radiation (Rs) [MJ m-2 day-1].
    :rtype: float
    """
    # Determine value of adjustment coefficient [deg C-0.5] for
    # coastal/interior locations
    if coastal:
        adj = 0.19
    else:
        adj = 0.16

    sol_rad = adj * math.sqrt(temperature_max - temperature_min) * et_radiation

    # The solar radiation value is constrained by the clear sky radiation
    return min(sol_rad, cs_radiation)


def get_solar_radiation_island(et_radiation):
    """
    Estimate incoming solar (or shortwave) radiation, *Rs* (radiation hitting
    a horizontal plane after scattering by the atmosphere) for an island
    location.

    An island is defined as a land mass with width perpendicular to the
    coastline <= 20 km. Use this method only if radiation data from
    elsewhere on the island is not available.

    **NOTE**: This method is only applicable for low altitudes (0-100 m)
    and monthly calculations.

    Based on FAO equation 51 in Allen et al (1998).

    :param et_radiation: Extraterrestrial radiation [MJ m-2 day-1]. Can be
        estimated using ``et_rad()``.
    :return: Incoming solar (or shortwave) radiation [MJ m-2 day-1].
    :rtype: float
    """
    return (0.7 * et_radiation) - 4.0


def get_sunset_hour_angle(latitude, solar_declination):
    """
    Calculate sunset hour angle (*Ws*) from latitude and solar
    declination.

    Based on FAO equation 25 in Allen et al (1998).

    :param latitude: Latitude [radians]. Note: *latitude* should be negative
        if it in the southern hemisphere, positive if in the northern
        hemisphere.
    :param solar_declination: Solar declination [radians]. Can be calculated using
        ``sol_dec()``.
    :return: Sunset hour angle [radians].
    :rtype: float
    """
    check_latitude(latitude)
    check_solar_declination(solar_declination)

    cos_sha = -math.tan(latitude) * math.tan(solar_declination)
    # If tmp is >= 1 there is no sunset, i.e. 24 hours of daylight
    # If tmp is <= 1 there is no sunrise, i.e. 24 hours of darkness
    # See https://www.itacanet.org/the-sun-as-a-source-of-energy/
    # part-3-calculating-solar-angles/
    # Domain of acos is -1 <= x <= 1 radians (this is not mentioned in FAO-56!)
    return math.acos(min(max(cos_sha, -1.0), 1.0))


def get_svp_from_temp(t):
    """
    Estimate saturation vapour pressure (*es*) from air temperature.

    Based on equations 11 and 12 in Allen et al (1998).

    :param t: Temperature [deg C]
    :return: Saturation vapour pressure [kPa]
    :rtype: float
    """
    return 0.6108 * math.exp((17.27 * t) / (t + 237.3))


def get_svp(svp_min, svp_max):
    return (svp_max + svp_min) / 2


def get_wind_speed_2m(ws, z):
    """
    Convert wind speed measured at different heights above the soil
    surface to wind speed at 2 m above the surface, assuming a short grass
    surface.

    Based on FAO equation 47 in Allen et al (1998).

    :param ws: Measured wind speed [m s-1]
    :param z: Height of wind measurement above ground surface [m]
    :return: Wind speed at 2 m above the surface [m s-1]
    :rtype: float
    """
    return ws * (4.87 / math.log((67.8 * z) - 5.42))


def get_monthly_mean_daylight_hours(latitude, year=None):
    """
    Calculate mean daylight hours for each month of the year for a given
    latitude.

    :param latitude: Latitude [radians]
    :param year: Year for the daylight hours are required. The only effect of
        *year* is to change the number of days in Feb to 29 if it is a leap
        year. If left as the default, None, then a normal (non-leap) year is
        assumed.
    :return: Mean daily daylight hours of each month of a year [hours]
    :rtype: List of floats.
    """
    check_latitude(latitude)

    if year is None or not calendar.isleap(year):
        month_days = _MONTHDAYS
    else:
        month_days = _LEAP_MONTHDAYS
    monthly_mean_dlh = []
    doy = 1  # Day of the year
    for mdays in month_days:
        dlh = 0.0  # Cumulative daylight hours for the month
        for daynum in range(1, mdays + 1):
            sd = get_solar_declination(doy)
            sha = get_sunset_hour_angle(latitude, sd)
            dlh += get_daylight_hours(sha)
            doy += 1
        # Calc mean daylight hours of the month
        monthly_mean_dlh.append(dlh / mdays)
    return monthly_mean_dlh
