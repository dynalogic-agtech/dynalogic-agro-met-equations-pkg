import calendar
import math
from AgroMetEquations.config import _MONTHDAYS, _LEAP_MONTHDAYS


# #############################################################################
#                         FAO-56 Penman-Monteith                              #
# #############################################################################
def fao56_penman_monteith(net_radiation, temperature_mean, ws, latent_ht, sat_vp, avp,
                          delta_sat_vp, psy, sol_rad, shf=0.0, time_period="15min"):
    """
    Estimate reference evapotranspiration (ETo) from a hypothetical
    short grass reference surface using the FAO-56 Penman-Monteith equation.

    Based on equation 6 in Allen et al. (1998).

    :param net_radiation: Net radiation at crop surface [MJ m-2 day-1]. If
        necessary this can be estimated using ``net_rad()``.
    :param temperature_mean: Air temperature at 2 m height [deg Celsius].
    :param ws: Wind speed at 2 m height [m s-1]. If not measured at 2m,
        convert using ``wind_speed_at_2m()``.
    :param latent_ht: Latent heat Can be calculated using ``latent_heat(temperature_mean)``.
    :param sat_vp: Saturation vapour pressure [kPa]. Can be estimated using
        ``svp_from_t()''.
    :param avp: Actual vapour pressure [kPa]. Can be estimated using a range
        of functions with names beginning with 'avp_from'.
    :param delta_sat_vp: Slope of saturation vapour pressure curve [kPa degC-1].
        Can be estimated using ``delta_svp()``.
    :param psy: Psychometric constant [kPa deg C]. Can be estimated using
        ``psy_const_of_psychrometer()`` or ``psy_const()``.
    :param sol_rad: Solar Radiation to calculate the day and night period
    :param shf: Soil heat flux (G) [MJ m-2 day-1] (default is 0.0, which is
        reasonable for a daily or 10-day time steps). For monthly time steps
        *shf* can be estimated using ``monthly_soil_heat_flux()`` or
        ``monthly_soil_heat_flux2()``.
    :param time_period The period of time that will be used to calculate the result
        ( Supported values: daily, hourly, half_hourly and 15min )
    :return: Reference evapotranspiration (ETo) from a hypothetical
        grass reference surface [mm day-1].
    :rtype: float
    """

    # time period conversion
    if time_period == "daily":
        time_period_conversion = 900
    elif time_period == "hourly":
        time_period_conversion = 37.5
    elif time_period == "half_hourly":
        time_period_conversion = 18.75
    else:
        time_period_conversion = 9.375  # 15min period

    cd = 0.24 if sol_rad > 1 else 0.96

    a1 = 1 / latent_ht
    a2 = (net_radiation - shf) * delta_sat_vp
    a3 = (time_period_conversion / (temperature_mean + 273)) * psy * ws * (sat_vp - avp)
    a4 = delta_sat_vp + (psy * (1 + cd * ws))

    return (a1 * a2 + a3) / a4


# #############################################################################
#                          Priestley-Taylor                                   #
# #############################################################################
def priestley_taylor(latent_ht, delta_sat_vp, psy, net_radiation, soil_heat_flux):
    """
    Estimate reference evapotranspiration (ETo) from a hypothetical
    :param latent_ht:
    :param delta_sat_vp:
    :param psy:
    :param net_radiation:
    :param soil_heat_flux:
    :return:
    """

    a1 = (1 / latent_ht) * 1.26
    a2 = delta_sat_vp / (delta_sat_vp + psy)
    a3 = (net_radiation - soil_heat_flux)

    return a1 * a2 * a3


# #############################################################################
#                 Hargreaves & Samani – With Solar Radiation                   #
# #############################################################################
def hargreaves_samani_with_solar_ratiation(latent_ht, temperature_mean, solar_radiation):
    a1 = (1 / latent_ht) * 0.0135
    a2 = (17.8 + temperature_mean) * solar_radiation

    return a1 * a2


# #############################################################################
#                 Hargreaves & Samani – Without Solar Radiation                   #
# #############################################################################
def hargreaves_samani_without_solar_ratiation(latent_ht, temperature_mean, temperature_max, temperature_min,
                                              extraterrestrial_radiation):
    a1 = (1 / latent_ht) * 0.0023
    a2 = temperature_mean + 17.8
    a3 = math.sqrt(temperature_max - temperature_min) * extraterrestrial_radiation

    return a1 * a2 * a3

# #############################################################################
#                     Thornthwaite (1948 method)                              #
# #############################################################################
"""
Calculate potential evapotranspiration using the Thornthwaite (1948 method)

:copyright: (c) 2015 by Mark Richards.
:license: BSD 3-Clause, see LICENSE.txt for more details.

References
----------

"""


def thornthwaite(monthly_t, monthly_mean_dlh, year=None):
    """
    Calculate potential evapotranspiration using the Thornthwaite CW (1948) An approach toward a rational
    classification of climate. Geographical Review, 38, 55-94.

    Estimate monthly potential evapotranspiration (PET) using the Thornthwaite (1948) method.

    Thornthwaite equation:
        *PET* = 1.6 (*L*/12) (*N*/30) (10*Ta* / *I*)***a*

    where:

    * *Ta* is the mean daily air temperature [deg C, if negative use 0] of the
      month being calculated
    * *N* is the number of days in the month being calculated
    * *L* is the mean day length [hours] of the month being calculated
    * *a* = (6.75 x 10-7)*I***3 - (7.71 x 10-5)*I***2 + (1.792 x 10-2)*I* + 0.49239
    * *I* is a heat index which depends on the 12 monthly mean temperatures and
      is calculated as the sum of (*Tai* / 5)**1.514 for each month, where
      Tai is the air temperature for each month in the year

    :param monthly_t: Iterable containing mean daily air temperature for each
        month of the year [deg C].
    :param monthly_mean_dlh: Iterable containing mean daily daylight
        hours for each month of the year [hours]. These can be calculated
        using ``monthly_mean_daylight_hours()``.
    :param year: Year for which PET is required. The only effect of year is
        to change the number of days in February to 29 if it is a leap year.
        If it is left as the default (None), then the year is assumed not to
        be a leap year.
    :return: Estimated monthly potential evaporation of each month of the year
        [mm/month]
    :rtype: List of floats
    """

    if len(monthly_t) != 12:
        raise ValueError(
            'monthly_t should be length 12 but is length {0}.'
            .format(len(monthly_t)))
    if len(monthly_mean_dlh) != 12:
        raise ValueError(
            'monthly_mean_dlh should be length 12 but is length {0}.'
            .format(len(monthly_mean_dlh)))

    if year is None or not calendar.isleap(year):
        month_days = _MONTHDAYS
    else:
        month_days = _LEAP_MONTHDAYS

    # Negative temperatures should be set to zero
    adj_monthly_t = [t * (t >= 0) for t in monthly_t]

    # Calculate the heat index (I)
    heat_index = 0.0
    for Tai in adj_monthly_t:
        if Tai / 5.0 > 0.0:
            heat_index += (Tai / 5.0) ** 1.514

    a = (6.75e-07 * heat_index ** 3) - (7.71e-05 * heat_index ** 2) + (1.792e-02 * heat_index) + 0.49239

    pet = []
    for Ta, L, N in zip(adj_monthly_t, monthly_mean_dlh, month_days):
        # Multiply by 10 to convert cm/month --> mm/month
        correction = (L / 12.0) * (N / 30.0)
        pet.append(
            1.6 * (L / 12.0) * (N / 30.0) * ((10.0 * Ta / heat_index) ** a) * 10.0)

    return pet
