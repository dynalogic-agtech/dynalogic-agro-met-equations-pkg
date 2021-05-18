Agrometeorological Equations Package
==========
##### Library for agrometeorological calculation like evapotranspiration, water balance, degree days, etc

|           |                                                    |
|-----------|----------------------------------------------------|
|Authors:   | Bruno Ducraux (bducraux@dynalogic.net)             |
|           | Mariana Gonçalves dos Reis (mreis@dynalogic.net)   |
| Based on: | [PyEto](https://github.com/woodcrafty/PyETo)       |
|           |                                                    |

About
-----------
This library was adapted to fit some needs of a private project, 
and I decided to disponibilize as open source since the original project 
that we used as base is.

All calculations and formulas were reviewed by the agrometeorologist Mariana Gonçalves dos Reis, based on the documents:

[Conversion factors and general equations applied in agricultural and forest meteorology](https://seer.sct.embrapa.br/index.php/agrometeoros/article/view/26527)

[Evapotranspiration_Equation.pdf](AgroMetEquations/docs/Evapotranspiration_Equation.pdf)

In case of questions related to the python code, bug fix ...
please contact Bruno Ducraux, who is the python developer responsible for the project.

Installation
------------
`pip install agro-met-equations-dynalogic`

Usage
-----

`from AgroMetEquations import`

# Steps to calculate FAO56

1. Day of year: days_passed_on_current_year()
2. Solar declination: sol_dec(day_of_year)
3. Sunset hour angle: sunset_hour_angle(latitude_rad, solar_declination)
4. Daylight hours: daylight_hours(sunset_hour_angle)
5. Inverse relative distance between earth and sun: inv_rel_dist_earth_sun(day_of_year)
6. Extraterrestrial radiation: et_rad(latitude_rad, solar_declination, sunset_hour_angle, inv_rel_dist_earth_sun)
7. Clear sky radiation: cs_rad(altitude, extraterrestrial_radiation)
8. Saturation vapour pressure Min: svp_from_t(temperature_min)
9. Saturation vapour pressure Max: svp_from_t(temperature_max)
10. Saturation vapour pressure: svp(svp_min, svp_max)
11. Actual vapour pressure: avp_from_rhmin_rhmax(svp_min, svp_max, relative_humidity_min, relative_humidity_max)
12. Net outgoing longwave radiation: net_out_lw_rad(temperature_min:kelvin, temperature_max:kelvin, solar_radiation, clear_sky_radiation, actual_vapour_pressure)
13. Net income solar radiation: net_in_sol_rad(solar_radiation)
14. Net radiation at the crop surface: net_rad(net_in_sol_rad, net_outgoing_longwave_radiation)
15. Soil heat flux: soil_heat_flux_by_nightday_period(net_rad)
16. Latent heat: latent_heat(temperature_mean)
17. Delta: delta_svp(temperature_mean)
18. Psychrometric constant: psy_const(atmosphere_pressure, latent_heat)
19. Wind speed measured at different heights: wind_speed_2m(wind_speed, sensor_height)
20. FAO56: fao56_penman_monteith(net_rad, temperature_mean, wind_speed, latent_heat, svp, actual_vapour_pressure, delta, psychrometric_constant, soil_heat_flux)


Testing
-------
To test the code you need to have pytest installed.

`pip install pytest`

Inside the AgroMetEquations folder run the command:

`pytest`
