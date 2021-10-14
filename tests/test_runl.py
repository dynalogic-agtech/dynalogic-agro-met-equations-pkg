from AgroMetEquations.convert import *
from AgroMetEquations.general import *
from AgroMetEquations.fao import *
from AgroMetEquations.thornthwaite import *

latitude_degree = 20
latitude_minute = 45
latitude_second = 0
latitude_direction = 'S'

longitude_degree = 47
longitude_minute = 38
longitude_second = 00
longitude_direction = 'W'

# sensor variables
temperature_min = 20
temperature_max = 25
temperature_mean = 20

relative_humidity_mean = 70

wind_speed = 3

# fixed variable by sensor
height_above_ground = 3

# end of sensors


latitude_dd = dms2dd(latitude_degree, latitude_minute, latitude_second, latitude_direction)
print("Latitude em graus:")
print(latitude_dd)

longitude_dd = dms2dd(longitude_degree, longitude_minute, longitude_second, longitude_direction)
print("Longitude em graus:")
print(longitude_dd)

latitude_rad = deg2rad(latitude_dd)
print("Latitude em radianos:")
print(latitude_rad)
longitute_rad = deg2rad(longitude_dd)
print("Longitute em radianos:")
print(longitute_rad)

altitude = 700
print("Altitude:")
print(altitude)

solar_radiation_15min = watt2mj15min(300)
print("Solar Radiation:")
print(solar_radiation_15min)

day_of_year = days_passed_on_current_year()

print("Day of year:")
print(day_of_year)

solar_declination = sol_dec(day_of_year)

print("Solar Declination rad:")
print(solar_declination)
print("Solar Declination deg:")
print(rad2deg(solar_declination))


sunset_hour_angle = sunset_hour_angle(latitude_rad, solar_declination)

print("Sunset Hour Angle in rad:")
print(sunset_hour_angle)
print("Sunset Hour Angle in degree:")
print(rad2deg(sunset_hour_angle))

sunrise_hour = sunrise_sunset_hour(rad2deg(sunset_hour_angle))
print("Sunrise / Sunset Hour:")
print(sunrise_hour)

daylight_h = daylight_hours(sunset_hour_angle)
print("Day light hours:")
print(daylight_h)

inv_rel_dist_earth_sun = inv_rel_dist_earth_sun(day_of_year)

print("Inverse relative distance between earth and sun:")
print(inv_rel_dist_earth_sun)

extraterrestrial_radiation_daily = et_rad(latitude_rad, solar_declination, sunset_hour_angle, inv_rel_dist_earth_sun)
print("Daily extraterrestrial radiation:")
print(extraterrestrial_radiation_daily)

extraterrestrial_radiation_15min = dailyto15min(extraterrestrial_radiation_daily)
print("15min extraterrestrial radiation:")
print(extraterrestrial_radiation_15min)

clear_sky_radiation = cs_rad(altitude, extraterrestrial_radiation_15min)
print("Clear sky radiation:")
print(clear_sky_radiation)

# es
# values to be extracted form sensors
tmax = 25.1
tmin = 19.1

svp_min = svp_from_t(tmin)
svp_max = svp_from_t(tmax)

svp = svp(svp_min, svp_max)

print("esMin:")
print(svp_min)

print("esMax:")
print(svp_max)

print("es:")
print(svp)


# values to be extracted from sensor
rh_min = 80
rh_max = 95
# ea
actual_vapour_pressure = avp_from_rhmin_rhmax(svp_min, svp_max, rh_min, rh_max)

print("ea")
print(actual_vapour_pressure)

# net outgoing longwave radiation

net_outgoing_longwave_radiation = net_out_lw_rad(celsius2kelvin(tmin), celsius2kelvin(tmax), solar_radiation_15min,
                                                 clear_sky_radiation, actual_vapour_pressure)
print("net_outgoing_longwave_radiation_daily")
print(net_outgoing_longwave_radiation)

print("net_outgoing_longwave_radiation_15min")
print(dailyto15min(net_outgoing_longwave_radiation))

# net incoming solar
net_in_sol_rad = net_in_sol_rad(solar_radiation_15min)
print("net_in_sol_rad:")
print(net_in_sol_rad)

# daily net radiation at the crop surface
net_rad = net_rad(net_in_sol_rad, dailyto15min(net_outgoing_longwave_radiation))
print("net_rad:")
print(net_rad)

# soil heat flux daylight
soil_heat_flux_daylight = soil_heat_flux_by_nightday_period(net_rad)
print("soil_heat_flux_daylight:")
print(soil_heat_flux_daylight)

# soil heat flux nightlight
soil_heat_flux_nightlight = soil_heat_flux_by_nightday_period(net_rad, False)
print("soil_heat_flux_nightlight:")
print(soil_heat_flux_nightlight)

# dew point temperature
# values to be extracted from sensor
dewpoint_temperature = dew_point(temperature_mean, relative_humidity_mean)
print("Dew Point: ")
print(dewpoint_temperature)

# frost point
frostpoint_temperature = frost_point(temperature_mean, dewpoint_temperature)
print("Frost Point: ")
print(frostpoint_temperature)

# latent heat
latent_heat = latent_heat(temperature_mean)
print("Latent Heat: ")
print(latent_heat)

# psychrometric constant
atmosphere_pressure = 100
psychrometric_constant = psy_const(atmosphere_pressure, latent_heat)
print("psychrometric constant")
print(psychrometric_constant)

# wind speed measured at different heights
wind_speed = wind_speed_2m(wind_speed, height_above_ground)
print("wind speed measured at different heights: ")
print(wind_speed)

# slope of the saturation vapour pressure curve
delta = delta_svp(temperature_mean)
print("slope of the saturation vapour pressure curve: ")
print(delta)

# fao56_penman_monteith
fao_56 = fao56_penman_monteith(net_rad, temperature_mean, wind_speed, latent_heat, svp, actual_vapour_pressure, delta,
                               psychrometric_constant, soil_heat_flux_daylight)
print("fao56_penman_monteith: ")
print(fao_56)


# degree day
test = degree_day(temperature_mean, 21, temperature_basal_superior=12, temperatura_max=25, temperature_min=11)


# Começo dos testes para cálculo do balanço hídrico

# Variaveis de entrada - Características físicas do solo
field_capacity = 35     # percent
permanent_wilting_point = 18    # percent
depth_of_soil = 0.20    # meters
latitude_degree = -20.75
latitude_rad = deg2rad(latitude_degree)

test_monthly_t = [22.1, 22.0, 21.4, 19.5, 17.1, 15.4, 15.0, 16.0, 18.2, 19.8, 20.7, 21.3]

mmdlh = monthly_mean_daylight_hours(latitude_rad)
print("monthly_mean_daylight_hours: ")
print(mmdlh)

pet = thornthwaite(test_monthly_t, mmdlh)
print("thornthwaite: ")
print(pet)
