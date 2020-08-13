import unittest
from AgroMetEquations import (sol_dec,
                              sunset_hour_angle,
                              daylight_hours,
                              inv_rel_dist_earth_sun,
                              et_rad,
                              dailyto15min,
                              cs_rad,
                              svp_from_t,
                              svp,
                              avp_from_rhmin_rhmax,
                              watt2mj15min,
                              net_out_lw_rad,
                              celsius2kelvin,
                              net_in_sol_rad,
                              net_rad,
                              soil_heat_flux_by_nightday_period,
                              latent_heat,
                              delta_svp,
                              psy_const,
                              wind_speed_2m,
                              fao56_penman_monteith
                              )


class TestFAO(unittest.TestCase):

    solar_declination = None
    sunset_h_angle = None
    inv_relative_distance_earth_sun = None
    extraterrestrial_radiation_15min = None
    svp_min = None
    svp_max = None
    clear_sky_radiation = None
    solar_radiation_15min = None
    net_income_solar_radiation = None
    net_outgoing_longwave_radiation_15min = None
    net_radiation = None
    wind_speed = None
    latent_heat = None
    svp = None
    actual_vapour_pressure = None
    delta = None
    gamma = None
    soil_heat_flux = None

    def setUp(self):
        # Sensor data that will come from the meteorologic station

        # temperature_min   - Minimum temperature
        # temperature_max   - Maximum temperature
        # mean_temp  - Mean temperature
        # relative_humidity_min     - Minimum relative humidity
        # relative_humidity_max     - Maximum relative humidity
        # mean_rh    - Mean relative humidity
        # solar_rad  - Solar radiation
        # atmosphere_pressure       - Atmospheric Pressure
        # wind_speed - Wind speed
        self.sensor_data = {'temperature_min': 36,
                            'temperature_max': 39,
                            'temperature_mean': 38,
                            'relative_humidity_min': 80,
                            'relative_humidity_max': 92,
                            'mean_rh': 85,
                            'solar_rad': 800,
                            'atmosphere_pressure': 95,
                            'wind_speed': 2.8}

        # Meteorologic Station Position Data
        # Latitude | Altitude | Sensor_Height
        self.station_position = {'latitude': -0.36215598,
                                 'altitude': 700,
                                 'sensor_height': 5}

        # Fixed Day of year
        self.doy = 84  # 25/03

    # Steps to calculate FAO56

    # 1 - Day of year: days_passed_on_current_year() -> For test purpose  this value is hardcoded to 25/03
    # 2 - Solar declination: sol_dec(day_of_year)
    # 3 - Sunset hour angle: sunset_hour_angle(latitude_rad, solar_declination)
    # 4 - Daylight hours: daylight_hours(sunset_hour_angle)
    # 5 - Inverse relative distance between earth and sun: inv_rel_dist_earth_sun(day_of_year)
    # 6 - Extraterrestrial radiation: et_rad(latitude_rad, solar_declination, sunset_hour_angle, inv_rel_dist_earth_sun)
    # 7 - Clear sky radiation: cs_rad(altitude, extraterrestrial_radiation)
    # 8 - Saturation vapour pressure Min: svp_from_t(temperature_min)
    # 9 - Saturation vapour pressure Max: svp_from_t(temperature_max)
    # 10 - Saturation vapour pressure: svp(svp_min, svp_max)
    # 11 - Actual vapour pressure: avp_from_rhmin_rhmax(svp_min, svp_max, relative_humidity_min, relative_humidity_max)
    # 12 - Net outgoing longwave radiation: net_out_lw_rad(temperature_min:kelvin, temperature_max:kelvin,
    #                                                      solar_radiation, clear_sky_radiation, actual_vapour_pressure)
    # 13 - Net income solar radiation: net_in_sol_rad(solar_radiation)
    # 14 - Net radiation at the crop surface: net_rad(net_in_sol_rad, net_outgoing_longwave_radiation)
    # 15 - Soil heat flux: soil_heat_flux_by_nightday_period(net_rad)
    # 16 - Latent heat: latent_heat(temperature_mean)
    # 17 - Delta: delta_svp(temperature_mean)
    # 18 - Psychrometric constant: psy_const(atmosphere_pressure, latent_heat)
    # 19 - Wind speed measured at different heights: wind_speed_2m(wind_speed, sensor_height)
    # 20 - FAO56: fao56_penman_monteith(net_rad, temperature_mean, wind_speed, latent_heat, svp, actual_vapour_pressure,
    #                                   delta, psychrometric_constant, soil_heat_flux)

    def test_01_solar_declination(self):
        expected = 0.022889
        TestFAO.solar_declination = sol_dec(self.doy)
        self.assertAlmostEqual(TestFAO.solar_declination, expected, delta=0.15)

    def test_02_sunset_hour_angle(self):
        expected = 1.562122837
        TestFAO.sunset_h_angle = sunset_hour_angle(self.station_position['latitude'], TestFAO.solar_declination)
        self.assertAlmostEqual(TestFAO.sunset_h_angle, expected, delta=0.15)

    def test_03_daylight_hours(self):
        expected = 11.93373942
        TestFAO.daylight_hrs = daylight_hours(TestFAO.sunset_h_angle)
        self.assertAlmostEqual(TestFAO.daylight_hrs, expected, delta=0.15)

    def test_04_inv_relative_distance_earth_sun(self):
        expected = 1.004107816
        TestFAO.inv_relative_distance_earth_sun = inv_rel_dist_earth_sun(self.doy)
        self.assertAlmostEqual(TestFAO.inv_relative_distance_earth_sun, expected, delta=0.15)

    def test_05_15min_extraterrestrial_radiation(self):
        expected = 0.3626729167
        daily_extraterrestrial_radiation = et_rad(self.station_position['latitude'],
                                                  TestFAO.solar_declination, TestFAO.sunset_h_angle,
                                                  TestFAO.inv_relative_distance_earth_sun)
        TestFAO.extraterrestrial_radiation_15min = dailyto15min(daily_extraterrestrial_radiation)
        self.assertAlmostEqual(TestFAO.extraterrestrial_radiation_15min, expected, delta=0.15)

    def test_06_claer_sky(self):
        expected = 0.2770
        TestFAO.clear_sky_radiation = cs_rad(self.station_position['altitude'],
                                             TestFAO.extraterrestrial_radiation_15min)
        self.assertAlmostEqual(TestFAO.clear_sky_radiation, expected, delta=0.15)

    def test_07_svp_min(self):
        expected = 5.940997702
        TestFAO.svp_min = svp_from_t(self.sensor_data['temperature_min'])
        self.assertAlmostEqual(TestFAO.svp_min, expected, delta=0.15)

    def test_08_svp_max(self):
        expected = 6.99146929
        TestFAO.svp_max = svp_from_t(self.sensor_data['temperature_max'])
        self.assertAlmostEqual(TestFAO.svp_max, expected, delta=0.15)

    def test_09_svp(self):
        expected = 6.466233496
        TestFAO.svp = svp(TestFAO.svp_min, TestFAO.svp_max)
        self.assertAlmostEqual(TestFAO.svp, expected, delta=0.15)

    def test_10_actual_vapour_pressure(self):
        expected = 5.529446659
        TestFAO.actual_vapour_pressure = avp_from_rhmin_rhmax(TestFAO.svp_min, TestFAO.svp_max,
                                                              self.sensor_data['relative_humidity_min'],
                                                              self.sensor_data['relative_humidity_max'])
        self.assertAlmostEqual(TestFAO.actual_vapour_pressure, expected, delta=0.15)

    def test_11_solar_radiation(self):
        expected = 0.72
        TestFAO.solar_radiation_15min = watt2mj15min(self.sensor_data['solar_rad'])
        self.assertEqual(TestFAO.solar_radiation_15min, expected)

    def test_12_net_outgoing_longwave_radiation(self):
        expected = 0.01621945
        kelvin_tmin = celsius2kelvin(self.sensor_data['temperature_min'])
        kelvin_tmax = celsius2kelvin(self.sensor_data['temperature_max'])
        net_outgoing_longwave_radiation_daily = net_out_lw_rad(kelvin_tmin, kelvin_tmax,
                                                               TestFAO.solar_radiation_15min,
                                                               TestFAO.clear_sky_radiation,
                                                               TestFAO.actual_vapour_pressure)
        TestFAO.net_outgoing_longwave_radiation_15min = dailyto15min(net_outgoing_longwave_radiation_daily)
        self.assertAlmostEqual(TestFAO.net_outgoing_longwave_radiation_15min, expected, delta=0.15)

    def test_13_net_income_solar_radiation(self):
        expected = 0.5544
        TestFAO.net_income_solar_radiation = net_in_sol_rad(TestFAO.solar_radiation_15min)
        self.assertAlmostEqual(TestFAO.net_income_solar_radiation, expected, delta=0.15)

    def test_14_daily_net_radiation(self):
        expected = 0.5382
        TestFAO.net_radiation = net_rad(TestFAO.net_income_solar_radiation,
                                        TestFAO.net_outgoing_longwave_radiation_15min)
        self.assertAlmostEqual(TestFAO.net_radiation, expected, delta=0.15)

    def test_15_soil_heat_flux(self):
        expected = 0.05382055
        TestFAO.soil_heat_flux = soil_heat_flux_by_nightday_period(TestFAO.net_radiation)
        self.assertAlmostEqual(TestFAO.soil_heat_flux, expected, delta=0.15)

    def test_16_latent_heat(self):
        expected = 2.40994
        TestFAO.latent_heat = latent_heat(self.sensor_data['temperature_mean'])
        self.assertAlmostEqual(TestFAO.latent_heat, expected, delta=0.15)

    def test_17_delta(self):
        expected = 0.358203
        TestFAO.delta = delta_svp(self.sensor_data['temperature_mean'])
        self.assertAlmostEqual(TestFAO.delta, expected, delta=0.15)

    def test_18_gamma(self):
        expected = 0.064200
        TestFAO.gamma = psy_const(self.sensor_data['atmosphere_pressure'], TestFAO.latent_heat)
        self.assertAlmostEqual(TestFAO.gamma, expected, delta=0.015)

    def test_19_wind_speed_2m(self):
        expected = 2.34
        TestFAO.wind_speed = wind_speed_2m(self.sensor_data['wind_speed'], self.station_position['sensor_height'])
        self.assertAlmostEqual(TestFAO.wind_speed, expected, delta=0.15)

    def test_20_et0(self):
        expected = 0.1584
        fao = fao56_penman_monteith(TestFAO.net_radiation,
                                    self.sensor_data['temperature_mean'],
                                    TestFAO.wind_speed,
                                    TestFAO.latent_heat,
                                    TestFAO.svp,
                                    TestFAO.actual_vapour_pressure,
                                    TestFAO.delta,
                                    TestFAO.gamma,
                                    TestFAO.soil_heat_flux
                                    )
        self.assertAlmostEqual(fao, expected, delta=0.05)


if __name__ == '__main__':
    unittest.main()
