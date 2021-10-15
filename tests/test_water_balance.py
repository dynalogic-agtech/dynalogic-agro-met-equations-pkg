from AgroMetEquations.thornthwaite_mather_1955_waterbalance import *


def main():
    cad = 78

    input_data = [
        {"et0": 42, "rain_fall": 127, "kc": 1},
        {"et0": 41, "rain_fall": 158, "kc": 1},
        {"et0": 44, "rain_fall": 189, "kc": 1},
        {"et0": 39, "rain_fall": 60, "kc": 1},
        {"et0": 38, "rain_fall": 41, "kc": 1},
        {"et0": 29, "rain_fall": 40, "kc": 1},
        {"et0": 36, "rain_fall": 133, "kc": 0.3},
        {"et0": 34, "rain_fall": 102, "kc": 0.4},
        {"et0": 34, "rain_fall": 71, "kc": 0.5},
        {"et0": 30, "rain_fall": 25, "kc": 0.6},
        {"et0": 28, "rain_fall": 14, "kc": 0.7},
        {"et0": 26, "rain_fall": 17, "kc": 0.9},
        {"et0": 24, "rain_fall": 7, "kc": 1},
        {"et0": 22, "rain_fall": 2, "kc": 1.2},
        {"et0": 23, "rain_fall": 8, "kc": 1.2},
        {"et0": 19, "rain_fall": 0, "kc": 1},
        {"et0": 17, "rain_fall": 0, "kc": 0.9},
        {"et0": 17, "rain_fall": 0, "kc": 0.8},
        {"et0": 18, "rain_fall": 0, "kc": 0.5},
        {"et0": 19, "rain_fall": 0, "kc": 1},
        {"et0": 23, "rain_fall": 0, "kc": 1},
        {"et0": 23, "rain_fall": 6, "kc": 1},
        {"et0": 24, "rain_fall": 10, "kc": 1},
        {"et0": 29, "rain_fall": 13, "kc": 1},
        {"et0": 28, "rain_fall": 0, "kc": 1},
        {"et0": 30, "rain_fall": 9, "kc": 1},
        {"et0": 32, "rain_fall": 8, "kc": 1},
        {"et0": 34, "rain_fall": 70, "kc": 0.3},
        {"et0": 36, "rain_fall": 25, "kc": 0.4},
        {"et0": 40, "rain_fall": 31, "kc": 0.5},
        {"et0": 38, "rain_fall": 120, "kc": 0.6},
        {"et0": 39, "rain_fall": 86, "kc": 0.7},
        {"et0": 40, "rain_fall": 38, "kc": 0.9},
        {"et0": 40, "rain_fall": 53, "kc": 1},
        {"et0": 41, "rain_fall": 70, "kc": 1.2},
        {"et0": 45, "rain_fall": 87, "kc": 1.2},
        {"et0": 42, "rain_fall": 84, "kc": 1},
        {"et0": 41, "rain_fall": 132, "kc": 0.9},
        {"et0": 44, "rain_fall": 115, "kc": 0.8},
        {"et0": 39, "rain_fall": 11, "kc": 0.5}
    ]

    cnt = 0
    result_list = []

    for data in input_data:
        if cnt == 0:
            last_negative_accumulated = 0
            last_water_storage = cad

        cnt = cnt + 1

        water_balance_info = waterbalance(
            data["et0"],
            data["kc"],
            data["rain_fall"],
            last_negative_accumulated,
            last_water_storage,
            cad)

        last_negative_accumulated = water_balance_info["negative_accumulated"]
        last_water_storage = water_balance_info["water_storage"]

        result_obj = {
            "id": cnt,
            "relative_water_storage": water_balance_info["relative_water_storage"],
            "water_storage": water_balance_info["water_storage"],
            "negative_accumulated": water_balance_info["negative_accumulated"],
            "water_deficit": water_balance_info["water_deficit"],
            "water_excess": water_balance_info["water_excess"],
            "alt": water_balance_info["alt"],
            "etr": water_balance_info["etr"],
            "etc": water_balance_info["etc"]
        }
        result_list.append(result_obj)

    return result_list


def test_day_1():
    result = main()
    assert result[0]["relative_water_storage"] == 100


if __name__ == "__main__":
    result = main()

    print(result)
