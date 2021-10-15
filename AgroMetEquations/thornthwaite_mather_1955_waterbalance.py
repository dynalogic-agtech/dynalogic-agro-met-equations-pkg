import math


def waterbalance(et0: float,
                 kc: float,
                 rain_fall: float,
                 last_negative_accumulated: int,
                 last_water_storage: int,
                 cad: int) -> dict:
    """
    Thornthwaite & Mather (1955)

    :param et0: float -> Reference Evapotranspiration
    :param kc: float -> Crop Coefficient
    :param rain_fall: float -> Sum of rail fall
    :param last_negative_accumulated: int -> Negative accumulated of water
    :param last_water_storage: int -> Last reported water storage
    :param cad:
    :return: Current water_storage: int,
    """
    # Step 1 -> Calculate etc
    etc = round(et0 * kc)

    # Step 2 - > Calculate p_etc
    p_etc = round(rain_fall - etc)

    # Step 3 -> Calculate water_storage and negative_accumulated
    water_storage, negative_accumulated = get_water_storage_negative_accumulated(last_water_storage,
                                                                                 cad,
                                                                                 last_negative_accumulated,
                                                                                 p_etc)

    # Step 4 -> Calculate difference between current water_storage and last water_storage
    alt = water_storage - last_water_storage

    # Step 5 -> Calculate etr
    etr = etc if p_etc >= 0 else round(rain_fall + abs(alt))

    # Step 6 -> Calculate water deficit
    water_deficit = etc - etr

    # Step 7 -> Calculate the water excess
    water_excess = p_etc - alt if water_storage == cad else 0

    # Step 8 -> Calculate Relative Water Storage
    relative_water_storage = round((water_storage / cad) * 100)

    return {
        "relative_water_storage": relative_water_storage,
        "water_storage": water_storage,
        "negative_accumulated": negative_accumulated,
        "water_deficit": water_deficit,
        "water_excess": water_excess,
        "alt": alt,
        "etr": etr,
        "etc": etc
    }


def get_water_storage_negative_accumulated(last_water_storage: int,
                                           cad: int,
                                           last_negative_accumulated: int,
                                           p_etc: int) -> tuple:

    water_storage = None
    negative_accumulated = None

    # First negative value for p_etc
    if p_etc < 0 and last_negative_accumulated == 0:
        negative_accumulated = p_etc
        water_storage = get_water_storage(p_etc, negative_accumulated, cad, last_water_storage)

    # p_etc negative
    if p_etc < 0 and last_negative_accumulated < 0:
        negative_accumulated = last_negative_accumulated + p_etc
        water_storage = get_water_storage(p_etc, negative_accumulated, cad, last_water_storage)

    # p_etc positive after negative sequence
    if p_etc >= 0 >= last_negative_accumulated:
        water_storage = get_water_storage(p_etc, last_negative_accumulated, cad, last_water_storage)
        negative_accumulated = round(cad * math.log(water_storage / cad)) if water_storage < cad else 0

    return water_storage, negative_accumulated


def get_water_storage(p_etc: int, negative_accumulated: int, cad: int, last_water_storage: int) -> int:
    if p_etc < 0:
        water_storage = round(cad * math.exp(negative_accumulated / cad))
    else:
        water_storage = round(last_water_storage + p_etc) if (last_water_storage + p_etc) < cad else cad

    return water_storage
