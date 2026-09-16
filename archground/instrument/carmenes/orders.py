import archground


def find_order(wavenumber, channel):
    """Find physical diffraction orders for a wavenumber in cm^-1.

    channel: 'VIS' or 'NIR' (case-insensitive).
    Returns an int for one match, a tuple of ints for overlapping orders,
    or float('nan') for no match (including nonpositive/nonfinite input).
    Bounds are inclusive and use the tabulated wavenumber limits.
    Invalid channel names raise ValueError.

    NIR internal detector gaps cannot be identified from these outer bounds.
    A match therefore means within an order's span, not necessarily a
    wavelength recorded by a detector.
    """
    import math

    channel = str(channel).strip().upper()
    if channel not in CARMENES_ORDERS:
        raise ValueError("channel must be 'VIS' or 'NIR'")

    wavenumber = float(wavenumber)
    if not math.isfinite(wavenumber) or wavenumber <= 0:
        return float('nan')

    matches = tuple(
        order for order, limits in CARMENES_ORDERS[channel].items()
        if limits['wavenumber_min_cm1'] <= wavenumber
        <= limits['wavenumber_max_cm1']
    )
    if not matches:
        return float('nan')
    return matches[0] if len(matches) == 1 else matches


######################################################################################################3

CARMENES_ORDERS = {
    "VIS": {
        118: {
            "wavelength_min_A": 5171.4,
            "wavelength_max_A": 5206.7,
            "wavenumber_min_cm1": 19206.023009,
            "wavenumber_max_cm1": 19337.12341,
        },
        117: {
            "wavelength_min_A": 5215.4,
            "wavelength_max_A": 5254.8,
            "wavenumber_min_cm1": 19030.219989,
            "wavenumber_max_cm1": 19173.984738,
        },
        116: {
            "wavelength_min_A": 5255.3,
            "wavelength_max_A": 5306.6,
            "wavenumber_min_cm1": 18844.457845,
            "wavenumber_max_cm1": 19028.409415,
        },
        115: {
            "wavelength_min_A": 5297.6,
            "wavelength_max_A": 5355.6,
            "wavenumber_min_cm1": 18672.044215,
            "wavenumber_max_cm1": 18876.472365,
        },
        114: {
            "wavelength_min_A": 5339.0,
            "wavelength_max_A": 5404.7,
            "wavenumber_min_cm1": 18502.414565,
            "wavenumber_max_cm1": 18730.09927,
        },
        113: {
            "wavelength_min_A": 5383.6,
            "wavelength_max_A": 5455.3,
            "wavenumber_min_cm1": 18330.797573,
            "wavenumber_max_cm1": 18574.931273,
        },
        112: {
            "wavelength_min_A": 5429.7,
            "wavelength_max_A": 5504.9,
            "wavenumber_min_cm1": 18165.634253,
            "wavenumber_max_cm1": 18417.223788,
        },
        111: {
            "wavelength_min_A": 5477.0,
            "wavelength_max_A": 5555.8,
            "wavenumber_min_cm1": 17999.208035,
            "wavenumber_max_cm1": 18258.170531,
        },
        110: {
            "wavelength_min_A": 5524.5,
            "wavelength_max_A": 5606.3,
            "wavenumber_min_cm1": 17837.076146,
            "wavenumber_max_cm1": 18101.185628,
        },
        109: {
            "wavelength_min_A": 5573.3,
            "wavelength_max_A": 5657.8,
            "wavenumber_min_cm1": 17674.714553,
            "wavenumber_max_cm1": 17942.691045,
        },
        108: {
            "wavelength_min_A": 5622.9,
            "wavelength_max_A": 5710.1,
            "wavenumber_min_cm1": 17512.828147,
            "wavenumber_max_cm1": 17784.417294,
        },
        107: {
            "wavelength_min_A": 5673.5,
            "wavelength_max_A": 5763.5,
            "wavenumber_min_cm1": 17350.568231,
            "wavenumber_max_cm1": 17625.804177,
        },
        106: {
            "wavelength_min_A": 5725.7,
            "wavelength_max_A": 5817.9,
            "wavenumber_min_cm1": 17188.33256,
            "wavenumber_max_cm1": 17465.113436,
        },
        105: {
            "wavelength_min_A": 5778.6,
            "wavelength_max_A": 5873.2,
            "wavenumber_min_cm1": 17026.493223,
            "wavenumber_max_cm1": 17305.22964,
        },
        104: {
            "wavelength_min_A": 5832.8,
            "wavelength_max_A": 5929.7,
            "wavenumber_min_cm1": 16864.259575,
            "wavenumber_max_cm1": 17144.424633,
        },
        103: {
            "wavelength_min_A": 5887.8,
            "wavelength_max_A": 5987.3,
            "wavenumber_min_cm1": 16702.019274,
            "wavenumber_max_cm1": 16984.272564,
        },
        102: {
            "wavelength_min_A": 5943.6,
            "wavelength_max_A": 6046.0,
            "wavenumber_min_cm1": 16539.861065,
            "wavenumber_max_cm1": 16824.819974,
        },
        101: {
            "wavelength_min_A": 6001.2,
            "wavelength_max_A": 6105.8,
            "wavenumber_min_cm1": 16377.870222,
            "wavenumber_max_cm1": 16663.334,
        },
        100: {
            "wavelength_min_A": 6061.3,
            "wavelength_max_A": 6166.8,
            "wavenumber_min_cm1": 16215.865603,
            "wavenumber_max_cm1": 16498.110966,
        },
        99: {
            "wavelength_min_A": 6122.5,
            "wavelength_max_A": 6229.1,
            "wavenumber_min_cm1": 16053.683518,
            "wavenumber_max_cm1": 16333.197223,
        },
        98: {
            "wavelength_min_A": 6185.0,
            "wavelength_max_A": 6292.7,
            "wavenumber_min_cm1": 15891.429752,
            "wavenumber_max_cm1": 16168.148747,
        },
        97: {
            "wavelength_min_A": 6248.7,
            "wavelength_max_A": 6357.5,
            "wavenumber_min_cm1": 15729.453401,
            "wavenumber_max_cm1": 16003.328692,
        },
        96: {
            "wavelength_min_A": 6313.8,
            "wavelength_max_A": 6423.7,
            "wavenumber_min_cm1": 15567.352149,
            "wavenumber_max_cm1": 15838.322405,
        },
        95: {
            "wavelength_min_A": 6380.3,
            "wavelength_max_A": 6491.3,
            "wavenumber_min_cm1": 15405.234699,
            "wavenumber_max_cm1": 15673.244205,
        },
        94: {
            "wavelength_min_A": 6448.2,
            "wavelength_max_A": 6560.4,
            "wavenumber_min_cm1": 15242.972989,
            "wavenumber_max_cm1": 15508.20384,
        },
        93: {
            "wavelength_min_A": 6517.5,
            "wavelength_max_A": 6630.9,
            "wavenumber_min_cm1": 15080.909077,
            "wavenumber_max_cm1": 15343.306483,
        },
        92: {
            "wavelength_min_A": 6588.4,
            "wavelength_max_A": 6703.0,
            "wavenumber_min_cm1": 14918.693122,
            "wavenumber_max_cm1": 15178.191974,
        },
        91: {
            "wavelength_min_A": 6660.8,
            "wavelength_max_A": 6776.6,
            "wavenumber_min_cm1": 14756.662633,
            "wavenumber_max_cm1": 15013.211626,
        },
        90: {
            "wavelength_min_A": 6734.8,
            "wavelength_max_A": 6851.9,
            "wavenumber_min_cm1": 14594.492039,
            "wavenumber_max_cm1": 14848.250876,
        },
        89: {
            "wavelength_min_A": 6810.5,
            "wavelength_max_A": 6928.9,
            "wavenumber_min_cm1": 14432.305272,
            "wavenumber_max_cm1": 14683.20975,
        },
        88: {
            "wavelength_min_A": 6887.9,
            "wavelength_max_A": 7007.6,
            "wavenumber_min_cm1": 14270.220903,
            "wavenumber_max_cm1": 14518.213098,
        },
        87: {
            "wavelength_min_A": 6967.0,
            "wavelength_max_A": 7088.1,
            "wavenumber_min_cm1": 14108.153102,
            "wavenumber_max_cm1": 14353.380221,
        },
        86: {
            "wavelength_min_A": 7048.1,
            "wavelength_max_A": 7170.5,
            "wavenumber_min_cm1": 13946.028868,
            "wavenumber_max_cm1": 14188.220939,
        },
        85: {
            "wavelength_min_A": 7131.0,
            "wavelength_max_A": 7254.9,
            "wavenumber_min_cm1": 13783.787509,
            "wavenumber_max_cm1": 14023.278643,
        },
        84: {
            "wavelength_min_A": 7215.9,
            "wavelength_max_A": 7341.2,
            "wavenumber_min_cm1": 13621.751212,
            "wavenumber_max_cm1": 13858.285176,
        },
        83: {
            "wavelength_min_A": 7302.8,
            "wavelength_max_A": 7429.7,
            "wavenumber_min_cm1": 13459.493654,
            "wavenumber_max_cm1": 13693.377882,
        },
        82: {
            "wavelength_min_A": 7391.9,
            "wavelength_max_A": 7520.3,
            "wavenumber_min_cm1": 13297.341861,
            "wavenumber_max_cm1": 13528.321541,
        },
        81: {
            "wavelength_min_A": 7483.1,
            "wavelength_max_A": 7613.1,
            "wavenumber_min_cm1": 13135.253707,
            "wavenumber_max_cm1": 13363.445631,
        },
        80: {
            "wavelength_min_A": 7576.7,
            "wavelength_max_A": 7708.2,
            "wavenumber_min_cm1": 12973.197374,
            "wavenumber_max_cm1": 13198.358124,
        },
        79: {
            "wavelength_min_A": 7672.6,
            "wavelength_max_A": 7805.8,
            "wavenumber_min_cm1": 12810.986702,
            "wavenumber_max_cm1": 13033.391549,
        },
        78: {
            "wavelength_min_A": 7771.0,
            "wavelength_max_A": 7905.8,
            "wavenumber_min_cm1": 12648.941284,
            "wavenumber_max_cm1": 12868.356711,
        },
        77: {
            "wavelength_min_A": 7871.9,
            "wavelength_max_A": 8008.5,
            "wavenumber_min_cm1": 12486.732846,
            "wavenumber_max_cm1": 12703.413407,
        },
        76: {
            "wavelength_min_A": 7975.5,
            "wavelength_max_A": 8113.8,
            "wavenumber_min_cm1": 12324.681407,
            "wavenumber_max_cm1": 12538.398846,
        },
        75: {
            "wavelength_min_A": 8081.8,
            "wavelength_max_A": 8222.0,
            "wavenumber_min_cm1": 12162.490878,
            "wavenumber_max_cm1": 12373.481155,
        },
        74: {
            "wavelength_min_A": 8191.0,
            "wavelength_max_A": 8333.1,
            "wavenumber_min_cm1": 12000.336009,
            "wavenumber_max_cm1": 12208.521548,
        },
        73: {
            "wavelength_min_A": 8303.2,
            "wavelength_max_A": 8447.2,
            "wavenumber_min_cm1": 11838.242258,
            "wavenumber_max_cm1": 12043.549475,
        },
        72: {
            "wavelength_min_A": 8418.6,
            "wavelength_max_A": 8564.5,
            "wavenumber_min_cm1": 11676.104851,
            "wavenumber_max_cm1": 11878.459601,
        },
        71: {
            "wavelength_min_A": 8537.2,
            "wavelength_max_A": 8685.1,
            "wavenumber_min_cm1": 11513.972205,
            "wavenumber_max_cm1": 11713.442346,
        },
        70: {
            "wavelength_min_A": 8659.1,
            "wavelength_max_A": 8809.2,
            "wavenumber_min_cm1": 11351.768606,
            "wavenumber_max_cm1": 11548.544306,
        },
        69: {
            "wavelength_min_A": 8784.6,
            "wavelength_max_A": 8936.8,
            "wavenumber_min_cm1": 11189.687584,
            "wavenumber_max_cm1": 11383.557589,
        },
        68: {
            "wavelength_min_A": 8913.8,
            "wavelength_max_A": 9068.2,
            "wavenumber_min_cm1": 11027.546812,
            "wavenumber_max_cm1": 11218.559986,
        },
        67: {
            "wavelength_min_A": 9046.9,
            "wavelength_max_A": 9203.5,
            "wavenumber_min_cm1": 10865.431629,
            "wavenumber_max_cm1": 11053.510042,
        },
        66: {
            "wavelength_min_A": 9183.9,
            "wavelength_max_A": 9342.9,
            "wavenumber_min_cm1": 10703.314817,
            "wavenumber_max_cm1": 10888.620303,
        },
        65: {
            "wavelength_min_A": 9325.2,
            "wavelength_max_A": 9486.6,
            "wavenumber_min_cm1": 10541.184407,
            "wavenumber_max_cm1": 10723.630592,
        },
        64: {
            "wavelength_min_A": 9471.0,
            "wavelength_max_A": 9634.8,
            "wavenumber_min_cm1": 10379.042637,
            "wavenumber_max_cm1": 10558.547144,
        },
    },
    "NIR": {
        63: {
            "wavelength_min_A": 9618.7,
            "wavelength_max_A": 9782.1,
            "wavenumber_min_cm1": 10222.753805,
            "wavenumber_max_cm1": 10396.415316,
        },
        62: {
            "wavelength_min_A": 9758.4,
            "wavelength_max_A": 9939.9,
            "wavenumber_min_cm1": 10060.463385,
            "wavenumber_max_cm1": 10247.581571,
        },
        61: {
            "wavelength_min_A": 9918.3,
            "wavelength_max_A": 10102.9,
            "wavenumber_min_cm1": 9898.148056,
            "wavenumber_max_cm1": 10082.372987,
        },
        60: {
            "wavelength_min_A": 10083.6,
            "wavelength_max_A": 10271.2,
            "wavenumber_min_cm1": 9735.960745,
            "wavenumber_max_cm1": 9917.093102,
        },
        59: {
            "wavelength_min_A": 10254.5,
            "wavelength_max_A": 10445.3,
            "wavenumber_min_cm1": 9573.683858,
            "wavenumber_max_cm1": 9751.816276,
        },
        58: {
            "wavelength_min_A": 10431.3,
            "wavelength_max_A": 10625.4,
            "wavenumber_min_cm1": 9411.410394,
            "wavenumber_max_cm1": 9586.532839,
        },
        57: {
            "wavelength_min_A": 10614.3,
            "wavelength_max_A": 10811.8,
            "wavenumber_min_cm1": 9249.153702,
            "wavenumber_max_cm1": 9421.252461,
        },
        56: {
            "wavelength_min_A": 10803.9,
            "wavelength_max_A": 11004.9,
            "wavenumber_min_cm1": 9086.861307,
            "wavenumber_max_cm1": 9255.916845,
        },
        55: {
            "wavelength_min_A": 11000.3,
            "wavelength_max_A": 11205.0,
            "wavenumber_min_cm1": 8924.587238,
            "wavenumber_max_cm1": 9090.661164,
        },
        54: {
            "wavelength_min_A": 11204.0,
            "wavelength_max_A": 11412.5,
            "wavenumber_min_cm1": 8762.322015,
            "wavenumber_max_cm1": 8925.383792,
        },
        53: {
            "wavelength_min_A": 11415.4,
            "wavelength_max_A": 11627.8,
            "wavenumber_min_cm1": 8600.079121,
            "wavenumber_max_cm1": 8760.096011,
        },
        52: {
            "wavelength_min_A": 11634.9,
            "wavelength_max_A": 11851.4,
            "wavenumber_min_cm1": 8437.821692,
            "wavenumber_max_cm1": 8594.831069,
        },
        51: {
            "wavelength_min_A": 11863.1,
            "wavelength_max_A": 12083.7,
            "wavenumber_min_cm1": 8275.610947,
            "wavenumber_max_cm1": 8429.499878,
        },
        50: {
            "wavelength_min_A": 12100.4,
            "wavelength_max_A": 12325.4,
            "wavenumber_min_cm1": 8113.326951,
            "wavenumber_max_cm1": 8264.189614,
        },
        49: {
            "wavelength_min_A": 12347.3,
            "wavelength_max_A": 12576.9,
            "wavenumber_min_cm1": 7951.084926,
            "wavenumber_max_cm1": 8098.93661,
        },
        48: {
            "wavelength_min_A": 12604.5,
            "wavelength_max_A": 12839.0,
            "wavenumber_min_cm1": 7788.768596,
            "wavenumber_max_cm1": 7933.674481,
        },
        47: {
            "wavelength_min_A": 12872.7,
            "wavelength_max_A": 13112.1,
            "wavenumber_min_cm1": 7626.543422,
            "wavenumber_max_cm1": 7768.37804,
        },
        46: {
            "wavelength_min_A": 13152.6,
            "wavelength_max_A": 13397.1,
            "wavenumber_min_cm1": 7464.301976,
            "wavenumber_max_cm1": 7603.059471,
        },
        45: {
            "wavelength_min_A": 13444.9,
            "wavelength_max_A": 13694.8,
            "wavenumber_min_cm1": 7302.041651,
            "wavenumber_max_cm1": 7437.764506,
        },
        44: {
            "wavelength_min_A": 13750.5,
            "wavelength_max_A": 14006.1,
            "wavenumber_min_cm1": 7139.746253,
            "wavenumber_max_cm1": 7272.46282,
        },
        43: {
            "wavelength_min_A": 14070.3,
            "wavelength_max_A": 14331.8,
            "wavenumber_min_cm1": 6977.490615,
            "wavenumber_max_cm1": 7107.169001,
        },
        42: {
            "wavelength_min_A": 14405.3,
            "wavelength_max_A": 14673.0,
            "wavenumber_min_cm1": 6815.238874,
            "wavenumber_max_cm1": 6941.889443,
        },
        41: {
            "wavelength_min_A": 14756.6,
            "wavelength_max_A": 15030.8,
            "wavenumber_min_cm1": 6653.005828,
            "wavenumber_max_cm1": 6776.628763,
        },
        40: {
            "wavelength_min_A": 15125.6,
            "wavelength_max_A": 15406.6,
            "wavenumber_min_cm1": 6490.724754,
            "wavenumber_max_cm1": 6611.307981,
        },
        39: {
            "wavelength_min_A": 15513.4,
            "wavelength_max_A": 15801.6,
            "wavenumber_min_cm1": 6328.473066,
            "wavenumber_max_cm1": 6446.040198,
        },
        38: {
            "wavelength_min_A": 15921.7,
            "wavelength_max_A": 16217.4,
            "wavenumber_min_cm1": 6166.216533,
            "wavenumber_max_cm1": 6280.736354,
        },
        37: {
            "wavelength_min_A": 16352.1,
            "wavelength_max_A": 16655.6,
            "wavenumber_min_cm1": 6003.986647,
            "wavenumber_max_cm1": 6115.422484,
        },
        36: {
            "wavelength_min_A": 16806.3,
            "wavelength_max_A": 17118.2,
            "wavenumber_min_cm1": 5841.735697,
            "wavenumber_max_cm1": 5950.149646,
        },
    },
}