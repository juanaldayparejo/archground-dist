#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# archGROUND - Python package for ground-based planetary spectroscopy and radiative transfer with archNEMESIS. 
# archnemesis.py - Set of functions to create the archNEMESIS inputs
#
# Copyright (C) 2026 Juan Alday
#
# This file is part of archGROUND.
#
# archGROUND is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import archnemesis as ans
import archground
import numpy as np
import matplotlib.pyplot as plt
import sys,os

from archground.paths import archnemesis_hitran24, archnemesis_tips

###########################################################################################################################

def create_measurement(wavenumber,radiance,error,lat,lon,emiss_ang,sol_ang,azi_ang,
                        lat_subobs=None,lon_subobs=None,v_doppler=0.0,resolving_power=90000.):
    """
    FUNCTION NAME : create_measurement()

    DESCRIPTION : Create a measurement instance with the given parameters

    INPUTS : 

        wavenumber :: Wavenumber grid (cm-1)
        radiance :: Radiance (erg s-1 cm-2 sr-1 (cm-1)-1)
        error :: Radiance error (erg s-1 cm-2 sr-1 (cm-1)-1)
        lat :: Latitude (degrees)
        lon :: Longitude (degrees)
        emiss_ang :: Emission angle (degrees)
        sol_ang :: Solar zenith angle (degrees)
        azi_ang :: Azimuth angle (degrees)

    OPTIONAL INPUTS:
    
        lat_subobs :: Latitude of the sub-observer point (degrees)
        lon_subobs :: Longitude of the sub-observer point (degrees)
            
    OUTPUTS : 
 
        Measurement :: An instance of the `Measurement` class containing the reference measurement from TEXES

    CALLING SEQUENCE:

        Measurement = create_measurement(wavenumber,radiance,error,lat,lon,emiss_ang,sol_ang,azi_ang)

    MODIFICATION HISTORY : Juan Alday (23/08/2026)
    """

    Measurement = ans.Measurement_0()
    Measurement.NGEOM = 1
    Measurement.NCONV = np.zeros(Measurement.NGEOM,dtype='int') + len(wavenumber)
    Measurement.edit_VCONV(wavenumber[:,None])
    Measurement.edit_MEAS(radiance[:,None] * 1.0e-7) #Converting from erg s-1 cm-2 sr-1 (cm-1)-1 to W cm-2 sr-1 (cm-1)-1
    Measurement.edit_ERRMEAS(error[:,None] * 1.0e-7) #Converting from erg s-1 cm-2 sr-1 (cm-1)-1 to W cm-2 sr-1 (cm-1)-1

    Measurement.NAV = np.ones(Measurement.NGEOM,dtype='int')
    Measurement.edit_FLAT(np.ones((Measurement.NGEOM,Measurement.NAV.max())) * lat )
    Measurement.edit_FLON(np.ones((Measurement.NGEOM,Measurement.NAV.max())) * lon )
    Measurement.edit_EMISS_ANG(np.ones((Measurement.NGEOM,Measurement.NAV.max())) * emiss_ang )
    Measurement.edit_SOL_ANG(np.ones((Measurement.NGEOM,Measurement.NAV.max())) * sol_ang )
    Measurement.edit_AZI_ANG(np.ones((Measurement.NGEOM,Measurement.NAV.max())) * azi_ang )
    Measurement.edit_WGEOM(np.ones((Measurement.NGEOM,Measurement.NAV.max())) )

    Measurement.V_DOPPLER = v_doppler

    if lat_subobs is not None:
        Measurement.SUBOBS_LAT = lat_subobs
    else:
        Measurement.SUBOBS_LAT = lat

    if lon_subobs is not None:
        Measurement.SUBOBS_LON = lon_subobs
    else:
        Measurement.SUBOBS_LON = lon

    #Building the ILS
    midwave = np.mean(Measurement.VCONV[:,0])
    fwhm = midwave / resolving_power
    Measurement.ISHAPE = 2 #Gaussian ILS
    Measurement.FWHM = fwhm
    Measurement.build_ils()

    Measurement.assess()

    return Measurement
