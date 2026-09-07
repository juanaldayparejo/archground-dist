#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# archGROUND - Python package for ground-based planetary spectroscopy and radiative transfer with archNEMESIS. 
# venus.py - Set of functions to compute archnemesis models for Venus
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

import numpy as np
from struct import *
import sys,os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import archground
import archnemesis as ans

vcddir = archground.paths.vcddir

###########################################################################################################################

def create_atmosphere_vcd(h,lat,lon,LST):
    """
    FUNCTION NAME : create_atmosphere_vcd()

    DESCRIPTION : Create the Atmosphere class from the VCD and using the Haus et al. (2016) model

    INPUTS : 

        h(nh) :: Altitude (km)
        lat :: Latitude (degrees)
        lon :: Longitude (degrees)
        LST :: Local time (h)

    OPTIONAL INPUTS:
    
        None
            
    OUTPUTS : 
 
        Atmosphere :: An instance of the `Atmosphere_0` class containing the reference atmosphere for Venus

    CALLING SEQUENCE:

        Atmosphere = create_atmosphere_vcd(lat,lon,LST)

    MODIFICATION HISTORY : Juan Alday (23/08/2026)
    """

    #Extracting p-T profile
    press,temp = archground.gcm.extract_vcd.get_pt_profile(h,lat,lon,LST,vcddir=vcddir,euv_scenario=1,albedo_scenario=1)

    #Extracting gaseous abundances
    gas_id, vmr = archground.gcm.extract_vcd.get_vmr_profile(h,lat,lon,LST,vcddir=vcddir,euv_scenario=1,albedo_scenario=1)
    
    #Building Atmosphere
    Atmosphere = ans.Atmosphere_0()
    Atmosphere.NLOCATIONS = 1
    Atmosphere.LATITUDE = lat
    Atmosphere.LONGITUDE = lon
    Atmosphere.IPLANET = 2
    Atmosphere.NP = len(h)
    Atmosphere.edit_H(h*1.0e3)
    Atmosphere.edit_P(press)
    Atmosphere.edit_T(temp)
    Atmosphere.NVMR = len(gas_id)
    Atmosphere.ID = np.array(gas_id,dtype='int')
    Atmosphere.ISO = np.zeros(gas_id.shape,dtype='int')
    Atmosphere.edit_VMR(vmr)
    Atmosphere.assess()

    #Adding the cloud
    n1,n2,n2p,n3 = archground.gcm.extract_vcd.get_haus_cloud_model(h,latitude=Atmosphere.LATITUDE,offset=0.)

    Atmosphere.NDUST = 4
    dust = np.zeros((Atmosphere.NP,Atmosphere.NDUST))
    dust[:,0] = n1 * 1.0e6  #cm-3 to m-3
    dust[:,1] = n2 * 1.0e6  #cm-3 to m-3
    dust[:,2] = n2p * 1.0e6  #cm-3 to m-3
    dust[:,3] = n3 * 1.0e6  #cm-3 to m-3
    Atmosphere.edit_DUST(dust)

    return Atmosphere

###########################################################################################################################

def create_surface_vcd(lat,lon,LST):
    """
    FUNCTION NAME : create_surface_vcd()

    DESCRIPTION : Create the Surface class from the VCD

    INPUTS : 

        lat :: Latitude (degrees)
        lon :: Longitude (degrees)
        LST :: Local time (h)

    OPTIONAL INPUTS:
    
        None
            
    OUTPUTS : 
 
        Surface :: An instance of the `Surface` class containing the reference surface for Venus

    CALLING SEQUENCE:

        Surface = create_surface_vcd(lat,lon,LST)

    MODIFICATION HISTORY : Juan Alday (23/08/2026)
    """

    #Extracting surface temperature
    tsurf = archground.gcm.extract_vcd.get_tsurf(lat,lon,LST,vcddir=vcddir,euv_scenario=1,albedo_scenario=1)

    #Creating surface class
    Surface = ans.Surface_0()
    Surface.NLOCATIONS = 1
    Surface.LATITUDE = lat
    Surface.LONGITUDE = lon
    Surface.TSURF = tsurf
    Surface.LOWBC = 0  #Thermal emission only, no reflection
    Surface.ISPACE = 0 #Wavenumber in cm-1
    Surface.NEM = 2
    Surface.VEM = np.array([0.,100000.]) #Wavenumber grid for the surface emissivity
    Surface.EMISSIVITY = np.ones(Surface.NEM)
    Surface.assess()

    return Surface
