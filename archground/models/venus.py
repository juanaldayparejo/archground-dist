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

def create_atmosphere_vira(h,lat,lon,LST):
    """
    FUNCTION NAME : create_atmosphere_vira()

    DESCRIPTION : Create the Atmosphere class from the VIRA and using the Haus et al. (2016) model

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

        Atmosphere = create_atmosphere_vira(lat,lon,LST)

    MODIFICATION HISTORY : Juan Alday (23/08/2026)
    """

    #Extracting p-T profile from VIRA
    #press,temp = archground.gcm.extract_vcd.get_pt_profile(h,lat,lon,LST,vcddir=vcddir,euv_scenario=1,albedo_scenario=1)
    press, temp, molwt = archground.gcm.extract_vcd.get_vira_profiles(h, latitude=lat)

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
    Atmosphere.MOLWT = molwt / 1.0e3  #g/mol to kg/mol
    Atmosphere.AMFORM = 0
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

###########################################################################################################################

def create_scatter_class(waven,iscat=0,nmu=5,nf=10,nphi=100):
    """
    FUNCTION NAME : create_cloud_model()

    DESCRIPTION : Create an instance of the `Scatter_0` class containing the cloud model for the Venus case

    INPUTS : 

        waven :: Wavenumber array of the measurement in cm-1

    OPTIONAL INPUTS:
    
        iscat :: Flag indicating the type of scattering calculations
                    iscat = 0 : No scattering
                    iscat = 1 : Multiple scattering
            
    OUTPUTS : 
 
        Scatter :: An instance of the `Scatter_0` class containing the cloud model for the Venus case

    CALLING SEQUENCE:

        Scatter = create_cloud_model(Atmosphere)

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    waven_min = waven.min() ; waven_max = waven.max()

    #Defining cloud modes
    #Mode 1
    r_g1 = 0.3
    sigma_g1 = 0.445

    #Mode 2
    r_g2 = 1.0
    sigma_g2 = 0.255

    #Mode 2'
    r_g2p = 1.40
    sigma_g2p = 0.207

    #Mode 3
    r_g3 = 3.65
    sigma_g3 = 0.245

    #First of all we need to define our wavelength or wavenumber array
    #and tell the class in what units we want these calculations (wavenumber in cm-1 (ISPACE=0) or wavelength in um (ISPACE=1))
    Scatter = ans.Scatter_0()
    Scatter.read_refind(4)  #Reading optical properties of the main condensates in the Venus atmosphere (sulfuric acid aerosols)
    Scatter.ISCAT = iscat
    Scatter.NF = nf
    Scatter.NMU = nmu
    Scatter.NPHI = nphi
    Scatter.ISPACE = 0 #Wavelength

    Scatter.IMIE = 2 #Legendre polynomial expansion of the phase function
    waven = np.arange(int(waven_min)-5.,waven_max+6.,1.)
    NDUST = 4      #Number of aerosol populations that we want to include in our atmosphere
    NWAVE = len(waven)    #Number of spectral points
    NTHETA = 361    #Number of phase angles for defining the phase function
    theta = np.linspace(0.,180.,NTHETA)

    #Now we initialise the arrays that will be filled with the calculations
    Scatter.initialise_arrays(NDUST,NWAVE,NTHETA,NLPOL=150)
    Scatter.WAVE = waven
    Scatter.THETA = theta

    #Mode 1
    print('Calculating Mode 1')
    iscat = 2  #Log-normal distribution
    pars = np.array([r_g1,sigma_g1])

    idust = 0    #The index of the aerosol populations in the class that this calculation corresponds to (from 0 to NDUST-1)
    Scatter.makephase(idust,iscat,pars)

    #Mode 2
    print('Calculating Mode 2')
    iscat = 2  #Log-normal distribution
    pars = np.array([r_g2,sigma_g2])

    idust = 1    #The index of the aerosol populations in the class that this calculation corresponds to (from 0 to NDUST-1)
    Scatter.makephase(idust,iscat,pars)

    #Mode 2'
    print('Calculating Mode 2p')
    iscat = 2  #Log-normal distribution
    pars = np.array([r_g2p,sigma_g2p])

    idust = 2    #The index of the aerosol populations in the class that this calculation corresponds to (from 0 to NDUST-1)
    Scatter.makephase(idust,iscat,pars)

    #Mode 3
    print('Calculating Mode 3')
    iscat = 2  #Log-normal distribution
    pars = np.array([r_g3,sigma_g3])

    idust = 3    #The index of the aerosol populations in the class that this calculation corresponds to (from 0 to NDUST-1)
    Scatter.makephase(idust,iscat,pars)

    return Scatter
