#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# archGROUND - Python package for ground-based planetary spectroscopy and radiative transfer with archNEMESIS. 
# rtm.general.py - Set of functions to run a generic forward model. This is typically useful for proposals.
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

import archground
import archnemesis as ans
import numpy as np
import matplotlib.pyplot as plt

######################################################################################################################################

def run_forward_model(outname,planet,observatory,date,time,
                        waven_min,waven_max,resolving_power,
                        id_gases,iso_gases,vmr_gases,
                        x0=0.,y0=0.,
                        iscat=1,
                        hitran_file=archground.paths.archnemesis_hitran24):
    """
    FUNCTION NAME : run_forward_model_subobserver_point()

    DESCRIPTION : Run a generic forward model at the sub-observer point 

    INPUTS : 

        outname :: Name of the output files
        planet :: Name of the planet 
        observatory :: Name of the observatory (see dictionary in geometry.py)
        date :: Observation date
        time :: Observation time (UTC)

        waven_min :: Minimum wavenumber of the measurement (cm-1)
        waven_max :: Maximum wavenumber of the measurement (cm-1)
        resolving_power :: Resolving power of the instrument

        id_gases(ngas) :: Radtran ID of the radiatively active gases to include
        iso_gases(ngas) :: Radtran isotope ID of the radiatively active gases to include
        vmr_gases(ngas) :: VMR of the gases to include. If None, it will use the reference ones from the GCMs

    OPTIONAL INPUTS:
    
        x0 :: Offset of the observing point with respect to the centre of the planet's disk (arcsec)
        y0 :: Offset of the observing point with respect to the centre of the planet's disk (arcsec)
        iscat :: Multiple scattering flag. 0 = thermal emission ; 1 = multiple scattering
        hitran_file :: Path of the spectroscopy file to use (default = HITRAN24 database)

    OUTPUTS : 
 
        waven(nconv) :: Wavenumber array (cm-1)
        speconv(nconv) :: Modelled spectrum

    CALLING SEQUENCE:

        waven, speconv = run_forward_model_subobserver_point(
                                        planet,observatory,date,time,
                                        wave_min,wave_max,resolving_power,
                                        id_gases,iso_gases,vmr_gases,
                                        )

    MODIFICATION HISTORY : Juan Alday (24/09/2026)
    """

    ############################################################################
    #CALCULATING GEOMETRY
    ############################################################################

    #Calculating the ephemerides
    ephemerides = archground.geometry.get_ephemerides_epoch(date,time,planet,observatory)

    #Calculating the observation maps
    maps = archground.geometry.build_maps_epoch(ephemerides,res_bin=0.01)

    #Getting the geometry at the required x,y
    geometry = archground.geometry.get_geometry_epoch(maps, x0, y0)

    elevation = ephemerides["EL"][0]  #elevation above local horizon (degrees)
    if elevation < 0.0:
        print("elevation above local horizon = ",elevation,"degrees")
        raise ValueError("error :: planet is not observable at that time")

    ############################################################################
    #CREATING REFERENCE CLASSES
    ############################################################################

    #Atmosphere
    ########################################################################

    if planet.upper() == "MARS":

        h = np.linspace(0.,100.,51)
        lat = geometry["lat"]
        lon = geometry["lon"]
        LST = geometry["lst"]
        Ls = ephemerides["App_Lon_Sun"][0]

        Atmosphere = archground.models.mars.create_atmosphere_mcd(h,lat,lon,Ls,LST,include_dust=True)

    elif planet.upper() == "VENUS":

        h = np.linspace(0.,120.,61)
        lat = geometry["lat"]
        lon = geometry["lon"]
        LST = geometry["lst"]
        Ls = ephemerides["App_Lon_Sun"][0]

        Atmosphere = archground.models.venus.create_atmosphere_vcd(h,lat,lon,LST)

    else:

        raise ValueError("error :: selected planet has not been implemented yet")

    
    #Measurement
    ########################################################################

    #Spectral parameters
    fwhm = 0.5*(waven_min + waven_max) / resolving_power #cm-1
    delv = fwhm / 5.
    wavenumber = np.arange(waven_min,waven_max+delv,delv)
    radiance = np.ones_like(wavenumber)
    error = np.ones_like(wavenumber)

    #Geometry
    lat = geometry["lat"]
    lon = geometry["lon"]
    emiss_ang = geometry["emiss_ang"]
    sol_ang = geometry["sol_ang"]
    azi_ang = geometry["azi_ang"]
    v_planet_observer = ephemerides["delta_rate"][0]

    #Calculating sub-observer latitude and longitude (planetocentric)
    equatorial_radius = archground.geometry.planet_dict[planet.upper()]["equatorial_radius_km"]
    polar_radius = archground.geometry.planet_dict[planet.upper()]["polar_radius_km"]
    lat_subobs = ephemerides["PDObsLat"][0]        #Planetodetic sub-observer latitude (degrees)
    lon_subobs = (-ephemerides["PDObsLon"][0] + 180.0) % 360.0 - 180.0        #Planetodetic sub-observer longitude (degrees) - East positive
    lat_subobs = archground.geometry.planetodetic_to_planetocentric(
        lat_subobs,
        equatorial_radius=equatorial_radius,
        polar_radius=polar_radius,
    )

    Measurement = archground.models.general.create_measurement(wavenumber,radiance,error,lat,lon,emiss_ang,sol_ang,azi_ang,
                            lat_subobs=lat_subobs,lon_subobs=lon_subobs,v_doppler=v_planet_observer,resolving_power=resolving_power)


    #############################################################

    if planet.upper() == "MARS":
        Scatter = archground.models.mars.create_scatter_class(wavenumber,iscat=iscat,nmu=5,nf=10,nphi=100,include_waterice=False)
    elif plannet.upper() == "VENUS":
        Scatter = archground.models.venus.create_scatter_class(wavenumber,iscat=iscat,nmu=5,nf=10,nphi=100)
    
    #Surface 
    #############################################################

    if planet.upper() == "MARS":
        Surface = archground.models.mars.create_surface_mcd(lat,lon,Ls,LST)
    elif planet.upper() == "VENUS":
        Surface = archground.models.venus.create_surface_vcd(lat,lon,LST)

    #Solar
    #############################################################

    #Getting relevant parameters from ephemerides
    distance_sun_planet = ephemerides["r"][0]   #AU
    velocity_sun_planet = ephemerides["r_rate"][0]  #km/s

    Stellar = archground.models.general.create_stellar_class(wavenumber,distance_sun_planet,velocity_sun_planet,t_sun=5772.0)

    #Telluric
    #############################################################

    elevation = ephemerides["EL"][0]  #elevation above local horizon (degrees)

    Telluric = archground.models.general.create_telluric_online(observatory,wavenumber,date,time,elevation,
                                                                era5=False,
                                                                resolving_power=resolving_power,id_gas_act=[1,2,3,4,5,6,7])

    #Spectroscopy
    #############################################################

    #Checking if isotopes must be split
    ngas = len(id_gases)
    flag_split_H2O_isotopes = False
    flag_split_CO2_isotopes = False
    flag_split_CO_isotopes = False
    for igas in range(ngas):
        if id_gases[igas] == 1:
            if iso_gases[igas] > 0:
                flag_split_H2O_isotopes = True
        if id_gases[igas] == 2:
            if iso_gases[igas] > 0:
                flag_split_CO2_isotopes = True
        if id_gases[igas] == 5:
            if iso_gases[igas] > 0:
                flag_split_CO_isotopes = True

    if planet.upper() == "MARS":
        dhratio=5.
        o18ratio=1.
        o17ratio=1.
        c13ratio=1.
    elif planet.upper() == "VENUS":
        dhratio=120.
        o18ratio=1.
        o17ratio=1.
        c13ratio=1.

    if flag_split_H2O_isotopes:
        Atmosphere = archground.models.general.split_H2O_isotopes(Atmosphere,dhratio=dhratio,o18ratio=o18ratio,o17ratio=o17ratio)
    if flag_split_CO2_isotopes:
        Atmosphere = archground.models.general.split_CO2_isotopes(Atmosphere,c13ratio=c13ratio,o18ratio=o18ratio,o17ratio=o17ratio)
    if flag_split_CO_isotopes:
        Atmosphere = archground.models.general.split_CO_isotopes(Atmosphere,c13ratio=c13ratio,o18ratio=o18ratio,o17ratio=o17ratio)

    #Checking that all species exist in Atmosphere and adding them otherwise
    for igas in range(ngas):

        igas_atm = Atmosphere.locate_gas(id_gases[igas],iso_gases[igas])

        if igas_atm == None:
            if vmr_gases[igas] == None:
                print("Gas and isotope IDs = ",id_gases[igas],iso_gases[igas])
                print("Default gases and isotopes in atmosphere = ",Atmosphere.ID,Atmosphere.ISO)
                raise ValueError("error :: vmr_gases must be a real number if the selected gas is not included by default")
            else:
                Atmosphere.add_gas(id_gases[igas],iso_gases[igas],np.ones(Atmosphere.NP)*vmr_gases[igas])
        
        else:
            if vmr_gases[igas] is not None:
                Atmosphere.VMR[:,igas_atm] = vmr_gases[igas]

    #Creating the Spectroscopy class
    Spectroscopy = archground.models.general.create_spectroscopy_class_online(wavenumber,id_gases,iso_gases,
                                                                                delv=0.001,
                                                                                hitran_file=hitran_file,
                                                                                resolving_power=resolving_power)

    #CIA
    #############################################################

    CIA = None 

    #Layer
    #############################################################

    nlay = int(Atmosphere.NP -1)

    Layer = archground.models.general.create_layer_class(Atmosphere,nlay=nlay)

    #Retrieval
    #############################################################

    Retrieval = archground.models.general.create_retrieval_class(-1,philimit=0.1,ncores=1)



    ############################################################################
    #WRITING INPUT FILES
    ############################################################################

    Measurement.write_hdf5(outname)
    Layer.write_hdf5(outname)
    Atmosphere.write_hdf5(outname)
    Telluric.write_hdf5(outname)
    Stellar.write_hdf5(outname)
    Surface.write_hdf5(outname)
    Spectroscopy.write_hdf5(outname)
    Retrieval.write_input_hdf5(outname)
    Scatter.write_hdf5(outname)

    #Writing dummy .apr file
    fapr = open(outname+'.apr','w')
    fapr.write('#Forward model with archground \n')
    fapr.write('\t %i \n' % (1))
    fapr.write('%i \t %i \t %i \n' % (Spectroscopy.ID[0],Spectroscopy.ISO[0],2))
    fapr.write('\t %7.4f \t %7.4f \n' % (1.0,0.1))
    fapr.close()


    ############################################################################
    #RUNNING FORWARD MODEL
    ############################################################################

    #NEMESIS = ans.Retrieval.from_hdf5(outname)
    #SPECONV = NEMESIS.run_forward_model()

    #return NEMESIS.Measurement.VCONV[:,0],SPECONV

