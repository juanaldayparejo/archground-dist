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

#########################################################################################################

def create_measurement_class(waven,meas,errmeas,ephemerides,x0=0.,y0=0.,
                             model_fov=False,
                             nir_channel=True,
                             normalise_spectrum=False,vnorm=None,norm_kernel_size=10):
    """
    Create the measurement class for archNEMESIS 

    Parameters
    ----------

    waven : ndarray (nwave)
        Wavenumber array (cm-1)
    meas : ndarray (nwave)
        Measured spectrum
    errmeas : ndarray (nwave)
        Uncertainty in measured spectrum
    ephemerides : dict
        Dictionary including the ephemerides information from the JPL Horizons System
    x0 : float
        Centre of the science fibre in N-S celestial direction (arcsec)
    y0 : float
        Centre of the science fibre in W-E celestial direction (arcsec)

    Returns
    -------
    
    Measurement : archNEMESIS class
        Instance of the archNEMESIS class

    lat, lon, Ls, LST : float
        Latitude, longitude, local time and subsolar longitude at the centre of the field of view

    """

    from scipy.special import voigt_profile

    #Sorting the wavenumber array
    isort = np.argsort(waven)
    waven = waven[isort]
    meas = meas[isort]
    errmeas = errmeas[isort]
    
    #Creating geometry maps
    res_bin = 0.15
    maps = maps = archground.geometry.build_maps_epoch(ephemerides,res_bin=res_bin)

    #Building measurement class
    Measurement = ans.Measurement_0()
    Measurement.NGEOM = 1
    Measurement.NCONV = np.zeros(Measurement.NGEOM,dtype="int32") + len(waven)

    if normalise_spectrum is True:

        if vnorm is None:
            raise ValueError("error :: if normalising spectrum, vnorm needs to be defined")

        if((vnorm < waven.min()) or (vnorm > waven.max())):
            print("vnorm = ",vnorm)
            print("waven limits = ",waven.min(),waven.max())
            raise ValueError("error :: vnorm must be within the wavenumber limits ")

        #Smoothing the transmission with a 10-pixel box kernel
        kernel = np.ones(norm_kernel_size) / norm_kernel_size
        smoothed = np.convolve(meas, kernel, mode='same')

        anorm = np.interp(vnorm,waven,meas)
        meas /= anorm
        errmeas /= anorm
        Measurement.IFORM = 5       
        Measurement.VNORM = vnorm 

    Measurement.edit_VCONV(waven[:,None])
    Measurement.edit_MEAS(meas[:,None])
    Measurement.edit_ERRMEAS(errmeas[:,None])

    #Spectral resolution
    if nir_channel is True:

        Measurement.FWHM = -1.0

        nwave = 101
        #Construct the NFIL,VFIL,AFIL arrays in the class
        nfil = np.zeros(Measurement.NCONV[0],dtype='int32') + nwave
        vfil = np.zeros((nwave,Measurement.NCONV[0]))
        afil = np.zeros((nwave,Measurement.NCONV[0]))
        for iwave in range(Measurement.NCONV[0]):

            #Parameters listed in Nagel et al. (2023), Section 3.5, Figure 3 and Table 2
            fwhm_gaussian = 1.18e-5 * Measurement.VCONV[iwave,0]
            fwhm_lorentzian = 0.17e-5 * Measurement.VCONV[iwave,0]

            # Parameters for scipy.special.voigt_profile:
            sigma = fwhm_gaussian / np.sqrt(8 * np.log(2))
            gamma = fwhm_lorentzian / 2

            dv = 2.0*fwhm_gaussian
            vwave = np.linspace(-dv,dv,nwave)
            vfil[:,iwave] = Measurement.VCONV[iwave,0] + vwave
            afil[:,iwave] = voigt_profile(vwave, sigma, gamma, out=None)

        Measurement.NFIL = nfil
        Measurement.VFIL = vfil
        Measurement.AFIL = afil

    else:
        raise ValueError("error :: the lineshape for the VIS channel has not been implemented yet")

    #Finding points within the fibre
    r_fibre = 1.5 / 2.
    ipoints = np.where( np.sqrt((maps["hx"] - x0)**2. + (maps["hy"]-y0)**2.) <= r_fibre )
    lats = maps["lat"][ipoints]
    lons = maps["lon"][ipoints]
    lsts = maps["lst"][ipoints]
    emiss_ang = maps["emiss_ang"][ipoints]
    sol_ang = maps["sol_ang"][ipoints]
    azi_ang = maps["azi_ang"][ipoints]

    #Calculating average latitude and longitude
    lat_mean = np.nanmean(lats)
    lon_mean = np.nanmean(lons)
    lst_mean = np.nanmean(lsts)

    #Setting up geometry in Measurement class
    if model_fov is True:
        Measurement.NAV = np.zeros(Measurement.NGEOM,dtype="int32") + len(lats)
        Measurement.edit_FLAT(lats[None,:])
        Measurement.edit_FLON(lons[None,:])
        Measurement.edit_EMISS_ANG(emiss_ang[None,:])
        Measurement.edit_SOL_ANG(sol_ang[None,:])
        Measurement.edit_AZI_ANG(azi_ang[None,:])
        Measurement.edit_WGEOM(np.ones((Measurement.NGEOM,Measurement.NAV[0])))
    else:
        Measurement.NAV = np.zeros(Measurement.NGEOM,dtype="int32") + 1
        Measurement.edit_FLAT(np.ones((Measurement.NGEOM,Measurement.NAV[0])) * np.mean(lats))
        Measurement.edit_FLON(np.ones((Measurement.NGEOM,Measurement.NAV[0])) * np.mean(lons))
        Measurement.edit_EMISS_ANG(np.ones((Measurement.NGEOM,Measurement.NAV[0])) * np.mean(emiss_ang))
        Measurement.edit_SOL_ANG(np.ones((Measurement.NGEOM,Measurement.NAV[0])) * np.mean(sol_ang))
        Measurement.edit_AZI_ANG(np.ones((Measurement.NGEOM,Measurement.NAV[0])) * np.mean(azi_ang))
        Measurement.edit_WGEOM(np.ones((Measurement.NGEOM,Measurement.NAV[0])))

    #Calculating latitude and longitude at centre of the FOV
    Measurement.LATITUDE = lat_mean
    Measurement.LONGITUDE = lon_mean

    #Calculating sub-observer latitude and longitude (planetocentric)
    target = str(ephemerides["targetname"][0]).split("(")[0].strip().upper()   #Name of planet
    equatorial_radius = archground.geometry.planet_dict[target]["equatorial_radius_km"]
    polar_radius = archground.geometry.planet_dict[target]["polar_radius_km"]
    lat_subobs = ephemerides["PDObsLat"][0]        #Planetodetic sub-observer latitude (degrees)
    lon_subobs = (-ephemerides["PDObsLon"][0] + 180.0) % 360.0 - 180.0        #Planetodetic sub-observer longitude (degrees) - East positive
    lat_subobs = archground.geometry.planetodetic_to_planetocentric(
        lat_subobs,
        equatorial_radius=equatorial_radius,
        polar_radius=polar_radius,
    )
    Measurement.SUBOBS_LAT = lat_subobs
    Measurement.SUBOBS_LON = lon_subobs

    #Calculating Doppler shift
    v_planet_observer = ephemerides["delta_rate"][0]
    Measurement.V_DOPPLER = v_planet_observer

    Ls = ephemerides["App_Lon_Sun"][0]

    Measurement.assess()

    return Measurement, lat_mean, lon_mean, Ls, lst_mean 

#########################################################################################################

def create_stellar_class(waven,ephemerides,t_sun=5772.0):
    """
    Create the stellar class for archNEMESIS

    Here, we are going to use a Planck function at 5772 K for the continuum, multiplied by the 
    Toon pseudo-linelist for the solar lines

    Parameters
    ----------

    waven : ndarray (nwave)
        Wavenumber array (cm-1)
    t_sun : float (default = 5772 K)
        Effective temperature of the Sun (K)

    Returns
    -------
    
    Stellar class
    """

    #Getting relevant parameters from ephemerides
    distance_sun_planet = ephemerides["r"][0]   #AU
    velocity_sun_planet = ephemerides["r_rate"][0]  #km/s

    #Reading Toon's linelist
    toon_file = os.path.join(archground.paths.archground_path,'../','data','solar','solar_merged_20240731_600_33300_100.out')
    data = np.loadtxt(toon_file, skiprows=3)
    waven_solar = data[:,0] 
    trans_solar = data[:,1]

    #Getting only relevant part of the spectrum
    waven_solarx = waven_solar[ (waven_solar>=waven.min() - 5.) & (waven_solar<=waven.max() + 5.) ]
    trans_solarx = trans_solar[ (waven_solar>=waven.min() - 5.) & (waven_solar<=waven.max() + 5.) ]
    del waven_solar, trans_solar

    #Calculating Luminosity from a blackbody
    Stellar = ans.Stellar_0()
    Stellar.ISPACE = 0
    Stellar.NWAVE = len(waven_solarx)
    Stellar.WAVE = waven_solarx
    Stellar.RADIUS = 6.957e8 / 1.0e3  # km
    Stellar.DIST = distance_sun_planet
    Stellar.calc_luminosity_blackbody(t_sun)

    Stellar.SOLSPEC *= trans_solarx

    #Applying Doppler shift between star and planet
    c_kms = 299792.458
    doppler_factor = 1.0 + velocity_sun_planet / c_kms
    Stellar.WAVE /= doppler_factor
    Stellar.assess()

    return Stellar

#########################################################################################################

def create_telluric_online(waven,ephemerides,
                            delv=0.001,hitran_file=archground.paths.archnemesis_hitran24,tips_file=archground.paths.archnemesis_tips,
                            era5=False,
                            resolving_power=90000.):
    """
    FUNCTION NAME : create_telluric_online()

    DESCRIPTION : Create an instance of the `Telluric_0` class containing the telluric atmosphere
                   and the spectroscopy calculated from the HITRAN24 and TIPS files

    INPUTS : 

        waven :: Wavenumber array of the measurement (cm-1)
        ephemerides :: Dictionary including the ephemerides calculated from the JPL Horizons System

    OPTIONAL INPUTS:
    
        delv :: Wavenumber step for the line-by-line opacity calculation in cm-1
        hitran_file :: Path to the HITRAN file to be used for the line-by-line opacity calculation (default: archground.paths.archnemesis_hitran24, which is a custom HITRAN24 file containing the main gases in the Venus atmosphere)
        tips_file :: Path to the TIPS file to be used for the partition function calculation (default: archground.paths.archnemesis_tips, which is a custom TIPS file containing the partition function data for the main gases in the Venus atmosphere)
        resolving_power :: Spectral resolving power of the measurement (default: 90000)
        era5 :: If True, then the profiles are downloaded from the ERA-5 reanalysis dataset

    OUTPUTS : 
 
        Telluric :: An instance of the `Telluric_0` class containing the telluric parameters for the Venus case

    CALLING SEQUENCE:

        Telluric = create_telluric_era5_online()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    from datetime import datetime

    #Extracting the relevant parameters from the ephemerides
    datetime_str = ephemerides["datetime_str"][0]
    dt = datetime.strptime(datetime_str, "%Y-%b-%d %H:%M:%S.%f")
    date_obs = dt.strftime("%Y-%m-%d")
    time_obs = dt.strftime("%H:%M:%S")  
    elevation = ephemerides["EL"][0]  #elevation above local horizon (degrees)

    #Extracting the observatory parameters
    lon_obs = archground.geometry.observatory_dict["CAHA"]["lon"]
    lat_obs = archground.geometry.observatory_dict["CAHA"]["lat"]
    alt_obs = archground.geometry.observatory_dict["CAHA"]["elevation"]

    #Defining the Telluric class
    Telluric = ans.Telluric_0()

    #Atmosphere
    ##################################################################################

    #Defining the inputs
    Telluric.DATE=date_obs               #UTC date of the observation
    Telluric.TIME=time_obs               #UTC time of the observation
    Telluric.LATITUDE=lat_obs            #Latitude of the observatory
    Telluric.LONGITUDE=lon_obs           #Longitude of the observatory
    Telluric.ALTITUDE=alt_obs            #Altitude of the observatory
    Telluric.EMISS_ANG=180.-elevation    #Observing angle, looking straight up to the zenith in this case

    if era5 is True:
        #Extracting the atmosphere from the ERA5 model
        Telluric.extract_atmosphere_era5()
    else:
        Telluric.extract_atmosphere_circ()

    Telluric.Atmosphere.NDUST = 0

    #Spectroscopy
    ##################################################################################

    Telluric.Spectroscopy = ans.Spectroscopy_0(ILBL=1)

    ids = [1,2,3,4,5,6,7]
    iso = [0,0,0,0,0,0,0]

    Telluric.Spectroscopy.NGAS = len(ids)
    Telluric.Spectroscopy.ID = ids
    Telluric.Spectroscopy.ISO = iso
    Telluric.Spectroscopy.ISPACE = 0
    Telluric.Spectroscopy.IPROC = np.zeros(Telluric.Spectroscopy.NGAS,dtype='int32')
    Telluric.Spectroscopy.LOCATION_LD = [hitran_file] * Telluric.Spectroscopy.NGAS
    Telluric.Spectroscopy.LOCATION_PF = [tips_file] * Telluric.Spectroscopy.NGAS
    Telluric.Spectroscopy.LOCATION_CD = [hitran_file] * Telluric.Spectroscopy.NGAS
    Telluric.Spectroscopy.LINE_DATA_PARAMS = [ans.MolLineDataParams()] * Telluric.Spectroscopy.NGAS

    #Calculating the spectral grid
    waven_min = waven.min() ; waven_max = waven.max()
    fwhm = np.mean([waven_min,waven_max]) / resolving_power
    waven_minx = waven_min - 5. * fwhm
    waven_maxx = waven_max + 5. * fwhm

    c = 299792458.0   #Speed of light (m/s)
    v_doppler_max = 50.
    waven_minx /= (1.0+v_doppler_max*1.0e3 / c)
    waven_maxx /= (1.0-v_doppler_max*1.0e3 / c)

    wavemin = np.floor(waven_minx/delv)*delv
    wavemax = np.ceil(waven_maxx/delv)*delv
    nwave = int(np.round((wavemax - wavemin) / delv))
    wave = np.linspace( wavemin , wavemax , nwave )
    Telluric.Spectroscopy.NWAVE = nwave
    Telluric.Spectroscopy.WAVE = wave

    return Telluric


