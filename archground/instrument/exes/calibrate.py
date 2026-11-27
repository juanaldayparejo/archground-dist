#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# archGROUND - Python package for ground-based planetary spectroscopy and radiative transfer with archNEMESIS. 
# calibrate.py - Set of functions to process and calibrate SOFIA/EXES data.
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
from datetime import datetime, timezone
from pathlib import Path
import h5py
import archground




##############################################################################################

def process_fits(fitsnames,x_ini_offset=None,y_ini_offset=None,force_slit_pa=None):
    """
    Read the EXES FITS file, extract the data, and compute the relevant ephemerides.

    Parameters
    ----------
    filename : list(str)
        Path to the EXES FITS file.
    x_ini_offset : list(float) (default = 0.0)
        Offset to apply to the EXES image in the x direction (East is negative) (arcsec)
    y_ini_offset : list(float) (default = 0.0)
        Offset to apply to the EXES image in the y direction (East is negative) (arcsec)
    force_slit_pa : list(float) (default = None)
        Force the slit position angle to a specific value (deg)
            
    Returns
    ------- 
    measurement : dict
        Dictionary including the main parameters about the measurement
    ephemerides : dict
        Dictionary including the ephemerides information from JPL Horizons
    """

    from astropy.table import Table, vstack

    if x_ini_offset is None:
        x_ini_offset = np.zeros(len(fitsnames))
    if y_ini_offset is None:
        y_ini_offset = np.zeros(len(fitsnames)) 
    if force_slit_pa is None:
        force_slit_pa = [None] * len(fitsnames)

    #Processing first file to get array sizes
    measurement, ephemerides = process_single_fits(fitsnames[0],x_ini_offset=y_ini_offset[0],y_ini_offset=y_ini_offset[0],force_slit_pa=force_slit_pa[0])

    data = measurement["meas"]
    nscans = data.shape[0]
    npixslit = data.shape[1]
    nwave = data.shape[2]
    norders = data.shape[3]
    nimg = len(fitsnames)
    slit_angle = ephemerides["slit_angle"][0]
    resolving_power = measurement["resolving_power"]

    #Processing all the files with no offsets
    x_vals = np.zeros((nscans,npixslit,nimg))
    y_vals = np.zeros((nscans,npixslit,nimg))
    spectra = np.zeros((nscans,npixslit,nwave,norders,nimg))
    error = np.zeros((nscans,npixslit,nwave,norders,nimg))
    ephemerides_list = []

    ix = 0
    for fitsname in fitsnames:

        measurement, ephemerides = process_single_fits(fitsname,x_ini_offset=x_ini_offset[ix],y_ini_offset=y_ini_offset[ix],force_slit_pa=force_slit_pa[ix])
        waven = measurement["waven"]
        x_vals[:,:,ix] = measurement["x_vals"]
        y_vals[:,:,ix] = measurement["y_vals"]
        spectra[:,:,:,:,ix] = measurement["meas"]
        error[:,:,:,:,ix] = measurement["errmeas"]
        ephemerides_list.append(ephemerides)
        
        #if ephemerides["slit_angle"][0] != slit_angle:
        #    print(ephemerides["slit_angle"],slit_angle)
        #    raise ValueError("error :: the slit angle is different for different observations")

        ix += 1

    #Building new dictionary
    measurement = {
        "nobs": nimg, 
        "x_vals": x_vals,
        "y_vals": y_vals,
        "waven": waven,
        "meas": spectra,
        "errmeas": error,
        "resolving_power": resolving_power
    }

    #Stacking all the ephemerides into a single table
    all_ephemerides = vstack(ephemerides_list, join_type="exact")

    return measurement, all_ephemerides

##############################################################################################

def process_single_fits(fitsname,x_ini_offset=0.0,y_ini_offset=0.0,force_slit_pa=None):
    """
    Read the EXES FITS file, extract the data, and compute the relevant ephemerides.

    Parameters
    ----------
    filename : str
        Path to the EXES FITS file.
    x_ini_offset : float (default = 0.0)
        Offset to apply to the EXES image in the x direction (East is negative) (arcsec)
    y_ini_offset : float (default = 0.0)
        Offset to apply to the EXES image in the y direction (East is negative) (arcsec)
    force_slit_pa : float (default = None)
        Force the slit position angle to a specific value (deg)
            
    Returns
    ------- 
    measurement : dict
        Dictionary including the main parameters about the measurement
    ephemerides : dict
        Dictionary including the ephemerides information from JPL Horizons
    """
    
    from astropy.io import fits
    
    #Reading the FITS file
    ####################################################################################################
    
    hdul = fits.open(fitsname)

    hdr = hdul[0].header

    #for key in hdr.keys():
    #    print(repr(hdr))

    #Observation characteristics
    ######################################################################

    planet = hdr['OBJECT']                       #Observed object
    datetime_obs = hdr["DATE-OBS"]
    date_obs, time_obs = datetime_obs.split("T", maxsplit=1)

    #Position of SOFIA
    sofia_lat = (hdr['LAT_STA']+hdr['LAT_END'])/2.0                     #Latitude of SOFIA at the midpoint of the observation
    sofia_lon = (hdr['LON_STA']+hdr['LON_END'])/2.0                     #Longitude of SOFIA at the midpoint of the observation
    sofia_alt = (hdr['ALTI_STA']+hdr['ALTI_END'])/2.0 * 0.3048          #Altitude of SOFIA at the midpoint of the observation (metres, it is feet on FITS file)

    #Resolving power
    resolving_power = hdr['RP']    #Resolving power of the observation

    #Slit geometry
    telescope_angle = hdr['TELVPA']         #Verical position angle - Degrees East of Mars North
    scan_stepx = hdr['MAPINTX']             #Step of the slit from scan to scan (arcsec) in X (perpendicular to SLIT_HEIGHT)
    scan_stepy = hdr['MAPINTY']             #Step of the slit from scan to scan (arcsec) in Y (perpendicular to SLIT_WIDTH)
    slit_width = hdr['SLTW_ARC']            #Slit width (arcsec)
    slit_width_pix = int(np.round(hdr['SLTW_PIX'],0))     #Slit width (pixels)
    slit_height = hdr['SLTH_ARC']           #Slit height (arcsec)
    slit_height_pix = int(np.round(hdr['SLTH_PIX'],0))    #Slit height (pixels)

    x_step, y_step = extract_values_telescope(hdr["TSTEP"])  #Step of the slit from scan to scan (arcsec) in N and E
    x_ini, y_ini = extract_values_telescope(hdr["OFFSET"])    #Initial position of the slit centre (arcsec) in N and E

    #Detector and diffraction order dimensions
    nx = hdr['NAXIS1']               #Detector size in x-direction
    ny = hdr['NAXIS2']               #Detector size in y-direction
    nscans = hdr['NAXIS3']           #Number of scans in observation
    ndiffor = hdr['NORDERS']              #Number of diffraction orders in observation

    x_diffstart1 = hdr['ORDR_B']      #Bottom pixel of each diffraction order
    x_diffend1 = hdr['ORDR_T']        #Top pixel of each diffraction order

    x_diffstart1 = x_diffstart1.split(',')
    x_diffend1 = x_diffend1.split(',')
    x_diffstart = np.zeros(len(x_diffstart1),dtype='int32')
    x_diffend = np.zeros(len(x_diffend1),dtype='int32')
    for i in range(len(x_diffstart1)):
        x_diffstart[i] = int(x_diffstart1[i])
    for i in range(len(x_diffend1)):
        x_diffend[i] = int(x_diffend1[i])

    nslitpix = np.unique(x_diffend-x_diffstart)    #Number of pixels in an order
    if len(nslitpix) > 1:
        raise ValueError("The number of pixels in the diffraction orders is not constant. Please check the FITS file.")
    else:
        nslitpix = nslitpix[0]
    plate_scale = hdr['PLTSCALE']                #arcseconds per pixel after focal reduction
    
    #Reading the data
    #########################################################################

    MEAS1 = hdul[0].data         #ergs s-1 cm-2 (cm-1)-1 sr-1  (nscans,nx,ny)
    ERRMEAS1 = hdul[1].data      #ergs s-1 cm-2 (cm-1)-1 sr-1  (nscans,nx,ny)
    VCONV1 = hdul[6].data        #cm-1  (nx,ny)
    SPATIAL1 = hdul[7].data      #arcsec 

    hdul.close()
    
    #Re-shaping the data 
    ###########################################################################

    #Images in the detector frame
    MEAS = np.transpose(MEAS1[:,:,:],axes=(1,2,0))        #(nx,ny,nscans)
    ERRMEAS= np.transpose(ERRMEAS1[:,:,:],axes=(1,2,0))   #(nx,ny,nscans)
    VCONV = VCONV1[:,:]                                   #(nx,ny)
    SPATIAL = SPATIAL1[:,:]                               #(nx,ny)

    #Coverting the detector images to spectra
    nconv = nx
    MEASX = np.zeros((nconv,nslitpix,nscans,ndiffor))
    ERRMEASX = np.zeros((nconv,nslitpix,nscans,ndiffor))
    VCONVX = np.zeros((nconv,ndiffor))

    for iDIFFOR in range(ndiffor):
        VCONVX[:,iDIFFOR] = VCONV[x_diffstart[iDIFFOR],:]
        for iSCAN in range(nscans):
            MEASX[:,:,iSCAN,iDIFFOR] = np.transpose(MEAS[x_diffstart[iDIFFOR]:x_diffend[iDIFFOR],:,iSCAN])
            ERRMEASX[:,:,iSCAN,iDIFFOR] = np.transpose(ERRMEAS[x_diffstart[iDIFFOR]:x_diffend[iDIFFOR],:,iSCAN])


    #Creating mask for assuring no nan values
    valid = ~np.isnan(VCONVX).any(axis=1)
    inotnan = np.flatnonzero(valid)

    # Apply the same pixel selection to all orders
    VCONVX = VCONVX[valid, :]
    MEASX = MEASX[valid, :, :, :]
    ERRMEASX = ERRMEASX[valid, :, :, :]

    #Reshaping so that they are in the format (nscans,nslitpix,nconv,norders)
    data = np.transpose(MEASX,axes=(2,1,0,3))
    noise = np.transpose(ERRMEASX,axes=(2,1,0,3))
    waven = VCONVX

    #Calculating the geometry of the observation from JPL horizons
    ####################################################################################################
    
    ephemerides = archground.get_ephemerides_epoch(date_obs, time_obs, planet, "CUSTOM",
                                                  latitude=sofia_lat, longitude=sofia_lon, elevation=sofia_alt/1000.)

    #Computing first guess of the geometry
    ####################################################################################################
    
    pa = np.deg2rad(telescope_angle)
    dx_slit = -np.sin(pa)
    dy_slit =  np.cos(pa)

    # Edge at the negative end of the slit, relative to its centre
    x_edge = x_ini - (slit_height / 2) * dx_slit
    y_edge = y_ini - (slit_height / 2) * dy_slit

    # Slit centre for each scan
    scan = np.arange(nscans)
    x0_scans = x_ini + scan * x_step
    y0_scans = y_ini + scan * y_step

    # Pixel-centre offsets along the slit
    s = (np.arange(nslitpix) - (nslitpix - 1) / 2) * plate_scale

    x_vals = x0_scans[:, None] + s * dx_slit
    y_vals = y0_scans[:, None] + s * dy_slit

    #Applying artificial offset if required
    x_vals += x_ini_offset
    y_vals += y_ini_offset

    #Removing the first and last pixel of the slit
    x_vals = x_vals[:,1:nslitpix-2]
    y_vals = y_vals[:,1:nslitpix-2]
    data = data[:,1:nslitpix-2,:,:]
    noise = noise[:,1:nslitpix-2,:,:]
    
    measurement = {
        "x_vals": x_vals,
        "y_vals": y_vals,
        "waven": waven,
        "meas": data,
        "errmeas": noise,
        "resolving_power": resolving_power,
    }

    #Adding geometry information to the ephemerides
    ephemerides["x0_offset"] = x_ini_offset
    ephemerides["y0_offset"] = y_ini_offset
    ephemerides["slit_angle"] = telescope_angle

    return measurement,ephemerides


##############################################################################################

def print_header_fits(fitsname):
    """
    Read the EXES FITS file, and print the header

    Parameters
    ----------
    filename : str
        Path to the EXES FITS file.
            
    Returns
    ------- 
    measurement : dict
        Dictionary including the main parameters about the measurement
    ephemerides : dict
        Dictionary including the ephemerides information from JPL Horizons
    """
    
    from astropy.io import fits
    
    hdul = fits.open(fitsname)
    hdr = hdul[0].header
    for key in hdr.keys():
        print(repr(hdr))
    hdul.close()

##############################################################################################

def construct_continuum_image(measurement,v0,vwidth=0.1):
    """
    Function to compute a continuum image using the spectral window specified

    Parameters
    ----------

    measurement : dict
        Dictionary including the spectra from the EXES observations
    v0 : float
        Wavenumber at which to compute the continuum (cm-1)
    vwidth : float
        Width of the window to use to compute the continuum (cm-1)

    Returns
    -------

    image_continuum : ndarray (nscans,nslitpix,nobs)
    """

    nobs = measurement["nobs"]
    nwave = measurement["waven"].shape[0]
    norder = measurement["waven"].shape[1]
    waven = measurement["waven"]

    #Finding the order corresponding to v0
    isel = np.nan
    for iorder in range(norder):
        if( (waven[:,iorder].min()<=v0) & (waven[:,iorder].max()>=v0) ):
            isel = iorder
    if isel == np.nan:
        raise ValueError("error :: no orders spectral the spectral range of v0")
    iorder = int(isel)

    icont = np.where( (waven[:,iorder]>=(v0-vwidth/2.)) & (waven[:,iorder]<=(v0+vwidth/2.)) )[0]

    image_continuum = np.mean(measurement["meas"][:,:,icont,iorder,:],axis=2)  #(nscans,nslitpix,nimg)

    return image_continuum

##############################################################################################

def construct_linedepth_image(measurement,line_window,cont1_window,cont2_window):
    """
    Function to compute a line depth image using the spectral window specified

    Parameters
    ----------
    measurement : dict
        Dictionary including the spectra from the EXES observations
    line_window : tuple(float,float)
        Lower and upper bounds of the window at which to compute the line depth (cm-1)
    cont1_window : tuple(float,float)
        Lower and upper bounds of the window at which to compute the continuum on the left (cm-1)
    cont2_window : tuple(float,float)
        Lower and upper bounds of the window at which to compute the continuum on the right (cm-1)

    Returns
    -------
    image_linedepth : ndarray (nscans,nslitpix,nobs)
    """

    nobs = measurement["nobs"]
    nwave = measurement["waven"].shape[0]
    norder = measurement["waven"].shape[1]
    waven = measurement["waven"]

    #Finding the order corresponding to v0
    v0 = np.mean(line_window)
    isel = np.nan
    for iorder in range(norder):
        if( (waven[:,iorder].min()<=v0) & (waven[:,iorder].max()>=v0) ):
            isel = iorder
    if isel == np.nan:
        raise ValueError("error :: no orders spectral the spectral range of v0")
    iorder = int(isel)

    wave = waven[:,iorder]
    iblue = (wave >= cont1_window[0]) & (wave <= cont1_window[1])
    ired  = (wave >= cont2_window[0])  & (wave <= cont2_window[1])
    iline = (wave >= line_window[0]) & (wave <= line_window[1])

    if not np.any(iblue) or not np.any(ired) or not np.any(iline):
        raise ValueError(
            "At least one selected window contains no wavelength samples."
        )

    # Representative wavelength of each continuum window
    wave_blue = np.nanmean(wave[iblue])
    wave_red = np.nanmean(wave[ired])

    # Mean flux in each continuum window for every spectrum
    flux_blue = np.nanmean(
        measurement["meas"][:, :, iblue, iorder, :],
        axis=2
    )  # (nscans, nslitpix, nimg)

    flux_red = np.nanmean(
        measurement["meas"][:, :, ired, iorder, :],
        axis=2
    )  # (nscans, nslitpix, nimg)

    # Continuum evaluated at every wavelength in the line window.
    # The added singleton dimension corresponds to wavelength.
    wave_line = wave[iline]
    fraction = (
        (wave_line - wave_blue) /
        (wave_red - wave_blue)
    )[None, None, :, None]

    continuum_line = (
        flux_blue[:, :, None, :]
        + (flux_red - flux_blue)[:, :, None, :] * fraction
    )

    flux_line = measurement["meas"][:, :, iline, iorder, :]

    # Positive EW for absorption and negative EW for emission
    #normalized_deficit = 1.0 - flux_line / continuum_line
    normalized_deficit = continuum_line - flux_line 

    image_continuum = np.trapezoid(
        normalized_deficit,
        x=wave_line,
        axis=2
    )  # (nscans, nslitpix, nimg)

    # Mask divisions by zero or invalid continuum estimates
    invalid = (
        ~np.all(np.isfinite(continuum_line), axis=2)
        | np.any(continuum_line == 0, axis=2)
    )
    image_continuum[invalid] = np.nan

    return image_continuum


import numpy as np


def construct_equivalent_width_image(
    measurement,
    line_window,
    cont1_window,
    cont2_window,
    *,
    min_continuum_snr=3.0,
):
    """Compute a continuum-normalized equivalent-width image in cm⁻¹.

    Each window is a pair of (lower, upper) wavenumber bounds in cm⁻¹.
    measurement["errmeas"] contains per-channel 1σ radiance uncertainties
    and has the same shape as measurement["meas"].

    The output has shape (nscans, nslitpix, nobs), with positive values
    for absorption and negative values for emission. A pixel returns NaN
    if any line sample has an invalid flux or a continuum below the S/N
    threshold. Uncertainty propagation assumes independent channels.
    """
    if not np.isfinite(min_continuum_snr) or min_continuum_snr < 0:
        raise ValueError("min_continuum_snr must be finite and non-negative.")

    waven = measurement["waven"]
    v0 = np.mean(line_window)

    # Select the last diffraction order that contains the line centre.
    iorder = None
    for order in range(waven.shape[1]):
        wave_order = waven[:, order]
        finite_wave = wave_order[np.isfinite(wave_order)]
        if finite_wave.size and finite_wave.min() <= v0 <= finite_wave.max():
            iorder = order

    if iorder is None:
        raise ValueError(f"No diffraction order contains {v0:g} cm⁻¹.")

    wave = waven[:, iorder]
    spectra = measurement["meas"][:, :, :, iorder, :]
    errors = measurement["errmeas"][:, :, :, iorder, :]

    icont1 = (wave >= cont1_window[0]) & (wave <= cont1_window[1])
    icont2 = (wave >= cont2_window[0]) & (wave <= cont2_window[1])
    iline = (wave >= line_window[0]) & (wave <= line_window[1])

    if not icont1.any() or not icont2.any():
        raise ValueError("Each continuum window must contain at least one sample.")
    if np.any(icont1 & icont2):
        raise ValueError("The continuum windows must use disjoint samples.")
    if np.count_nonzero(iline) < 2:
        raise ValueError("The line window must contain at least two samples.")

    def continuum_anchor(window):
        """Return the mean wavenumber, mean radiance, and its uncertainty."""
        flux = spectra[:, :, window, :]
        sigma = errors[:, :, window, :]
        wave_samples = wave[window][None, None, :, None]
        valid = np.isfinite(flux) & np.isfinite(sigma) & (sigma > 0)

        count = valid.sum(axis=2)
        wave_sum = np.sum(np.where(valid, wave_samples, 0.0), axis=2)
        flux_sum = np.sum(np.where(valid, flux, 0.0), axis=2)
        variance_sum = np.sum(np.where(valid, sigma, 0.0) ** 2, axis=2)

        mean_wave = np.full(count.shape, np.nan, dtype=float)
        mean_flux = np.full(count.shape, np.nan, dtype=float)
        mean_sigma = np.full(count.shape, np.nan, dtype=float)

        np.divide(wave_sum, count, out=mean_wave, where=count > 0)
        np.divide(flux_sum, count, out=mean_flux, where=count > 0)
        np.divide(np.sqrt(variance_sum), count, out=mean_sigma, where=count > 0)

        return mean_wave, mean_flux, mean_sigma

    wave_cont1, flux_cont1, sigma_cont1 = continuum_anchor(icont1)
    wave_cont2, flux_cont2, sigma_cont2 = continuum_anchor(icont2)

    # Integrate in increasing wavenumber so absorption gives positive EW.
    line_indices = np.flatnonzero(iline)
    line_indices = line_indices[np.argsort(wave[line_indices])]
    wave_line = wave[line_indices]
    flux_line = spectra[:, :, line_indices, :]

    # Interpolate the continuum and propagate its uncertainty.
    wave_offset = wave_line[None, None, :, None] - wave_cont1[:, :, None, :]
    anchor_separation = (wave_cont2 - wave_cont1)[:, :, None, :]
    weight2 = np.full(flux_line.shape, np.nan, dtype=float)
    np.divide(
        wave_offset,
        anchor_separation,
        out=weight2,
        where=np.isfinite(anchor_separation) & (anchor_separation != 0),
    )
    weight1 = 1.0 - weight2

    continuum_line = (
        weight1 * flux_cont1[:, :, None, :]
        + weight2 * flux_cont2[:, :, None, :]
    )
    sigma_continuum_line = np.hypot(
        weight1 * sigma_cont1[:, :, None, :],
        weight2 * sigma_cont2[:, :, None, :],
    )

    valid = (
        np.isfinite(flux_line)
        & np.isfinite(continuum_line)
        & np.isfinite(sigma_continuum_line)
        & (continuum_line > 0)
        & (continuum_line >= min_continuum_snr * sigma_continuum_line)
    )

    relative_flux = np.full(flux_line.shape, np.nan, dtype=float)
    np.divide(flux_line, continuum_line, out=relative_flux, where=valid)

    return np.trapezoid(1.0 - relative_flux, x=wave_line, axis=2)

##############################################################################################

def calibrate_geometry(x_vals,y_vals,image,ephemerides,
                        fitting_method='uniform',
                        max_radiance=None,
                        max_off=4.,res_off=0.1,
                        hx=None,hy=None,model_image=None,
                        convol=True,
                        fwhm_psf=1.0,
                        MakePlot=False):
    """
    Function to calibrate the geometry of the EXES data

    Parameters
    ----------
    xvals : ndarray (nacq,nslitpix)
        Spatial position of each pixel of the image in the x-direction (arcsec)
    yvals : ndarray (nacq,nslitpix)
        Spatial position of each pixel of the image in the y-direction (arcsec) 
    image : ndarray (nacq,nslitpix)
        Measured image that will be used for the calibration
    max_radiance : float
        If not 

    Returns
    -------
    x_ini_offset : float
        Offset to be applied to the location of the pixels in the x-direction (arcsec)
    y_ini_offset : float
        Offset to be applied to the location of the pixels in the y-direction (arcsec)
    """

    from astropy.convolution import convolve, Gaussian2DKernel
    from scipy.interpolate import RectBivariateSpline
    
    #Calculating the geometry of the observation using the geometry parameters
    ####################################################################################################
    
    angsize = np.atleast_1d(ephemerides["ang_width"])[0]
    if fitting_method != 'custom':
        angmax = float(int(np.atleast_1d(ephemerides["ang_width"])[0])) ; res = 0.1
        maps = archground.geometry.build_maps_epoch(ephemerides,res_bin=0.1,angmax=angmax)
        hx = maps["hx"] ; hy = maps["hy"] ; sol_ang = maps["sol_ang"]
            
    #Creating the model to fit the geometry
    ##############################################################################################################

    if fitting_method == 'sza':
        
        bsurf = np.cos( sol_ang / 180. * np.pi )
        bsurf[np.isnan(bsurf)==True] = 0.0
        
    elif fitting_method == 'uniform':
        
        sza = np.cos( sol_ang / 180. * np.pi )
        bsurf = np.ones(sza.shape)
        bsurf[np.isnan(sza)==True] = 0.0

    elif fitting_method == 'custom':

        if model_image is None or hx is None or hy is None:
            raise ValueError("For fitting_method='custom', you must provide a hx, hy and model_image.")
        bsurf = model_image


    if convol is True:
        
        #Convolving the image with a 2D gaussian
        fwhm_exes = fwhm_psf  #arcsec
        sigma_exes = fwhm_exes / 2.355  # σ = FWHM / 2.355
        
        # Convert sigma from arcsec to pixels
        sigma_pixels = sigma_exes / res
        
        # Create 2D Gaussian kernel
        gaussian_kernel = Gaussian2DKernel(x_stddev=sigma_pixels)
        
        #Convolving image
        bsurf_conv = convolve(bsurf, gaussian_kernel, normalize_kernel=True)
        
    else:
        
        bsurf_conv = bsurf    
        
    if max_radiance is not None:
        bsurf_conv = bsurf_conv / bsurf_conv.max() * max_radiance
    
        
    if MakePlot is True:
        
        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(8,3))
        
        cmap = 'turbo'
        
        im = ax1.contourf(hx,hy,bsurf*1.0e7,cmap=cmap,levels=21)
        ax1.set_xlabel('Angular distance (arcsec)')
        ax1.set_ylabel('Angular distance (arcsec)')
        Drawing_uncolored_circle = plt.Circle( (0.0, 0.0), angsize/2. ,fill = False, color='white' )
        ax1.add_artist( Drawing_uncolored_circle )
        ax1.set_title('Model image')
        ax1.grid()
        ax1.set_facecolor('lightgray')
        
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label('Radiance (erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$ (cm$^{-1}$)$^{-1}$)')
        
        im = ax2.contourf(hx,hy,bsurf_conv*1.0e7,cmap=cmap,levels=21)
        ax2.set_xlabel('Angular distance (arcsec)')
        ax2.set_ylabel('Angular distance (arcsec)')
        Drawing_uncolored_circle = plt.Circle( (0.0, 0.0), angsize/2. ,fill = False, color='white' )
        ax2.add_artist( Drawing_uncolored_circle )
        ax2.set_title('Model image')
        ax2.grid()
        ax2.set_facecolor('lightgray')
        
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label('Radiance (erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$ (cm$^{-1}$)$^{-1}$)')
                
        plt.tight_layout()
        
        
    #Fitting the geometry of the EXES image
    ##########################################################################################
    
    xoff = np.arange(-max_off,max_off+res_off,res_off)
    yoff = np.arange(-max_off,max_off+res_off,res_off)

    x = hx[0,:]
    y = hy[:,0]

    s = RectBivariateSpline(x,y,bsurf_conv.T)

    rms = np.zeros((len(xoff),len(yoff)))
    for ioff in range(len(xoff)):
        for joff in range(len(yoff)):

            fm = np.zeros(x_vals.shape)
            for i in range(x_vals.shape[0]):
                for j in range(x_vals.shape[1]):
            
                    fm[i,j] = s(x_vals[i,j]+xoff[ioff],y_vals[i,j]+yoff[joff])[0,0]

            if max_radiance is not None:
                rms[ioff,joff] = np.sum( (fm - image)**2. )
            else:
                rms[ioff,joff] = np.sum( (fm/fm.max() - image/image.max())**2. )
    
    imax,jmax = np.unravel_index(rms.argmin(), rms.shape)
    x_ini_offset = xoff[imax] ; y_ini_offset = yoff[jmax]
    
    if MakePlot==True:
        
        fig,ax1 = plt.subplots(1,1,figsize=(3,3))
        ax1.contourf(xoff,yoff,rms.T,levels=101,cmap='turbo')
        imax,jmax = np.unravel_index(rms.argmin(), rms.shape)
        ax1.scatter(xoff[imax],yoff[jmax],c='tab:red')
        ax1.set_xlabel('x offset (arcsec)')
        ax1.set_ylabel('y offset (arcsec)')
        plt.tight_layout()
        

        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(6,3))

        psize = 5.

        ax1.scatter(x_vals+x_ini_offset,y_vals+y_ini_offset,s=psize,c=image,cmap='turbo',vmin=0.0)
        Drawing_uncolored_circle = plt.Circle( (0.0, 0.0), angsize/2. ,fill = False, color='white' )
        ax1.add_artist( Drawing_uncolored_circle )
        ax1.set_title('EXES map')
        ax1.set_xlabel('Angular distance (arcsec)')
        ax1.set_ylabel('Angular distance (arcsec)')
        ax1.grid()
        archground.geometry.add_compass(ax1)

        lim = angsize
        ax1.set_xlim(-lim,lim)
        ax1.set_ylim(-lim,lim)
        ax1.set_facecolor('lightgray')

        fm = np.zeros(x_vals.shape)
        for i in range(x_vals.shape[0]):
            for j in range(x_vals.shape[1]):
                fm[i,j] = s(x_vals[i,j]+x_ini_offset,y_vals[i,j]+y_ini_offset)[0,0]

        ax2.scatter(x_vals+x_ini_offset,y_vals+y_ini_offset,s=psize,c=fm,cmap='turbo',vmin=0.0)
        Drawing_uncolored_circle = plt.Circle( (0.0, 0.0), angsize/2. ,fill = False, color='white' )
        ax2.add_artist( Drawing_uncolored_circle )
        ax2.set_title('Modelled map')
        ax2.set_xlabel('Angular distance (arcsec)')
        ax2.set_ylabel('Angular distance (arcsec)')
        ax2.grid()
        archground.geometry.add_compass(ax2)

        lim = angsize
        ax2.set_xlim(-lim,lim)
        ax2.set_ylim(-lim,lim)
        ax2.set_facecolor('lightgray')

        plt.tight_layout()
    
    return x_ini_offset,y_ini_offset
    
########################################################################################################################

def save_geometry_calibration_results(fitsname,x_ini_offset,y_ini_offset):
    """
    FUNCTION NAME : save_geometry_calibration_results()

    DESCRIPTION : Save the results of the geometry calibration in a text file

    INPUTS : 

        fitsname :: Name of the FITS file

    OPTIONAL INPUTS:
    
        x_ini_offset :: Offset to apply to the EXES image in the x direction (arcsec)
        y_ini_offset :: Offset to apply to the EXES image in the y direction (arcsec)
            
    OUTPUTS : 
 
        None

    CALLING SEQUENCE:

        save_geometry_calibration_results(fitsname,x_ini_offset,y_ini_offset)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    """
    
    #Saving the results in a text file
    with open(fitsname+'.txt', 'w') as f:
        f.write('x_ini_offset, y_ini_offset\n')
        f.write('{:.6f}, {:.6f}\n'.format(x_ini_offset, y_ini_offset))

########################################################################################################################

def read_geometry_calibration_results(fitsname):
    """
    FUNCTION NAME : read_geometry_calibration_results()

    DESCRIPTION : Read the results of the geometry calibration from a text file

    INPUTS : 

        fitsname :: Name of the FITS file

    OPTIONAL INPUTS:
    
        None
        
    OUTPUTS : 

        x_ini_offset :: Offset to apply to the EXES image in the x direction (arcsec)
        y_ini_offset :: Offset to apply to the EXES image in the y direction (arcsec)

    CALLING SEQUENCE:

        x_ini_offset,y_ini_offset = read_geometry_calibration_results(fitsname)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    """
    
    #Reading the results from a text file
    with open(fitsname+'.txt', 'r') as f:
        lines = f.readlines()
        x_ini_offset, y_ini_offset = map(float, lines[1].strip().split(','))
    
    return x_ini_offset, y_ini_offset





##############################################################################################
##############################################################################################
#                                     COMPOSITE IMAGE
##############################################################################################
##############################################################################################

def calculate_composite_image(measurement,ephemerides,bin_res=0.5):
    """
    Function to calculate a composite disk image by averaging many EXES maps

    Parameters
    ----------

    measurement : dict
        Dictionary including the spectra from the EXES observations
    ephemerides : dict
        Dictionary including the ephemerides information for the EXES observations
    bin_res : float
        Spatial resolution of the composite image (arcsec) 

    Returns
    -------

    composite_measurement : dict
        Dictionary including the composite spectral map
    composite_ephemerides : dict
        Dictionary including the ephemerides at the mid-time
    composite_geometry : dict
        Dictionary includijng the geometry of the observation at the mid-time
    """

    from astropy.time import Time

    #Calculating the composite image
    ##############################################################################

    #Calculating the composite x-y grid
    x_res = bin_res
    y_res = bin_res
    angsize = np.mean(ephemerides['ang_width'][:])
    angmax = np.ceil(angsize)

    # Create pixel centers with zero exactly included
    nx = int(np.ceil(angsize / x_res))
    ny = int(np.ceil(angsize / y_res))

    x_grid = np.arange(-nx, nx + 1, dtype=float) * x_res
    y_grid = np.arange(-ny, ny + 1, dtype=float) * y_res

    # Pixel boundaries
    x_edges = np.r_[x_grid - x_res / 2,
                    x_grid[-1] + x_res / 2]

    y_edges = np.r_[y_grid - y_res / 2,
                    y_grid[-1] + y_res / 2]

    # spectra shape:
    # (original_x, original_y, spectral_dim_1, spectral_dim_2, nimg)
    spectral_shape = measurement["meas"].shape[2:4]

    # Output shape:
    # (new_y, new_x, spectral_dim_1, spectral_dim_2)
    output_shape = (
        len(y_grid),
        len(x_grid),
        *spectral_shape,
    )

    spectra_sum = np.zeros(output_shape, dtype=float)
    error_variance_sum = np.zeros(output_shape, dtype=float)
    spectra_count = np.zeros(output_shape, dtype=np.int64)

    limpix_y = 2
    nimg = measurement["meas"].shape[-1]
    for ix in range(nimg):
        x = measurement["x_vals"][:, limpix_y:-limpix_y, ix].ravel()
        y = measurement["y_vals"][:, limpix_y:-limpix_y, ix].ravel()

        current_spectra = measurement["meas"][
            :, limpix_y:-limpix_y, :, :, ix
        ]
        
        #current_error = np.broadcast_to(
        #    measurement["errmeas"][None, limpix_y:-limpix_y, :, :, ix],
        #    current_spectra.shape,
        #)
        current_error = measurement["errmeas"][
            :, limpix_y:-limpix_y, :, :, ix
        ]

        current_spectra = current_spectra.reshape(
            -1,
            *spectral_shape,
        )
        
        current_error = current_error.reshape(
            -1,
            *spectral_shape,
        )

        # Find the output pixel for each original spatial point
        x_bin = np.searchsorted(x_edges, x, side="right") - 1
        y_bin = np.searchsorted(y_edges, y, side="right") - 1

        valid_position = (
            np.isfinite(x)
            & np.isfinite(y)
            & (x_bin >= 0)
            & (x_bin < len(x_grid))
            & (y_bin >= 0)
            & (y_bin < len(y_grid))
        )

        x_bin = x_bin[valid_position]
        y_bin = y_bin[valid_position]
        current_spectra = current_spectra[valid_position]
        current_error = current_error[valid_position]

        # Account for possible NaNs within individual spectra
        valid_values = (
            np.isfinite(current_spectra)
            & np.isfinite(current_error)
        )

        np.add.at(
            spectra_sum,
            (y_bin, x_bin),
            np.where(valid_values, current_spectra, 0),
        )

        np.add.at(
            error_variance_sum,
            (y_bin, x_bin),
            np.where(
                valid_values,
                current_error**2,
                0.0,
            ),
        )

        np.add.at(
            spectra_count,
            (y_bin, x_bin),
            valid_values,
        )

    # Mean spectrum in every new spatial pixel
    composite_spectra = np.full(
        output_shape,
        np.nan,
        dtype=float,
    )

    composite_error = np.full(
        output_shape,
        np.nan,
        dtype=float,
    )

    has_data = spectra_count > 0

    # Arithmetic mean
    np.divide(
        spectra_sum,
        spectra_count,
        out=composite_spectra,
        where=has_data,
    )

    # Error of the arithmetic mean:
    # sqrt(sigma_1² + ... + sigma_N²) / N
    np.divide(
        np.sqrt(error_variance_sum),
        spectra_count,
        out=composite_error,
        where=has_data,
    )

    composite_measurement = {
        "orders" : np.arange(0,measurement["waven"].shape[1]),
        "waven" : measurement["waven"],
        "x_grid" : x_grid,
        "y_grid" : y_grid,
        "meas" : composite_spectra,
        "errmeas" : composite_error,
        "resolving_power" : measurement["resolving_power"]
    }

    #Calculating the ephemerides at the mid-time
    ##############################################################################

    #Calculating mean date and time
    times = Time(ephemerides["datetime_jd"], format="jd", scale="utc")
    mean_datetime = times.mean()
    mean_date, mean_time = mean_datetime.iso.split()

    #Extracting the planet information
    planet = ephemerides["targetname"][0].split(" (")[0]
 
    composite_ephemerides = archground.geometry.get_ephemerides_epoch(mean_date,mean_time,planet,"IRTF")

    #Calculating the geometry maps
    ###############################################################################

    composite_maps = archground.geometry.build_maps_epoch(composite_ephemerides,res_bin=bin_res,angmax=x_grid.max())

    return composite_measurement, composite_ephemerides, composite_maps


def _decode_hdf5_value(value):
    """Convert HDF5 byte strings into normal Python strings."""

    if isinstance(value, bytes):
        return value.decode("utf-8")

    if isinstance(value, np.ndarray) and value.dtype.kind == "S":
        return np.char.decode(value, "utf-8")

    return value


#######################################################################################

def _read_attributes(hdf5_object):
    """Read all attributes attached to a group or dataset."""

    return {
        name: _decode_hdf5_value(value)
        for name, value in hdf5_object.attrs.items()
    }


#######################################################################################

def _read_group(group):
    """Read the datasets and metadata from one HDF5 group."""

    data = {}
    dataset_attributes = {}

    for name, item in group.items():

        if not isinstance(item, h5py.Dataset):
            continue

        # asstr() handles scalar and array string datasets.
        if h5py.check_string_dtype(item.dtype) is not None:
            value = item.asstr()[()]
        else:
            value = item[()]

        data[name] = _decode_hdf5_value(value)
        dataset_attributes[name] = _read_attributes(item)

    return {
        "data": data,
        "attributes": _read_attributes(group),
        "dataset_attributes": dataset_attributes,
    }


#######################################################################################

def create_composite_hdf5(
    output_filename,
    composite_measurement,
    composite_ephemerides,
    composite_geometry,
    input_files=None,
):
    """
    Function to write the composite image information into a HDF5 file

    Parameters
    ----------

    composite_measurement : dict
        Dictionary including the composite spectral map from the TEXES observations
    composite_ephemerides : dict
        Dictionary including the ephemerides information for the composite of the TEXES observations
    composite_geometry : dict
        Dictionary including the geometry information for the composite of the TEXES observations
    input_files : list(str)
        List of files that were used to generate the composite image (for the records)

    Returns
    -------
        HDF5 file including all the information
    """

    from astropy.time import Time
    
    output_filename = Path(output_filename)
    output_filename.parent.mkdir(parents=True, exist_ok=True)

    spectra = np.asarray(composite_measurement["meas"])
    error = np.asarray(composite_measurement["errmeas"])
    x_grid = np.asarray(composite_measurement["x_grid"])
    y_grid = np.asarray(composite_measurement["y_grid"])
    waven = np.asarray(composite_measurement["waven"])
    orders = np.asarray(composite_measurement["orders"])
    resolving_power = composite_measurement["resolving_power"]

    # Validate the main dimensions.
    if spectra.ndim != 4:
        raise ValueError(
            "spectra must have dimensions "
            "(y, x, spectral_pixel, order)"
        )

    ny, nx, nwave, norders = spectra.shape

    if x_grid.shape != (nx,):
        raise ValueError(
            f"x_grid has shape {x_grid.shape}; expected {(nx,)}"
        )

    if y_grid.shape != (ny,):
        raise ValueError(
            f"y_grid has shape {y_grid.shape}; expected {(ny,)}"
        )

    if waven.shape != (nwave, norders):
        raise ValueError(
            f"waven has shape {waven.shape}; "
            f"expected {(nwave, norders)}"
        )

    if orders.shape != (norders,):
        raise ValueError(
            f"orders has shape {orders.shape}; "
            f"expected {(norders,)}"
        )

    # Validate all per-pixel geometry maps.
    expected_map_shape = (ny, nx)

    for name, values in composite_geometry.items():
        values = np.asarray(values)

        if values.shape != expected_map_shape:
            raise ValueError(
                f"composite_geometry['{name}'] has shape {values.shape}; "
                f"expected {expected_map_shape}"
            )

    with h5py.File(output_filename, "w") as h5:

        # -------------------------------------------------------------
        # File-level metadata
        # -------------------------------------------------------------
        h5.attrs["title"] = "Composite TEXES spectral map"
        h5.attrs["product_type"] = "composite_spectral_map"
        h5.attrs["format_version"] = "1.0"
        h5.attrs["creation_time_utc"] = datetime.now(
            timezone.utc
        ).isoformat()

        h5.attrs["spectra_dimension_order"] = (
            "y, x, spectral_pixel, order"
        )

        h5.attrs["geometry_map_dimension_order"] = "y, x"

        # -------------------------------------------------------------
        # Science data and coordinates
        # -------------------------------------------------------------
        science = h5.create_group("science")

        ds = science.create_dataset(
            "meas",
            data=spectra.astype(np.float32),
            compression="gzip",
            compression_opts=4,
            shuffle=True,
            chunks=True,
        )
        ds.attrs["long_name"] = "Spectral radiance"
        ds.attrs["units"] = (
            "erg s-1 cm-2 sr-1 (cm-1)-1"
        )
        ds.attrs["dimensions"] = (
            "y, x, spectral_pixel, order"
        )
        ds.attrs["missing_value"] = np.nan


        ds = science.create_dataset(
            "errmeas",
            data=error.astype(np.float32),
            compression="gzip",
            compression_opts=4,
            shuffle=True,
            chunks=True,
        )
        ds.attrs["long_name"] = "Uncertainty in spectral radiance"
        ds.attrs["units"] = (
            "erg s-1 cm-2 sr-1 (cm-1)-1"
        )
        ds.attrs["dimensions"] = (
            "y, x, spectral_pixel, order"
        )
        ds.attrs["missing_value"] = np.nan
        
        
        ds = science.create_dataset(
            "waven",
            data=waven.astype(np.float64),
            compression="gzip",
            compression_opts=4,
            shuffle=True,
        )
        ds.attrs["long_name"] = "Spectral wavenumber"
        ds.attrs["units"] = "cm-1"
        ds.attrs["dimensions"] = "spectral_pixel, order"

        ds = science.create_dataset(
            "order",
            data=orders,
        )
        ds.attrs["long_name"] = "Diffraction-order index"
        ds.attrs["units"] = "1"
        ds.attrs["dimensions"] = "order"

        ds = science.create_dataset(
            "resolving_power",
            data=resolving_power,
        )
        ds.attrs["long_name"] = "Resolving power"
        ds.attrs["units"] = "1"
        ds.attrs["dimensions"] = ""

        # -------------------------------------------------------------
        # Per-pixel geometry stored alongside the science data
        # -------------------------------------------------------------
        map_metadata = {
            "hx": {
                "output_name": "hx",
                "long_name": "Horizontal angular coordinate",
                "units": "arcsec",
            },
            "hy": {
                "output_name": "hy",
                "long_name": "Vertical angular coordinate",
                "units": "arcsec",
            },
            "lat": {
                "output_name": "latitude",
                "long_name": "Planetocentric latitude",
                "units": "degree",
            },
            "lon": {
                "output_name": "longitude",
                "long_name": "East-positive longitude",
                "units": "degree",
            },
            "emiss_ang": {
                "output_name": "emission_angle",
                "long_name": "Emission angle",
                "units": "degree",
            },
            "sol_ang": {
                "output_name": "solar_zenith_angle",
                "long_name": "Solar zenith angle",
                "units": "degree",
            },
            "azi_ang": {
                "output_name": "azimuth_angle",
                "long_name": "Azimuth angle",
                "units": "degree",
            },
            "lst": {
                "output_name": "local_solar_time",
                "long_name": "Local solar time",
                "units": "hour",
            },
        }
        
        for input_name, values in composite_geometry.items():
            values = np.asarray(values)
        
            if values.shape != (ny, nx):
                raise ValueError(
                    f"composite_geometry['{input_name}'] has shape "
                    f"{values.shape}; expected {(ny, nx)}"
                )
        
            metadata = map_metadata.get(
                input_name,
                {
                    "output_name": input_name,
                    "long_name": input_name,
                    "units": "unknown",
                },
            )
        
            ds = science.create_dataset(
                metadata["output_name"],
                data=values.astype(np.float32),
                compression="gzip",
                compression_opts=4,
                shuffle=True,
            )
        
            ds.attrs["long_name"] = metadata["long_name"]
            ds.attrs["units"] = metadata["units"]
            ds.attrs["dimensions"] = "y, x"
            ds.attrs["missing_value"] = np.nan

        
        # -------------------------------------------------------------
        # Scalar/global geometry
        # -------------------------------------------------------------
        ephemerides_group = h5.create_group("ephemerides")

        row = composite_ephemerides[0]  #Getting the table with only one row
        
        date_obs, time_obs = Time(
            row["datetime_jd"], format="jd", scale="utc"
        ).iso.split()

        lat_subobs = row["PDObsLat"]                              #Planetodetic sub-observer latitude (degrees)
        lon_subobs = (-row["PDObsLon"] + 180.0) % 360.0 - 180.0   #Planetodetic sub-observer longitude (degrees) - East positive
        lat_subsol = row["PDSunLat"]                              #Planetodetic sub-solar latitude (degrees)
        lon_subsol = (-row["PDSunLon"] + 180.0) % 360.0 - 180.0   #Planetodetic sub-solar longitude (degrees) - East positive

        #Changing from planetodetic to planetocentric
        target = str(row["targetname"]).split("(")[0].strip().upper()   #Name of planet
        print(target)
        equatorial_radius = archground.geometry.planet_dict[target]["equatorial_radius_km"]
        polar_radius = archground.geometry.planet_dict[target]["polar_radius_km"]
        lat_subobs = archground.geometry.planetodetic_to_planetocentric(
            lat_subobs,
            equatorial_radius=equatorial_radius,
            polar_radius=polar_radius,
        )
        lat_subsol = archground.geometry.planetodetic_to_planetocentric(
            lat_subsol,
            equatorial_radius=equatorial_radius,
            polar_radius=polar_radius,
        )
        
        ephemerides = {
            "date_obs": date_obs,
            "time_obs": time_obs,
            "planet": row["targetname"].split(" (")[0],
            "ra": float(row["RA"]),
            "dec": float(row["DEC"]),
            "pa": float(row["NPole_ang"]),
            "elevation": float(row["EL"]),
            "angsize": float(row["ang_width"]),
            "lat_subsol": float(lat_subsol),
            "lon_subsol": float(lon_subsol),
            "lat_subobs": float(lat_subobs),
            "lon_subobs": float(lon_subobs),
            "phase": float(row["alpha"]),
            "v_doppler_target_earth": float(row["delta_rate"]),
            "v_doppler_target_sun": float(row["r_rate"]),
            "dist_target_sun": float(row["r"]),
            "dist_target_earth": float(row["delta"]),
            "Ls": float(row["App_Lon_Sun"]),
        }
        
        ephemerides_units = {
            "date_obs": "ISO 8601 calendar date",
            "time_obs": "UTC",
            "planet": "1",
            "ra": "degree",
            "dec": "degree",
            "pa": "degree",
            "elevation": "degree",
            "angsize": "arcsec",
            "lat_subsol": "degree",
            "lon_subsol": "degree",
            "lat_subobs": "degree",
            "lon_subobs": "degree",
            "phase": "degree",
            "v_doppler_target_earth": "km s-1",
            "v_doppler_target_sun": "km s-1",
            "dist_target_sun": "AU",
            "dist_target_earth": "AU",
            "Ls": "degrees",
        }

        for name, value in ephemerides.items():
            if value is None:
                continue

            if isinstance(value, str):
                dtype = h5py.string_dtype(encoding="utf-8")
                ds = ephemerides_group.create_dataset(
                    name,
                    data=value,
                    dtype=dtype,
                )
            else:
                ds = ephemerides_group.create_dataset(
                    name,
                    data=value,
                )

            if name in ephemerides_units:
                ds.attrs["units"] = ephemerides_units[name]

        # Useful combined timestamp
        if "date_obs" in ephemerides and "time_obs" in ephemerides:
            timestamp = (
                f"{ephemerides['date_obs']}T"
                f"{ephemerides['time_obs']}"
            )

            dtype = h5py.string_dtype(encoding="utf-8")
            ds = ephemerides_group.create_dataset(
                "datetime_obs",
                data=timestamp,
                dtype=dtype,
            )
            ds.attrs["long_name"] = (
                "Representative composite observation time"
            )
            ds.attrs["timescale"] = "UTC"

        # Document coordinate conventions explicitly.
        ephemerides_group.attrs["latitude_type"] = "planetocentric"
        ephemerides_group.attrs["longitude_direction"] = "east_positive"
        ephemerides_group.attrs["longitude_range"] = "[-180, 180)"
        ephemerides_group.attrs["position_angle_convention"] = (
            "zero at north, positive toward east"
        )
    
        # -------------------------------------------------------------
        # Provenance
        # -------------------------------------------------------------
        if input_files is not None:
            provenance = h5.create_group("provenance")

            string_dtype = h5py.string_dtype(
                encoding="utf-8"
            )

            provenance.create_dataset(
                "input_files",
                data=np.asarray(input_files, dtype=object),
                dtype=string_dtype,
            )

    return output_filename

#######################################################################################

def read_composite_hdf5(filename):
    """
    Read the science and geometry sections of a composite HDF5 file.

    Returns
    -------
    result : dict
        Dictionary containing:

        result["science"]["data"]
        result["ephemerides"]["data"]
        result["science"]["dataset_attributes"]
        result["ephemerides"]["dataset_attributes"]
        result["file_attributes"]
    """

    with h5py.File(filename, "r") as h5:

        if "science" not in h5:
            raise KeyError(
                "The HDF5 file does not contain a 'science' group."
            )

        if "ephemerides" not in h5:
            raise KeyError(
                "The HDF5 file does not contain a 'ephemerides' group."
            )

        result = {
            "science": _read_group(h5["science"]),
            "ephemerides": _read_group(h5["ephemerides"]),
            "file_attributes": _read_attributes(h5),
        }

    return result





##############################################################################################
##############################################################################################
#                                           UTILS
##############################################################################################
##############################################################################################
    
def extract_values_telescope(value_fits):
    '''
    Change the format of N and E co-ordinates to numerical values (North is positive y, East is negative x).
    Supports values in decimal degrees, optionally followed by a quote (") symbol.
    '''
    
    import re

    # Regex now allows optional trailing double-quote (") after numbers
    pattern = r'([-+]?\d+(?:\.\d+)?)"?\s*([ENWS])'

    matches = re.findall(pattern, value_fits)

    x, y = None, None  # Default to None to distinguish missing values

    for value, direction in matches:
        value = float(value)  # Convert extracted string to float

        if direction == 'E':
            x = -value  # East is negative
        elif direction == 'W':
            x = value  # West is positive
        elif direction == 'N':
            y = value  # North is positive
        elif direction == 'S':
            y = -value  # South is negative

    # Return 0 if x or y is missing
    return (x if x is not None else 0, y if y is not None else 0)