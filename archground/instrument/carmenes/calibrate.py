#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# archGROUND - Python package for ground-based planetary spectroscopy and radiative transfer with archNEMESIS. 
# calibrate.py - Set of functions to process and calibrate CARMENES data.
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
from astropy.io import fits

#########################################################################################################

def process_fits(filename,print_file_structure=False):
    """
    Read the CARMENES FITS file, extract the data, and compute the relevant ephemerides.

    Parameters
    ----------
    filename : str
        Path to the CARMENES FITS file.

    Returns
    -------
    
    wavelength : ndarray (nwave, norders)
        Wavelength in microns.
    wavenumber : ndarray (nwave, norders)
        Wavenumber in cm-1.
    flux : ndarray (nwave, norders)
        Measured spectra in flux units.
    continuum : ndarray (nwave, norders)
        Estimated continuum in flux units.
    uncertainty : ndarray (nwave, norders)
        Uncertainty in flux units.
    ephemerides : dict
        Dictionary containing the computed ephemerides.
    """
    
    # ------------------------------------------------------------
    # Open and inspect the FITS file
    # ------------------------------------------------------------

    with fits.open(filename, memmap=True) as hdul:
        
        primary_header = hdul[0].header

        if print_file_structure is True:
            print("\nFITS structure:")
            hdul.info()
            print("\nPrimary-header information:")
            print("Object:    ", primary_header.get("OBJECT"))
            print("Date:      ", primary_header.get("DATE-OBS"))
            print("Instrument:", primary_header.get("INSTRUME"))
            print("Subsystem: ", primary_header.get("SUBSYS"))

        datetime = primary_header.get("DATE-OBS")
        planet = primary_header.get("OBJECT")

        # Copy data so it remains usable after the FITS file is closed.
        flux = hdul["SPEC"].data.astype(float)
        continuum = hdul["CONT"].data.astype(float)
        uncertainty = hdul["SIG"].data.astype(float)
        wavelength = hdul["WAVE"].data.astype(float)  #Angstroms

    #Converting wavelength from Angstroms to microns and calculating wavenumber in cm-1
    wavelength /= 1.0e4  #microns
    wavenumber = 1. / wavelength * 1.0e4  #cm-1

    #Transposing the data arrays to have the shape (nwave,norders)
    flux = flux.T
    continuum = continuum.T
    uncertainty = uncertainty.T
    wavelength = wavelength.T
    wavenumber = wavenumber.T

    #Calculating the ephemerides for the observation
    date = datetime.split("T")[0]
    time = datetime.split("T")[1]
    ephemerides = archground.get_ephemerides_epoch(date, time, planet, "CAHA")

    return wavelength, wavenumber, flux, continuum, uncertainty, ephemerides

#########################################################################################################