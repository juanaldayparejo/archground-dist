#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# archGROUND - Python package for ground-based planetary spectroscopy and radiative transfer with archNEMESIS. 
# spectral_windows.py - Dictionary for each spectral window to be used in the retrievals of CARMENES data
#
# Copyright (C) 2026 Juan Alday
#
# ACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

process_orders_info = {

    ###########################################################################################
    # POSITION 6 - DERIVATION OF 18O/16O in CO2 
    ###########################################################################################
    
    "alday_12co2_nir": {

        #Measurement
        "diffor" : 23,                     #Diffraction order (from 0 to 27)
        "waven_min" : 6519.6,              #Minimum wavenumber (cm-1)
        "waven_max" : 6540.0,              #Maximum wavenumber (cm-1)
        "normalise" : True,                #Flag to indicate if radiance must be normalised
        "vnorm" : 6521.75,                 #Wavenumber at which to normalise the spectra

        #Atmosphere
        "id_act" : [2],                    #Radtran ID of active gases in the atmosphere
        "iso_act" : [1],                   #Radtran isotope ID of active gases in the atmosphere
        "split_CO_iso": False,             #Flag indicating whether CO must be separated into its 4 main isotopes
        "split_H2O_iso": True,             #Flag indicating whether H2O must be separated into its 4 main isotopes
        "split_CO2_iso": True,             #Flag indicating whether CO2 must be separated into its 4 main isotopes

        #Telluric
        "id_act_tel" : [1,2],              #Radtran ID of active gases in the telluric atmosphere
        "iso_act_tel" : [0,0],             #Radtran isotope ID of active gases in the telluric atmosphere

        #Retrieval
        "scaling_factor_id" : [0,2],       #IDs of the atmospheric parameters to be retrieved through scaling factors
        "scaling_factor_iso" : [0,1],      #IDs of the atmospheric parameters to be retrieved through scaling factors
        "scaling_factor_apr" : [1.0,1.0],  #A priori scaling factors
        "scaling_factor_err" : [0.5,0.5],  #A priori uncertainty in scaling factors

        "tel_scaling_factor_id" : [1,2],       #IDs of the telluric atmospheric parameters to be retrieved through scaling factors
        "tel_scaling_factor_iso" : [0,0],      #IDs of the telluric atmospheric parameters to be retrieved through scaling factors
        "tel_scaling_factor_apr" : [1.0,1.0],  #A priori scaling factors
        "tel_scaling_factor_err" : [0.5,0.5],  #A priori uncertainty in scaling factors 

        "flag_baseline" : True,            #Flag indicating whether a baseline must be retrieved
    },

    "alday_13co2_nir": {

        #Measurement
        "diffor" : 22,                     #Diffraction order (from 0 to 27)
        "waven_min" : 6736.5,              #Minimum wavenumber (cm-1)
        "waven_max" : 6775.5,              #Maximum wavenumber (cm-1)
        "normalise" : True,                #Flag to indicate if radiance must be normalised
        "vnorm" : 6739.75,                 #Wavenumber at which to normalise the spectra

        #Atmosphere
        "id_act" : [2,1],                  #Radtran ID of active gases in the atmosphere
        "iso_act" : [2,1],                 #Radtran isotope ID of active gases in the atmosphere
        "split_CO_iso": False,             #Flag indicating whether CO must be separated into its 4 main isotopes
        "split_H2O_iso": True,             #Flag indicating whether H2O must be separated into its 4 main isotopes
        "split_CO2_iso": True,             #Flag indicating whether CO2 must be separated into its 4 main isotopes

        #Telluric
        "id_act_tel" : [1,2],              #Radtran ID of active gases in the telluric atmosphere
        "iso_act_tel" : [0,0],             #Radtran isotope ID of active gases in the telluric atmosphere

        #Retrieval
        "scaling_factor_id" : [1,2],       #IDs of the atmospheric parameters to be retrieved through scaling factors
        "scaling_factor_iso" : [1,2],      #IDs of the atmospheric parameters to be retrieved through scaling factors
        "scaling_factor_apr" : [1.0,1.0],  #A priori scaling factors
        "scaling_factor_err" : [0.5,0.5],  #A priori uncertainty in scaling factors

        "tel_scaling_factor_id" : [1,2],       #IDs of the telluric atmospheric parameters to be retrieved through scaling factors
        "tel_scaling_factor_iso" : [0,0],      #IDs of the telluric atmospheric parameters to be retrieved through scaling factors
        "tel_scaling_factor_apr" : [1.0,1.0],  #A priori scaling factors
        "tel_scaling_factor_err" : [0.5,0.5],  #A priori uncertainty in scaling factors 

        "flag_baseline" : True,            #Flag indicating whether a baseline must be retrieved
    },

}