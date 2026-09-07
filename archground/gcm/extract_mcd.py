#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# archGROUND - Python package for ground-based planetary spectroscopy and radiative transfer with archNEMESIS. 
# extract_mcd.py - Set of functions to extract data from the Mars Climate Database.
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
from archground.gcm.fmcd import mcd
import numpy as np
import sys,os

mcddir = archground.paths.mcddir


#MARS CLIMATE DATABASE
#################################################################################

# MCD scenarios:
# 1=Climatology Scenario, solarEUVaverageconditions 
# 2=Climatology Scenario, solarEUVminimumconditions 
# 3=Climatology Scenario, solarEUVmaximumconditions 
# 4=dust storm τ=5, solarminimumconditions 
# 5=dust storm τ=5, solaraveragedconditions 
# 6=dust storm τ=5, solarmaximumconditions 
# 7=warm scenario: dusty atmosphere, solarmax 
# 8=cold scenario: low-dust conditions, solarmin. 
# 24=Mars Year 24, with associated solar EUV conditions. 
# 25=Mars Year 25, with associated solar EUV conditions. 
# 26=Mars Year 26, with associated solar EUV conditions. 
# 27=Mars Year 27, with associated solar EUV conditions. 
# 28=Mars Year 28, with associated solar EUV conditions. 
# 29=Mars Year 29, with associated solar EUV conditions. 
# 30=Mars Year 30, with associated solar EUV conditions. 
# 31=Mars Year 31, with associated solar EUV conditions. 
# 32=Mars Year 32, with associated solar EUV conditions. 
# 33=Mars Year 33, with associated solar EUV conditions. 
# 34=Mars Year 34, with associated solar EUV conditions. 
# 35=Mars Year 35, with associated solar EUV conditions.


########################################################################################################################

def get_profile(h,gasID,lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,zkey=2,hrkey=1):
    '''
    Function to get the vertical profiles of the volume mixing ratios 

    Inputs
    ------

    h(nh) :: Altitude above the Martian areoid (km)
    gasID(ngas) :: Gas ID
    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    zkey :: Altitudes defined above Martian areoid (2) or surface (3)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    press(nh) :: Pressure (Pa)
    temp(nh) :: Temeprature (K)
    rho(nh) :: Specific density (kg/m3)
    vmr(nh,ngas) :: Volume mixing ratios of the different gases
    '''

    #MCD inputs
    ##########################################

    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    nh = len(h)
    ngas = len(gasID)


    vmr = np.zeros((nh,ngas))
    press = np.zeros(nh)
    temp = np.zeros(nh)
    rho = np.zeros(nh)
    for i in range(nh):
        (press[i], rho[i], temp[i], zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,h[i]*1000.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100,dtype='int32'))


        for igas in range(ngas):

            if gasID[igas]==1: #H2O
                vmr[i,igas] = extvar[48-1]
            elif gasID[igas]==2: #CO2
                vmr[i,igas] = extvar[64-1]
            elif gasID[igas]==22: #N2
                vmr[i,igas] = extvar[65-1]
            elif gasID[igas]==76: #Ar
                vmr[i,igas] = extvar[66-1]
            elif gasID[igas]==5: #CO
                vmr[i,igas] = extvar[67-1]
            elif gasID[igas]==45: #O
                vmr[i,igas] = extvar[68-1]
            elif gasID[igas]==7: #O2
                vmr[i,igas] = extvar[69-1]
            elif gasID[igas]==3: #O3
                vmr[i,igas] = extvar[70-1]
            elif gasID[igas]==48: #H
                vmr[i,igas] = extvar[71-1]
            elif gasID[igas]==39: #H2
                vmr[i,igas] = extvar[72-1]
            elif gasID[igas]==40: #He
                vmr[i,igas] = extvar[73-1]
            else:
                sys.exit('error in get_mcd_vmr_profile :: gasID is not included in the MCD')

    return press,temp,rho,vmr

########################################################################################################################

def get_dust_profile(h,lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,zkey=2,hrkey=1):
    '''
    Function to get the vertical profiles of the dust density and effective radius

    Inputs
    ------

    h(nh) :: Altitude above the Martian areoid (km)
    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    zkey :: Altitudes defined above Martian areoid (2) or surface (3)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    dust(nh) :: Dust density (kg/kg_air)
    rdust(nh) :: Dust effective radius (m)

    '''

    #MCD inputs
    ##########################################

    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    nh = len(h)

    press = np.zeros(nh)
    temp = np.zeros(nh)
    dust = np.zeros(nh)
    rdust = np.zeros(nh)
    for i in range(nh):
        (press[i], dens, temp[i], zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,h[i]*1000.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100))
        dust[i] = extvar[41-1]
        rdust[i] = extvar[42-1]
    
    return dust,rdust


########################################################################################################################

def get_waterice_profile(h,lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,zkey=2,hrkey=1):
    '''
    Function to get the vertical profiles of the water ice density and effective radius

    Inputs
    ------

    h(nh) :: Altitude above the Martian areoid (km)
    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    zkey :: Altitudes defined above Martian areoid (2) or surface (3)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    waterice(nh) :: Water ice density (mol/mol_air)
    rdust(nh) :: Dust effective radius (m)

    '''

    #MCD inputs
    ##########################################

    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    nh = len(h)

    press = np.zeros(nh)
    temp = np.zeros(nh)
    waterice = np.zeros(nh)
    rwaterice = np.zeros(nh)
    for i in range(nh):
        (press[i], dens, temp[i], zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,h[i]*1000.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100))
        waterice[i] = extvar[50-1]
        rwaterice[i] = extvar[51-1]
    
    return waterice,rwaterice

########################################################################################################################

def get_sun(lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,hrkey=1):
    '''
    Function to get the Sun-Mars distance and the solar zenith angle

    Inputs
    ------

    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    dist_sun :: Sun-Mars distance (AU)
    sza :: Solar zenith angle (deg)

    '''

    #MCD inputs
    ##########################################

    zkey = 3
    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    (press, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,0.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100))
    dist_sun = extvar[8-1]
    sza = extvar[13-1]
    
    return dist_sun,sza


########################################################################################################################

def get_tsurf(lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,hrkey=1):
    '''
    Function to get the surface temperature 

    Inputs
    ------

    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    tsurf :: Surface Temperature (K)
    '''

    #MCD inputs
    ##########################################

    zkey = 3
    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    (press, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,0.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100))
    
    tsurf = extvar[14-1]

    return tsurf

########################################################################################################################

def get_taudust(lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,hrkey=1):
    '''
    Function to get the Daily mean dust column visible optical depth above
    surface. From local surface to the top of the atmosphere, at wavelength
    0.67 µm.

    Inputs
    ------

    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    taudust :: Dust visible optical depth
    '''

    #MCD inputs
    ##########################################

    zkey = 3
    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    (press, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,0.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100))
    
    taudust = extvar[40-1]

    return taudust


########################################################################################################################

def get_icecol(lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,hrkey=1):
    '''
    Function to get the water ice column in kg/m2

    Inputs
    ------

    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    icecol :: Water ice column (kg m-2)
    '''

    #MCD inputs
    ##########################################

    zkey = 3
    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    (press, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,0.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100))
    
    icecol = extvar[49-1]

    return icecol


########################################################################################################################

def get_galb(lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,hrkey=1):
    '''
    Function to get the GCM surface bare ground albedo

    Inputs
    ------

    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    taudust :: Dust visible optical depth
    '''

    #MCD inputs
    ##########################################

    zkey = 3
    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    (press, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,0.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100))
    
    galb = extvar[38-1]

    return galb

########################################################################################################################

def get_o3col(lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,hrkey=1):
    '''
    Function to get the ozone column density 

    Inputs
    ------

    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    o3col :: Ozone column density (m-2)
    '''

    #MCD inputs
    ##########################################

    zkey = 3
    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    (press, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,0.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100))
    
    o3_mars = 0.048 #kg mol-1
    NA = 6.022e23 #molecules mol-1
    
    o3col = extvar[80-1] / o3_mars * NA   #m-2

    return o3col


########################################################################################################################

def get_ps(lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,hrkey=1):
    '''
    Function to get the surface pressure 

    Inputs
    ------

    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    ps :: Surface pressure (Pa)
    '''

    #MCD inputs
    ##########################################

    zkey = 3
    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    (press, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,0.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100))
    
    surpress = extvar[15-1] #Pa

    return surpress

########################################################################################################################

def get_waterice_surface(lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,hrkey=1):
    '''
    Function to get the surface pressure 

    Inputs
    ------

    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    waterice_nonperennial :: Non-perennial surface water ice (kg/m2)
    waterice_perennial :: Perennial surface water ice (1/0 - on/off)
    '''

    #MCD inputs
    ##########################################

    zkey = 3
    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    (press, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,0.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100))
    
    waterice_perennial = extvar[46-1]
    waterice_nonperennial = extvar[45-1]

    return waterice_nonperennial,waterice_perennial

########################################################################################################################

def get_vmr_profile(h,gasID,lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,zkey=2,hrkey=1):
    '''
    Function to get the vertical profiles of the volume mixing ratios 

    Inputs
    ------

    h(nh) :: Altitude above the Martian areoid (km)
    gasID(ngas) :: Gas ID
    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    zkey :: Altitudes defined above Martian areoid (2) or surface (3)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    vmr(ngas,nh) :: Volume mixing ratios of the different gases
    '''

    #MCD inputs
    ##########################################

    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    nh = len(h)
    ngas = len(gasID)


    vmr = np.zeros((ngas,nh))
    for i in range(nh):
        (pres, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,h[i]*1000.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100,dtype='int32'))


        for igas in range(ngas):

            if gasID[igas]==1: #H2O
                vmr[igas,i] = extvar[48-1]
            elif gasID[igas]==2: #CO2
                vmr[igas,i] = extvar[64-1]
            elif gasID[igas]==22: #N2
                vmr[igas,i] = extvar[65-1]
            elif gasID[igas]==76: #Ar
                vmr[igas,i] = extvar[66-1]
            elif gasID[igas]==5: #CO
                vmr[igas,i] = extvar[67-1]
            elif gasID[igas]==45: #O
                vmr[igas,i] = extvar[68-1]
            elif gasID[igas]==7: #O2
                vmr[igas,i] = extvar[69-1]
            elif gasID[igas]==3: #O3
                vmr[igas,i] = extvar[70-1]
            elif gasID[igas]==48: #H
                vmr[igas,i] = extvar[71-1]
            elif gasID[igas]==39: #H2
                vmr[igas,i] = extvar[72-1]
            elif gasID[igas]==40: #He
                vmr[igas,i] = extvar[73-1]
            else:
                sys.exit('error in get_mcd_vmr_profile :: gasID is not included in the MCD')

    return vmr


########################################################################################################################

def get_pt_profile(h,lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,zkey=2,hrkey=1):
    '''
    Function to get the vertical profiles of the pressure and temperature 

    Inputs
    ------

    h(nh) :: Altitude above the Martian areoid (km)
    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    zkey :: Altitudes defined above Martian areoid (2) or surface (3)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    press(nh) :: Pressure (Pa)
    temp(nh) :: Temperature (K)
    '''

    #MCD inputs
    ##########################################

    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    nh = len(h)

    press = np.zeros(nh)
    temp = np.zeros(nh)
    for i in range(nh):
        (press[i], dens, temp[i], zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,h[i]*1000.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.zeros(100))

    return press,temp


########################################################################################################################

def get_molwt_profile(h,lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,zkey=2,hrkey=1):
    '''
    Function to get the vertical profiles of the molecular weight (kg mol-1)

    Inputs
    ------

    h(nh) :: Altitude above the Martian areoid (km)
    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    zkey :: Altitudes defined above Martian areoid (2) or surface (3)
    hrkey :: High resolution flag (1) or not (0)
    
    Outputs
    --------

    molwt(nh) :: Molecular weight (kg mol-1)
    '''

    #MCD inputs
    ##########################################

    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    nh = len(h)

    #rho = np.zeros(nh)
    #press = np.zeros(nh)
    #temp = np.zeros(nh)
    k_B = 1.380649e-23
    molwt = np.zeros(nh)
    for i in range(nh):
        (press, rho, temp, zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,h[i]*1000.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.zeros(100))
        molwt[i] = rho * k_B * temp / press

    #numdens = press / (1.380649e-23*temp)

    return molwt


########################################################################################################################

def get_h2o_column(lat,lon,Ls,Loct,mcddir=mcddir,scenario=1):
    '''
    Function to get the vertical profiles of the molecular weight (kg mol-1)

    Inputs
    ------

    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    
    Outputs
    --------

    h2o_col :: Water vapour column (kg/m2)
                 To convert to pr .um just multiple h2o_col * 1000.
    '''
    
    #MCD inputs
    ##########################################

    zkey = 3
    hrkey = 1     #High-resolution (1 - yes)
    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    (press, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,0.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100))
    
    h2ocol = extvar[47-1]
    
    return h2ocol
    

########################################################################################################################

def get_dens_profile(h,lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,zkey=2):
    '''
    Function to get the vertical profiles of the density and number density

    Inputs
    ------

    h(nh) :: Altitude above the Martian areoid (km)
    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    zkey :: Altitudes defined above Martian areoid (2) or surface (3)
    
    Outputs
    --------

    rho(nh) :: Atmospheric density (kg m-3)
    numdens(nh) :: Number density (m-3)
    '''

    #MCD inputs
    ##########################################

    hrkey = 1     #High-resolution (1 - yes)
    datekey = 1   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

    nh = len(h)

    rho = np.zeros(nh)
    press = np.zeros(nh)
    temp = np.zeros(nh)
    k_B = 1.380649e-23
    for i in range(nh):
        (press[i], rho[i], temp[i], zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,h[i]*1000.,lon,lat,hrkey,datekey,Ls,Loct,mcddir,scenario,perturkey,seedin,gwlength,np.zeros(100))

    numdens = press / (1.380649e-23*temp)

    return rho,numdens


########################################################################################################################

def convert_earth_mars_time(year,month,day,hour,minute,second,mcddir=mcddir,longitude=0.):
    '''
    Function to convert Earth time (UTC) to Mars time (Ls, sol and local time) 

    Inputs
    ------

    year, month, day :: Date (UTC)
    hour, minute, second :: Time (UTC)

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    longitude :: If 0, it provides the Universal Solar Time, otherwise it will provide Local Solar Time (LST)
    
    Outputs
    --------

    Ls :: Solar longitude 
    time :: UST or LST
    '''

    #Calculating Julian Date
    (ier,xdate)=mcd.julian(month,day,year,hour,minute,second)

    #MCD inputs
    ##########################################

    zkey = 3
    hrkey = 1     #High-resolution (1 - yes)
    datekey = 0   #Dates in Ls
    perturkey = 1 # default to no perturbation
    seedin = 0 # perturbation seed (unused if perturkey=1)
    gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)
    lat = 0.
    loct = 0.
    scenario = 1

    (press, dens, temp, zonwind, merwind, meanvar, extvar, seedout, ierr)  = mcd(zkey,0.,longitude,lat,hrkey,datekey,xdate,loct,mcddir,scenario,perturkey,seedin,gwlength,np.ones(100))
    Ls = extvar[9-1]
    LST = extvar[10-1]

    return Ls, LST
    
########################################################################################################################

def get_homopause(lat,lon,Ls,Loct,mcddir=mcddir,scenario=1,MakePlot=False):
    '''
    Function to calculate the homopause altitude predicted from the Mars Climate Database.
    This altitude is calculated by computing the ratios of the Ar/N2 ratio.

    Inputs
    ------
    
    lat :: Latitude
    lon :: Longitude
    Ls :: Solar longitude
    Loct :: Local time

    Optional inputs
    ----------------

    mcddir :: directory where the MCD dataset is stored
    scenario :: MCD scenario to use (1 - Climatological)
    MakePlot :: If True, a summary plot is made
    
    Outputs
    --------

    h_homopause :: Altitude of the homopause (km)
    '''
    
    from scipy.interpolate import interp1d

    #Getting the vertical profiles for Ar and N2    
    h = np.linspace(20.,250.,231)
    gasID = [76,22]
    press,temp,rho,vmr = get_mcd_profile(h,gasID,lat,lon,Ls,Loct,mcddir=mcddir,scenario=scenario)
    
    #Calculating the Ar/N2 ratio
    y = np.zeros(len(h))
    y[:] = vmr[:,0]/vmr[:,1]
    
    #Calculating the lower atmospheric value
    ilo = np.where((h>=30.) & (h<=40.))[0]
    ylo = np.mean(y[ilo])
    
    #Calculating the height where the Ar/N2 ratio is 0.95 of that in the lower atmosphere
    s = interp1d(y,h)
    h_homopause = s(ylo*0.95)
    
    #Calculating the Ar/N2 ratio and its first and second derivatives
    #dy_dx = np.gradient(y, h)
    #dy2_dx2 = np.gradient(dy_dx, h)
    
    #Smoothing the second derivative (a value of 10 seems to work fine)
    #dy2_dx2s = smooth_array(dy2_dx2, 10)
    
    #Finding the value with minimum second derivative
    #h_homopause = h[np.argmin(dy2_dx2s)]
    
    if MakePlot==True:
        
        fig,ax1 = plt.subplots(1,1,figsize=(4,6))
        ax1.plot(y,h,c='black')
        ax1.axhline(h_homopause)
        ax1.set_xlabel('Ar/N$_2$ ratio')
        ax1.set_ylabel('Altitude (km)')
        ax1.grid()
        plt.tight_layout()
    
    return h_homopause
    

def smooth_array(input_array, window_size):
    """
    Smooths an array using a simple moving average.

    Parameters:
    - input_array: The input array to be smoothed.
    - window_size: The size of the moving average window.

    Returns:
    - smoothed_array: The smoothed array.
    """
    weights = np.ones(window_size) / window_size
    smoothed_array = np.convolve(input_array, weights, mode='same')
    return smoothed_array
