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
from archground.gcm.fvcd import vcd
import numpy as np
import sys,os

vcddir = archground.paths.vcddir

########################################################################################################################

def get_sun(lat,lon,Loct,vcddir=vcddir,euv_scenario=1,albedo_scenario=1):
    '''
    Function to get the Sun-Venus distance and the solar zenith angle

    Inputs
    ------

    lat :: Latitude
    lon :: Longitude
    Loct :: Local time at longitude lon

    Optional inputs
    ----------------

    vcddir :: directory where the VCD dataset is stored
    
    euv_scenario :: EUV scenario to use
                    (1) Standard solar EUV conditions
                    (2) Minimum solar EUV conditions
                    (3) Maximum solar EUV conditions
                    (4) Solar EUV from a given Julian date
                    (5) Solar EUV from a given input (varE107)
                    
    albedo_scenario :: Cloud albedo scenario
                       (1) Standard cloud albedo
                       (2) Low cloud albedo
                       (3) High cloud albedo
    
    Outputs
    --------

    dist_sun :: Sun-Venus distance (AU)
    sza :: Solar zenith angle (deg)

    '''

    #VCD inputs
    ##########################################

    z_key = 3      #(0 - Pressure in Pa ; 1 - Distance to planet center in m ; 2 - Altitude above reference sphere in m ; 3 - Altitude above surface in m)
    hires_key = 1     #High-resolution (1 - yes)
    date_key = 1   #0 - Earth time in Julian days ; 1 - Venus time localtime is the value of the local true solar time at longitude lon, given in martian hours.
    perturb_key = 0 # 0 - no perturbations ; 1 - add small scale perturbations (gravity waves)
    perturb_seed = 0 # perturbation seed (unused if perturb_key=0)
    perturb_gw_length = 0 # Gravity Wave length for perturbations (unused if perturb_key=0)


    if date_key==0:
        localtime=0.0
        juliandate=Loct
    if date_key==1:
        juliandate=0.0
        localtime=Loct
        
    if euv_scenario!=5:
        varE107=0.0
    
    z=0.
    (zon_wind,mer_wind,vert_wind,temp,pres,dens,extvar,seed_out,ier) = vcd(z_key,z,lon,lat,hires_key, \
                                                                            date_key,juliandate,localtime,vcddir,euv_scenario,albedo_scenario,varE107, \
                                                                            perturb_key,perturb_seed,perturb_gw_length,np.ones(100))
    
    dist_sun = extvar[7-1]
    sza = extvar[10-1]
    
    return dist_sun,sza


########################################################################################################################

def get_tsurf(lat,lon,Loct,vcddir=vcddir,euv_scenario=1,albedo_scenario=1):
    '''
    Function to get the surface temperature

    Inputs
    ------

    lat :: Latitude
    lon :: Longitude
    Loct :: Local time at longitude lon

    Optional inputs
    ----------------

    vcddir :: directory where the VCD dataset is stored
    
    euv_scenario :: EUV scenario to use
                    (1) Standard solar EUV conditions
                    (2) Minimum solar EUV conditions
                    (3) Maximum solar EUV conditions
                    (4) Solar EUV from a given Julian date
                    (5) Solar EUV from a given input (varE107)
                    
    albedo_scenario :: Cloud albedo scenario
                       (1) Standard cloud albedo
                       (2) Low cloud albedo
                       (3) High cloud albedo
    
    Outputs
    --------

    tsurf :: Surface temperature (K)

    '''

    #VCD inputs
    ##########################################

    z_key = 3      #(0 - Pressure in Pa ; 1 - Distance to planet center in m ; 2 - Altitude above reference sphere in m ; 3 - Altitude above surface in m)
    hires_key = 1     #High-resolution (1 - yes)
    date_key = 1   #0 - Earth time in Julian days ; 1 - Venus time localtime is the value of the local true solar time at longitude lon, given in martian hours.
    perturb_key = 0 # 0 - no perturbations ; 1 - add small scale perturbations (gravity waves)
    perturb_seed = 0 # perturbation seed (unused if perturb_key=0)
    perturb_gw_length = 0 # Gravity Wave length for perturbations (unused if perturb_key=0)


    if date_key==0:
        localtime=0.0
        juliandate=Loct
    if date_key==1:
        juliandate=0.0
        localtime=Loct
        
    if euv_scenario!=5:
        varE107=0.0
    
    z=0.
    (zon_wind,mer_wind,vert_wind,temp,pres,dens,extvar,seed_out,ier) = vcd(z_key,z,lon,lat,hires_key, \
                                                                            date_key,juliandate,localtime,vcddir,euv_scenario,albedo_scenario,varE107, \
                                                                            perturb_key,perturb_seed,perturb_gw_length,np.ones(100))
    
    tsurf = extvar[17-1]
    
    return tsurf



########################################################################################################################

def get_pt_profile(h,lat,lon,Loct,vcddir=vcddir,euv_scenario=1,albedo_scenario=1):
    '''
    Function to get the temperature profile

    Inputs
    ------

    h(nh) :: Altitude above the Venusian surface (km)
    lat :: Latitude
    lon :: Longitude
    Loct :: Local time at longitude lon

    Optional inputs
    ----------------

    vcddir :: directory where the VCD dataset is stored
    
    euv_scenario :: EUV scenario to use
                    (1) Standard solar EUV conditions
                    (2) Minimum solar EUV conditions
                    (3) Maximum solar EUV conditions
                    (4) Solar EUV from a given Julian date
                    (5) Solar EUV from a given input (varE107)
                    
    albedo_scenario :: Cloud albedo scenario
                       (1) Standard cloud albedo
                       (2) Low cloud albedo
                       (3) High cloud albedo
    
    Outputs
    --------

    press(nh) :: Pressure / Pa
    temp(nh) :: Temperature / K

    '''

    #VCD inputs
    ##########################################

    z_key = 3      #(0 - Pressure in Pa ; 1 - Distance to planet center in m ; 2 - Altitude above reference sphere in m ; 3 - Altitude above surface in m)
    hires_key = 1     #High-resolution (1 - yes)
    date_key = 1   #0 - Earth time in Julian days ; 1 - Venus time localtime is the value of the local true solar time at longitude lon, given in martian hours.
    perturb_key = 0 # 0 - no perturbations ; 1 - add small scale perturbations (gravity waves)
    perturb_seed = 0 # perturbation seed (unused if perturb_key=0)
    perturb_gw_length = 0 # Gravity Wave length for perturbations (unused if perturb_key=0)


    if date_key==0:
        localtime=0.0
        juliandate=Loct
    if date_key==1:
        juliandate=0.0
        localtime=Loct
        
    if euv_scenario!=5:
        varE107=0.0
    
    nh = len(h)
    press = np.zeros(nh)
    temp = np.zeros(nh)
    for i in range(nh):
        (zon_wind,mer_wind,vert_wind,temp[i],press[i],dens,extvar,seed_out,ier) = vcd(z_key,h[i]*1.0e3,lon,lat,hires_key, \
                                                                            date_key,juliandate,localtime,vcddir,euv_scenario,albedo_scenario,varE107, \
                                                                            perturb_key,perturb_seed,perturb_gw_length,np.ones(100))
    
    return press,temp


########################################################################################################################

def get_vmr_profile(h,lat,lon,Loct,vcddir=vcddir,euv_scenario=1,albedo_scenario=1,gasID=[2,5,7,45,48,39,1,79,9,78,19,3,15,22,40]):
    '''
    Function to get the temperature profile

    Inputs
    ------

    h(nh) :: Altitude above the Venusian surface (km)
    lat :: Latitude
    lon :: Longitude
    Loct :: Local time at longitude lon

    Optional inputs
    ----------------

    vcddir :: directory where the VCD dataset is stored
    
    euv_scenario :: EUV scenario to use
                    (1) Standard solar EUV conditions
                    (2) Minimum solar EUV conditions
                    (3) Maximum solar EUV conditions
                    (4) Solar EUV from a given Julian date
                    (5) Solar EUV from a given input (varE107)
                    
    albedo_scenario :: Cloud albedo scenario
                       (1) Standard cloud albedo
                       (2) Low cloud albedo
                       (3) High cloud albedo
                       
    gasID(ngas) :: Radtran ID of the gases to be extracted
    
    Outputs
    --------

    vmr(ngas,npro) :: Volume mixing ratio of each gas at each altitude

    '''

    #VCD inputs
    ##########################################

    z_key = 3      #(0 - Pressure in Pa ; 1 - Distance to planet center in m ; 2 - Altitude above reference sphere in m ; 3 - Altitude above surface in m)
    hires_key = 1     #High-resolution (1 - yes)
    date_key = 1   #0 - Earth time in Julian days ; 1 - Venus time localtime is the value of the local true solar time at longitude lon, given in martian hours.
    perturb_key = 0 # 0 - no perturbations ; 1 - add small scale perturbations (gravity waves)
    perturb_seed = 0 # perturbation seed (unused if perturb_key=0)
    perturb_gw_length = 0 # Gravity Wave length for perturbations (unused if perturb_key=0)


    if date_key==0:
        localtime=0.0
        juliandate=Loct
    if date_key==1:
        juliandate=0.0
        localtime=Loct
        
    if euv_scenario!=5:
        varE107=0.0
    
    nh = len(h)
    ngas = len(gasID)
    vmr = np.zeros((nh,ngas))
    
    gasID = np.array(gasID)
    
    for i in range(nh):
        (zon_wind,mer_wind,vert_wind,temp,press,dens,extvar,seed_out,ier) = vcd(z_key,h[i]*1.0e3,lon,lat,hires_key, \
                                                                            date_key,juliandate,localtime,vcddir,euv_scenario,albedo_scenario,varE107, \
                                                                            perturb_key,perturb_seed,perturb_gw_length,np.ones(100))
    
        for j in range(ngas):
             
            if gasID[j]==2: #CO2
                vmr[i,j] = extvar[50-1]
            elif gasID[j]==5: #CO
                vmr[i,j] = extvar[51-1]
            elif gasID[j]==7: #O2
                vmr[i,j] = extvar[52-1]
            elif gasID[j]==45: #O
                vmr[i,j] = extvar[53-1]
            elif gasID[j]==48: #H
                vmr[i,j] = extvar[54-1]
            elif gasID[j]==39: #H2
                vmr[i,j] = extvar[55-1]                
            elif gasID[j]==1: #H2O
                vmr[i,j] = extvar[56-1]
            elif gasID[j]==79: #H2SO4
                vmr[i,j] = extvar[57-1]
            elif gasID[j]==9: #SO2
                vmr[i,j] = extvar[58-1]
            elif gasID[j]==78: #SO
                vmr[i,j] = extvar[59-1]
            elif gasID[j]==19: #OCS
                vmr[i,j] = extvar[60-1]
            elif gasID[j]==3: #O3
                vmr[i,j] = extvar[61-1]
            elif gasID[j]==15: #HCl
                vmr[i,j] = extvar[62-1]
            elif gasID[j]==22: #N2
                vmr[i,j] = extvar[63-1]
            elif gasID[j]==40: #He
                vmr[i,j] = extvar[64-1]
            else:
                print('warning :: gasID',gasID[j],'not included in Venus Climate Database. Filling with zeros.')

    gasID = np.array(gasID)
    return gasID,vmr


########################################################################################################################

def convert_earth_venus_time(year,month,day,hour,minute,second,vcddir=vcddir,longitude=0.):
    '''
    Function to convert Earth time (UTC) to Venus time (Ls, sol and local time) 

    Inputs
    ------

    year, month, day :: Date (UTC)
    hour, minute, second :: Time (UTC)

    Optional inputs
    ----------------

    vcddir :: directory where the VCD dataset is stored
    longitude :: If 0, it provides the Universal Solar Time, otherwise it will provide Local Solar Time (LST)
    
    Outputs
    --------

    Ls :: Solar longitude 
    time :: UST or LST
    '''

    #Calculating Julian Date
    (ier,xdate)=fvcd.julian(month,day,year,hour,minute,second)

    #Calculating Ls out of Julian date
    #Ls,dist_sun = vcd.orbit(date)


    #VCD inputs
    ##########################################

    z_key = 3      #(0 - Pressure in Pa ; 1 - Distance to planet center in m ; 2 - Altitude above reference sphere in m ; 3 - Altitude above surface in m)
    hires_key = 1     #High-resolution (1 - yes)
    date_key = 0   #0 - Earth time in Julian days ; 1 - Venus time localtime is the value of the local true solar time at longitude lon, given in martian hours.
    perturb_key = 0 # 0 - no perturbations ; 1 - add small scale perturbations (gravity waves)
    perturb_seed = 0 # perturbation seed (unused if perturb_key=0)
    perturb_gw_length = 0 # Gravity Wave length for perturbations (unused if perturb_key=0)

    if date_key==0:
        localtime=0.0
        juliandate=xdate
        
    euv_scenario=1
    albedo_scenario=1
        
    if euv_scenario!=5:
        varE107=0.0
    
    z = 0.
    latitude = 0.
    (zon_wind,mer_wind,vert_wind,temp,press,dens,extvar,seed_out,ier) = vcd(z_key,0.,longitude,latitude,hires_key, \
                                                                        date_key,juliandate,localtime,vcddir,euv_scenario,albedo_scenario,varE107, \
                                                                        perturb_key,perturb_seed,perturb_gw_length,np.ones(100))

    Ls = extvar[6-1]
    LST = extvar[8-1]

    return Ls, LST

#####################################################################################################

def get_haus_cloud_model(h,latitude=0.,offset=0.):
    '''
    Function to compute the vertical structure of the Venus cloud using the model proposed in Haus et al. (2016)

    Inputs
    ------

    h(nh) :: Altitude (km)
    
    Optional inputs
    ----------------
    
    offset :: Altitude offset to apply to all layers simultaneously (km)
    
    Outputs
    --------

    n1(nh) :: Particle number density of cloud mode 1 (cm-3)
    n2(nh) :: Particle number density of cloud mode 2 (cm-3)
    n2p(nh) :: Particle number density of cloud mode 2' (cm-3)
    n3(nh) :: Particle number density of cloud mode 3 (cm-3)

    '''
    
    #Cheking that latitude is within boundaries
    latitude = abs(float(latitude))
    if not 0.0 <= latitude <= 90.0:
        raise ValueError(
            "LATITUDE must lie between -90 and 90 degrees"
        )
    

    #Get Table 4 from Haus+2016 to get latitude variations
    lat_mf = np.array(
        [
            0.0, 15.0, 20.0, 25.0,
            30.0, 35.0, 40.0, 45.0,
            50.0, 55.0, 60.0, 65.0,
            70.0, 75.0, 80.0, 90.0,
        ]
    )

    MF12_table = np.array(
        [
            0.98, 0.98, 0.99, 1.00,
            0.98, 0.94, 0.86, 0.81,
            0.73, 0.67, 0.64, 0.61,
            0.59, 0.47, 0.36, 0.36,
        ]
    )

    MF3_table = np.array(
        [
            1.30, 1.30, 1.26, 1.23,
            1.17, 1.13, 1.06, 1.03,
            1.04, 1.09, 1.22, 1.51,
            1.82, 2.02, 2.09, 2.09,
        ]
    )

    MF12 = np.interp(
        latitude,
        lat_mf,
        MF12_table,
    )

    MF3 = np.interp(
        latitude,
        lat_mf,
        MF3_table,
    )


    #Get Table 3 from Haus+2016 for cloud mode 2 latitude variation
    lat_mode2 = np.array(
        [
            0.0, 45.0, 50.0, 55.0,
            60.0, 65.0, 70.0, 75.0,
            80.0, 90.0,
        ]
    )

    zb2_table = np.array(
        [
            65.0, 65.0, 65.0, 65.0,
            64.5, 63.8, 63.1, 62.5,
            62.0, 62.0,
        ]
    )

    Hup2_table = np.array(
        [
            3.5, 3.5, 3.4, 3.2,
            2.6, 2.0, 1.0, 0.6,
            0.5, 0.5,
        ]
    )

    nh = len(h)

    #Cloud mode 1
    ###################################################
    
    zb1 = 49. + offset   #Lower base of peak altitude (km)
    zc1 = 16.    #Layer thickness of constant peak particle (km)
    Hup1 = 3.5   #Upper scale height (km)
    Hlo1 = 1.    #Lower scale height (km)
    #n01 = 193.5  #Particle number density at zb (cm-3)
    n01 = 193.5 * MF12
    
    #N1 = 3982.04e5 #Total column particle density (cm-2)
    #tau1 = 3.88    #Total column optical depth at 1 um
    
    n1 = np.zeros(nh)
    
    ialt1 = np.where(h<zb1)
    ialt2 = np.where((h<=(zb1+zc1)) & (h>=zb1))
    ialt3 = np.where(h>(zb1+zc1))
    
    n1[ialt1] = n01 * np.exp( -(zb1-h[ialt1])/Hlo1 )
    n1[ialt2] = n01
    n1[ialt3] = n01 * np.exp( -(h[ialt3]-(zb1+zc1))/Hup1 )
    
    #Cloud mode 2
    ###################################################
    
    zb2 = np.interp(
        latitude,
        lat_mode2,
        zb2_table,
    ) + offset    #Lower base of peak altitude (km)

    zc2 = 1.0   #Layer thickness of constant peak particle (km)

    Hup2 = np.interp(
        latitude,
        lat_mode2,
        Hup2_table,
    )         #Upper scale height (km)

    Hlo2 = 3.   #Lower scale height (km)
    #n02 = 100.  #Particle number density at zb (cm-3)
    n02 = 100.0 * MF12

    #N2 = 748.54e5  #Total column particle density (cm-2)
    #tau2 = 7.62    #Total column optical depth at 1 um
    
    n2 = np.zeros(nh)
    
    ialt1 = np.where(h<zb2)
    ialt2 = np.where((h<=(zb2+zc2)) & (h>=zb2))
    ialt3 = np.where(h>(zb2+zc2))
    
    n2[ialt1] = n02 * np.exp( -(zb2-h[ialt1])/Hlo2 )
    n2[ialt2] = n02
    n2[ialt3] = n02 * np.exp( -(h[ialt3]-(zb2+zc2))/Hup2 )
    
    #Cloud mode 2'
    ###################################################
    
    zb2p = 49. + offset   #Lower base of peak altitude (km)
    zc2p = 11.   #Layer thickness of constant peak particle (km)
    Hup2p = 1.0  #Upper scale height (km)
    Hlo2p = 0.1  #Lower scale height (km)
    n02p = 50.   #Particle number density at zb (cm-3)
    
    #N2p = 613.71e5  #Total column particle density (cm-2)
    #tau2p = 9.35    #Total column optical depth at 1 um
    
    n2p = np.zeros(nh)
    
    ialt1 = np.where(h<zb2p)
    ialt2 = np.where((h<=(zb2p+zc2p)) & (h>=zb2p))
    ialt3 = np.where(h>(zb2p+zc2p))
    
    n2p[ialt1] = n02p * np.exp( -(zb2p-h[ialt1])/Hlo2p )
    n2p[ialt2] = n02p
    n2p[ialt3] = n02p * np.exp( -(h[ialt3]-(zb2p+zc2p))/Hup2p )
    
    #Cloud mode 3
    ###################################################
    
    zb3 = 49. + offset  #Lower base of peak altitude (km)
    zc3 = 8.    #Layer thickness of constant peak particle (km)
    Hup3 = 1.0  #Upper scale height (km)
    Hlo3 = 0.5  #Lower scale height (km)
    #n03 = 14.   #Particle number density at zb (cm-3)
    n03 = 14.0 * MF3
    
    #N3 = 133.86e5   #Total column particle density (cm-2)  
    #tau3 = 14.14    #Total column optical depth at 1 um
    
    n3 = np.zeros(nh)
    
    ialt1 = np.where(h<zb3)
    ialt2 = np.where((h<=(zb3+zc3)) & (h>=zb3))
    ialt3 = np.where(h>(zb3+zc3))
    
    n3[ialt1] = n03 * np.exp( -(zb3-h[ialt1])/Hlo3 )
    n3[ialt2] = n03
    n3[ialt3] = n03 * np.exp( -(h[ialt3]-(zb3+zc3))/Hup3 )
    
    return n1,n2,n2p,n3



