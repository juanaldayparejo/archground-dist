#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# archGROUND - Python package for ground-based planetary spectroscopy and radiative transfer with archNEMESIS. 
# mars.py - Set of functions to compute archnemesis models for Mars
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

mcddir = archground.paths.mcddir

###########################################################################################################################

def create_atmosphere_mcd(h,lat,lon,Ls,LST,include_dust=True):
    """
    FUNCTION NAME : create_atmosphere_mcd()

    DESCRIPTION : Create the Atmosphere class from the MCD

    INPUTS : 

        h(nh) :: Altitude (km)
        lat :: Latitude (degrees)
        lon :: Longitude (degrees)
        Ls :: Solar longitude (degrees)
        LST :: Local time (h)

    OPTIONAL INPUTS:
    
        include_dust :: If True, it includes the dust profile from the MCD. The dust in the scatter class must be normalised to 1 at 0.67 um.
            
    OUTPUTS : 
 
        Atmosphere :: An instance of the `Atmosphere_0` class containing the reference atmosphere for Mars

    CALLING SEQUENCE:

        Atmosphere = create_atmosphere_mcd(lat,lon,Ls,LST)

    MODIFICATION HISTORY : Juan Alday (23/08/2026)
    """

    #Extracting p-T profile
    gasID = np.array([1,2,3,5,7,22,39,45,48,76]) #RADTRAN ID for each gas
    press,temp,rho,vmr = archground.gcm.extract_mcd.get_profile(h,gasID,lat,lon,Ls,LST,scenario=1,zkey=2,hrkey=1)

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
    Atmosphere.NVMR = len(gasID)
    Atmosphere.ID = np.array(gasID,dtype='int')
    Atmosphere.ISO = np.zeros(gasID.shape,dtype='int')
    Atmosphere.edit_VMR(vmr)
    Atmosphere.AMFORM = 0
    Atmosphere.calc_molwt()

    Atmosphere.NDUST = 0
    if include_dust is True:
        Atmosphere.NDUST += 1
        
        dust,rdust = archground.gcm.extract_mcd.get_dust_profile(h,lat,lon,Ls,LST,scenario=1,zkey=2,hrkey=1)
        #dust is in kg/kg_air, so we convert it to dust density in kg/m3
        dust *= rho
        Atmosphere.edit_DUST(dust[:,None])

        #normalising the profile so that it integrates to 1
        Atmosphere.normalise_dust(idust=0)

        #multiplying the dust profile by the visible column optical depth
        taudust = archground.gcm.extract_mcd.get_taudust(lat,lon,Ls,LST,scenario=1,hrkey=1)
        Atmosphere.DUST[:,0] *= taudust
    
    Atmosphere.assess()

    return Atmosphere

###########################################################################################################################

def create_scatter_class(waven,iscat=0,nmu=5,nf=10,nphi=100,include_waterice=False):
    """
    FUNCTION NAME : create_scatter_cass()

    DESCRIPTION : Create an instance of the `Scatter_0` class containing the cloud model for Mars

    INPUTS : 

        waven :: Wavenumber array of the measurement in cm-1

    OPTIONAL INPUTS:
    
        iscat :: Flag indicating the type of scattering calculations
                    iscat = 0 : No scattering
                    iscat = 1 : Multiple scattering
            
    OUTPUTS : 
 
        Scatter :: An instance of the `Scatter_0` class containing the cloud model for the Mars case

    CALLING SEQUENCE:

        Scatter = create_cloud_model(Atmosphere)

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    waven_min = waven.min() ; waven_max = waven.max()

    #Defining aerosol modes
    #Mode 1 - Dust
    reff1 = 1.0
    veff1 = 0.1
    r_g1 = reff1/(1.+veff1)**(5./2.)
    sigma_g1 = np.sqrt(np.log(1.0+veff1))

    #Mode 2 - Water ice
    reff2 = 3.0
    veff2 = 0.1
    r_g2 = reff2/(1.+veff2)**(5./2.)
    sigma_g2 = np.sqrt(np.log(1.0+veff2))

    #First of all we need to define our wavelength or wavenumber array
    #and tell the class in what units we want these calculations (wavenumber in cm-1 (ISPACE=0) or wavelength in um (ISPACE=1))
    Scatter = ans.Scatter_0()
    Scatter.ISCAT = iscat
    Scatter.NF = nf
    Scatter.NMU = nmu
    Scatter.NPHI = nphi
    Scatter.ISPACE = 0 #Wavenumber

    #Calculating maximum values for the Doppler shift
    c = 299792458.0   #Speed of light (m/s)
    v_doppler_max = 50.
    waven_min /= (1.0+v_doppler_max*1.0e3 / c)
    waven_max /= (1.0-v_doppler_max*1.0e3 / c)

    Scatter.IMIE = 2 #Legendre polynomial expansion of the phase function
    waven = np.arange(int(waven_min)-1.,waven_max+2.,1.)
    NDUST = 1      #Number of aerosol populations that we want to include in our atmosphere
    if include_waterice is True:
        NDUST += 1
    NWAVE = len(waven)    #Number of spectral points
    NTHETA = 361    #Number of phase angles for defining the phase function
    theta = np.linspace(0.,180.,NTHETA)

    #Now we initialise the arrays that will be filled with the calculations
    Scatter.initialise_arrays(NDUST,NWAVE,NTHETA,NLPOL=150)
    Scatter.WAVE = waven
    Scatter.THETA = theta

    #Mode 1
    print('Calculating Mode 1 - Dust')
    Scatter.read_refind(1)  #Reading optical properties of the Mars dust
    iscat = 2  #Log-normal distribution
    pars = np.array([r_g1,sigma_g1])

    idust = 0    #The index of the aerosol populations in the class that this calculation corresponds to (from 0 to NDUST-1)
    Scatter.makephase(idust,iscat,pars)

    if include_waterice is True:
        #Mode 2
        print('Calculating Mode 2 - Water ice')
        Scatter.read_refind(2)  #Reading optical properties of the water ice
        iscat = 2  #Log-normal distribution
        pars = np.array([r_g2,sigma_g2])

        idust = 1    #The index of the aerosol populations in the class that this calculation corresponds to (from 0 to NDUST-1)
        Scatter.makephase(idust,iscat,pars)


    #The dust profiles are normalised to the optical depth in the visible. Here we therefore need to normalise the KEXT to 0.67 um
    Scatter_NORM = ans.Scatter_0()
    Scatter_NORM.ISCAT = iscat
    Scatter_NORM.NF = nf
    Scatter_NORM.NMU = nmu
    Scatter_NORM.NPHI = nphi
    Scatter_NORM.ISPACE = 1 #Wavelength

    Scatter_NORM.IMIE = 2 #Legendre polynomial expansion of the phase function
    wavel = np.array([0.67])
    NDUST = 1      #Number of aerosol populations that we want to include in our atmosphere
    if include_waterice is True:
        NDUST += 1
    NWAVE = len(wavel)    #Number of spectral points
    NTHETA = 361    #Number of phase angles for defining the phase function
    theta = np.linspace(0.,180.,NTHETA)

    #Now we initialise the arrays that will be filled with the calculations
    Scatter_NORM.initialise_arrays(NDUST,NWAVE,NTHETA,NLPOL=150)
    Scatter_NORM.WAVE = wavel
    Scatter_NORM.THETA = theta

    #Mode 1
    Scatter_NORM.read_refind(1)  #Reading optical properties of the Mars dust
    iscat = 2  #Log-normal distribution
    pars = np.array([r_g1,sigma_g1])

    idust = 0    #The index of the aerosol populations in the class that this calculation corresponds to (from 0 to NDUST-1)
    Scatter_NORM.makephase(idust,iscat,pars)

    if include_waterice is True:
        #Mode 2
        Scatter_NORM.read_refind(2)  #Reading optical properties of the water ice
        iscat = 2  #Log-normal distribution
        pars = np.array([r_g2,sigma_g2])

        idust = 1    #The index of the aerosol populations in the class that this calculation corresponds to (from 0 to NDUST-1)
        Scatter_NORM.makephase(idust,iscat,pars)

    Scatter.KEXT /= Scatter_NORM.KEXT

    return Scatter

###########################################################################################################################

def create_surface_mcd(lat,lon,Ls,LST,galb=None):
    """
    FUNCTION NAME : create_surface_mcd()

    DESCRIPTION : Create the Surface class from the MCD

    INPUTS : 

        lat :: Latitude (degrees)
        lon :: Longitude (degrees)
        Ls :: Solar longitude (degrees)
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
    tsurf = archground.gcm.extract_mcd.get_tsurf(lat,lon,Ls,LST,mcddir=mcddir,scenario=1,hrkey=1)

    #Creating surface class
    Surface = ans.Surface_0()
    Surface.NLOCATIONS = 1
    Surface.LATITUDE = lat
    Surface.LONGITUDE = lon
    Surface.TSURF = tsurf
    Surface.LOWBC = 1  #Lambertian surface
    Surface.ISPACE = 0 #Wavenumber in cm-1

    #Albedo model (should be more accurate that this really)
    Surface.NEM = len(albedo_model)
    wavenumber_albedo_model = 1. / wavelength_albedo_model * 1.0e4
    isort = np.argsort(wavenumber_albedo_model)
    Surface.VEM = wavenumber_albedo_model[isort]
    Surface.EMISSIVITY = 1. - albedo_model[isort]
    Surface.assess()

    return Surface

###########################################################################################################################

def write_apr_prof(filename,Atmosphere,Measurement,varID,var2ID,error,clen=1.5,ScaleFactor=True,xprof=None):

    """

    FUNCTION NAME : write_apr_prof()

    DESCRIPTION : Write the file for retrieving a continuous profile of gas VMR, temperature or dust density

    INPUTS : 

        filename :: Name of the file (must be the same as indicated in .apr file)
        Atmosphere :: Python class defining the atmosphere
        varID :: Variable to be retrieved
                  - varID = 0 indicates that temperature is to be retrieved
                  - varID > 0 indicates that a gas vmr profile is to be retrieved. varID
                              must then be equal to the gasID of the gas that wants to be retrieved
                  - varID < 0 indicates that an aerosol density profile is to be retrieved.
                              varID must then indicate the aerosol population number that
                              wants to be retrieved (e.g. -2 indicates that the second 
                              population of aerosols in the aerosol.ref wants to be retrieved)

        var2ID :: Isotopologue to be retrieved. If varID = 0 or varID < 0 then this variables,
                  must be set to 0. If a gas vmr profile is to be retrieved, then var2ID must
                  be equal to the isotopologue ID of the gas that wants to be retrieved. 
        err :: A priori uncertainty on the profile to be retrieved, expressed as an absolute
               magnitude in the same units as the profile is defined. 

    OPTIONAL INPUTS:
    
        clen :: Correlation length, expressed as a factor of scale height
        ScaleFactor :: If True, then error represents a fractional error rather than absolute
        xprof(npro) :: If we do not want to use the profile in the Atmosphere as the initial file
                        we can define xprof to use it as the a priori profile (default is None)
            
    OUTPUTS : 
 
        A priori file

    CALLING SEQUENCE:

        write_apr_prof(filename,Atmosphere,gasID,isoID)

    MODIFICATION HISTORY : Juan Alday (06/05/2023)

    """
    
    npro = Atmosphere.NP
    
    #Temperature
    if varID==0:
        xref = np.zeros(Atmosphere.NP)
        if xprof is None:
            xref[:] = Atmosphere.T
        else:
            if len(xprof)!=npro:
                raise ValueError('error :: xprof must have the same length as NPRO')
            else:
                xref[:] = xprof
    
    #Dust    
    elif varID<0:
        caero = abs(varID)
        if caero>Atmosphere.NDUST:
            raise ValueError('error :: The aerosol population that wants to be retrieved does not exist in aerosol.ref')
        
        xref = np.zeros(npro)
        if xprof is None:
            xref[0:npro] = Atmosphere.DUST[0:npro,caero-1]
        else:
            if len(xprof)!=npro:
                raise ValueError('error :: xprof must have the same length as NPRO')
            else:
                xref[:] = xprof
                
    #Gas abundance      
    elif varID>0:
        cgas1 = np.where((Atmosphere.ID==varID) & (Atmosphere.ISO==var2ID))
        cgas = cgas1[0]

        xref = np.zeros(npro)
        if xprof is None:
            xref[0:npro] = Atmosphere.VMR[0:npro,cgas[0]]
        else:
            if len(xprof)!=npro:
                raise ValueError('error :: xprof must have the same length as NPRO')
            else:
                xref[:] = xprof
                
    #Creating error array
    errarr = np.zeros(npro)
    errarr[0:npro] = error

    if ScaleFactor==True:
        errarr = xref*errarr

    #Write file
    fref = open(filename,'w')
    fref.write('\t %i \t %10.3f \n' % (npro,clen))
    for i in range(npro):
        fref.write('\t %10.6e \t %10.6e \t %10.6e \n' % (Atmosphere.P[i],xref[i],errarr[i]))
    fref.close()


###############################################################################################

def create_apr_file(runname,
        Atmosphere, 
        Measurement,
        retrieve_continuous_prof=False,cont_varID1=None,cont_varID2=None,cont_clen=None,cont_err=None,
        retrieve_temp=False,flag_temp_analytic=True,temp_clen=1.5,temp_err=None,
        retrieve_press=False,htan=40.,ptan=None,ptanerr=0.1,
        retrieve_scaling_factor=False,scale_varID1=None,scale_varID2=None,scale_apr=None,scale_err=None,
        retrieve_scaling_factor_telluric=False,tel_varID1=None,tel_varID2=None,tel_apr=None,tel_err=None,
        retrieve_baseline=False,baseline_degree=2,
        ):

    """
    FUNCTION NAME : create_apr_file()

    DESCRIPTION : Function to create the archNEMESIS .apr file to define the retrieved parameters

    INPUTS : 

        runname :: Name of the archNEMESIS run
        Atmosphere :: Atmosphere class
        Measurement :: Measurement class

    OPTIONAL INPUTS:
    
        retrieve_temp :: Flag to retrieve temperature profile with numerical calculation of Jacobian (i.e., hydrostatic approach)
            flag_temp_analytic :: Flag indicating whether the temperature retrieval approach must be from hydrostatic (False) or spectroscopy (True)
            temp_clen :: Correlation length expressed as a fraction of the scale height (required if retrieve_temp=True)
            temp :: A priori uncertainty in the temperature (K) (required if retrieve_temp=True)
        retrieve_press :: Flag to retrieve the pressure at a given tangent height
            htan :: Tangent height at which the pressure must be retrieved in km (required if retrieve_press=True)
            ptan :: Pressure at the given tangent height (optional if retrieve_press=True. If None then it will be taken from the Atmosphere class)
            ptanerr :: Fractional error in the a priori pressure (required if retrieve_press=True)
        retrieve_continuous_prof :: Flag to retrieve a continuous profile of gas VMR or dust number density
            cont_varID1(ngas) :: ID of the gases or aerosols to be retrieved (see archnemesis documentation. required if retrieve_continuous_prof=True)
            cont_varID2(ngas) :: ID of the isotopes to be retrieved (see archnemesis documentation. required if retrieve_continuous_prof=True)
            cont_clen :: Correlation length expressed as a fraction of the scale height (required if retrieve_continuous_prof=True)
            cont_err(ngas) :: Fractional error in the gas vmrs or dust densities (required if retrieve_continuous_prof=True)
        retrieve_scaling_factor :: Flag to retrieve a scaling factor of gas VMR or dust number density
            scale_varID1(ngas) :: ID of the gases or aerosols to be retrieved (see archnemesis documentation. required if retrieve_scaling_factor=True)
            scale_varID2(ngas) :: ID of the isotopes to be retrieved (see archnemesis documentation. required if retrieve_scaling_factor=True)
            scale_err(ngas) :: Fractional error in the gas vmrs or dust densities (required if retrieve_scaling_factor=True)
        retrieve_scaling_factor_telluric :: Flag to retrieve a scaling factor of gas VMR or dust number density in telluric atmosphere
            tel_varID1(ngas) :: ID of the gases or aerosols to be retrieved (see archnemesis documentation. required if retrieve_scaling_factor_telluric=True)
            tel_varID2(ngas) :: ID of the isotopes to be retrieved (see archnemesis documentation. required if retrieve_scaling_factor_telluric=True)
            tel_err(ngas) :: Fractional error in the gas vmrs or dust densities (required if retrieve_scaling_factor_telluric=True)
        retrieve_baseline :: Flag to retrieve the baseline at each tangent height with a polynomial function
            baseline_degree :: Degree of the polynomial that must be used to fit the baseline (required if retrieve_baseline=True)

    OUTPUTS : 
 
        archNEMESIS .apr file

    MODIFICATION HISTORY : Juan Alday (29/07/2026)

    """

    #Counting number of retrieved parameters
    ####################################################################################################

    nvar = 0

    #Counting number of variables from continuous profiles
    if retrieve_continuous_prof is True:
        nvar += len(cont_varID1)

    if retrieve_scaling_factor is True:
        nvar += len(scale_varID1)

    #Counting other variables
    nvar += int(retrieve_temp) + int(retrieve_press) + int(retrieve_baseline)

    #Writing the file
    ####################################################################################################

    fapr = open(runname+'.apr','w')
    fapr.write('#Ground-based retrieval \n')
    fapr.write('\t'+str(nvar)+' \n')

    #Going through the continuous profiles
    if retrieve_continuous_prof is True:

        ngas = len(cont_varID1)
        if len(cont_varID1) != len(cont_varID2):
            raise ValueError("error while writing continuous profiles in .apr file :: varID1 and varID2 must be of the same length")

        if cont_err is None:
            raise ValueError("error while writing continuous profiles in .apr file :: cont_err must be defined for each of the retrieved species")

        if cont_clen is None:
            cont_clen = np.ones(ngas) * 1.5 

        for igas in range(ngas):

            gasid = cont_varID1[igas] ; isoid = cont_varID2[igas]

            fapr.write('\t %i \t %i \t %i \n' % (gasid,isoid,0))
            filename = 'contprof'+str(igas)+'.dat'
            fapr.write(filename+' \n')

            write_apr_prof(filename,Atmosphere,gasid,isoid,cont_err[igas],clen=cont_clen[igas],ScaleFactor=True)


    #Going through the scaling factors
    if retrieve_scaling_factor is True:

        ngas = len(scale_varID1)
        if len(scale_varID1) != len(scale_varID2):
            raise ValueError("error while writing scaling factors in .apr file :: varID1 and varID2 must be of the same length")

        if scale_err is None or scale_apr is None:
            raise ValueError("error while writing scaling factors in .apr file :: scale_apr and scale_err must be defined for each of the retrieved species")

        for igas in range(ngas):

            gasid = scale_varID1[igas] ; isoid = scale_varID2[igas]

            fapr.write('\t %i \t %i \t %i \n' % (gasid,isoid,3))
            fapr.write('\t %7.4e \t %7.4e \n' % (scale_apr[igas],scale_err[igas]))


    #Going through the scaling factors in telluric atmosphere
    if retrieve_scaling_factor_telluric is True:

        ngas = len(tel_varID1)
        if len(tel_varID1) != len(tel_varID2):
            raise ValueError("error while writing scaling factors in .apr file :: varID1 and varID2 must be of the same length")

        if tel_err is None or tel_apr is None:
            raise ValueError("error while writing scaling factors in .apr file :: tel_apr and tel_err must be defined for each of the retrieved species")

        for igas in range(ngas):

            gasid = tel_varID1[igas] ; isoid = tel_varID2[igas]

            fapr.write('\t %i \t %i \t %i \n' % (gasid,isoid,103))
            fapr.write('\t %7.4e \t %7.4e \n' % (tel_apr[igas],tel_err[igas]))



    #Going through the retrieval of the baseline
    if retrieve_baseline is True:
        
        fapr.write('\t %i \t %i \t %i \n' % (231,0,231))
        fapr.write('transapr.dat \n')

        #Getting an estimate number for the continuum
        ngeom = Measurement.NGEOM
        tau = np.zeros(ngeom)
        for it in range(ngeom):
            tau[it] = np.mean(Measurement.MEAS[0:Measurement.NCONV[it],it])

        #Defining number of windows, measurements and polynomial degree
        ftau = open('transapr.dat','w')
        ftau.write('\t '+str(ngeom)+' \t '+str(baseline_degree)+' \n')

        if baseline_degree == 0:

            data = np.zeros((2,ngeom))
            data[0,:] = tau
            data[1,:] = tau * 0.5
            for it in range(ngeom):
                ftau.write('\t %10.5e \t %10.5e \n' % (data[0,it],data[1,it]))

        elif baseline_degree == 1:

            data = np.zeros((4,ngeom))
            for it in range(ngeom):
                
                data[0,it] = tau[it]          #A0
                data[1,it] = tau[it]*0.5      #A0_ERR
                data[2,it] = 0.001 / (Measurement.VCONV[Measurement.NCONV[0]-1,0]-Measurement.VCONV[0,0])         #A1
                data[3,it] = 0.05 / (Measurement.VCONV[Measurement.NCONV[0]-1,0]-Measurement.VCONV[0,0])        #A1_ERR
                
                ftau.write('\t %10.5e \t %10.5e \t %10.5e \t %10.5e \n' % (data[0,it],data[1,it],data[2,it],data[3,it]))

        elif baseline_degree==2:

            data = np.zeros((6,ngeom))
            for it in range(ngeom):
                data[0,it] = tau[it]          #A0
                data[1,it] = tau[it]*0.5      #A0_ERR
                data[2,it] = 0.01 / (Measurement.VCONV[Measurement.NCONV[0]-1,0]-Measurement.VCONV[0,0])         #A1
                data[3,it] = 0.05 / (Measurement.VCONV[Measurement.NCONV[0]-1,0]-Measurement.VCONV[0,0])        #A1_ERR
                data[4,it] = -0.01 / (Measurement.VCONV[Measurement.NCONV[0]-1,0]-Measurement.VCONV[0,0])**2.    #A2
                data[5,it] = 0.05 / (Measurement.VCONV[Measurement.NCONV[0]-1,0]-Measurement.VCONV[0,0])**2.    #A2_ERR
                
                ftau.write('\t %10.5e \t %10.5e \t %10.5e \t %10.5e \t %10.5e \t %10.5e \n' % (data[0,it],data[1,it],data[2,it],data[3,it],data[4,it],data[5,it]))
        
        ftau.close()


    #Going through the temperature retrieval
    if retrieve_temp is True:

        if flag_temp_analytic is True:
            fapr.write('\t %i \t %i \t %i \n' % (0,0,0))
        else:
            fapr.write('\t %i \t %i \t %i \n' % (0,-1,0))
        fapr.write('tempapr.dat \n')

        filename = 'tempapr.dat'
        write_apr_prof(filename,Atmosphere,0,0,temp_err,clen=temp_clen,ScaleFactor=False)

    
    #Going through the pressure retrieval
    if retrieve_press is True:

        if htan is None:
            raise ValueError("error while writing the pressure retrieval in the .apr file :: htan must be defined if we want to retrieve pressure")


        fapr.write('\t %i \t %i \t %i \n' % (666,0,666))
        ipress0 = np.argmin(np.abs(Atmosphere.H/1.0e3-htan))
        tanhe0 = Atmosphere.H[ipress0]/1.0e3

        fapr.write('\t %10.7f \n' % (tanhe0))

        if ptan is None:
            ptan = Atmosphere.P[ipress0]/101325.

        perr = ptan * ptanerr
        fapr.write('\t %10.7e \t %10.7e \n' % (ptan,perr))


    fapr.close()



#ALBEDO MODEL

wavelength_albedo_model = np.array([0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600, 1.700, 1.800, 1.900, 2.000, 2.100, 2.200, 2.300, 2.400, 2.500, 2.600, 2.700, 2.800,
 2.900, 3.000, 3.100, 3.200, 3.300, 3.400, 3.500, 3.600, 3.700, 3.800, 3.900, 4.000, 4.100, 4.200, 4.300, 4.400, 4.500, 4.600, 4.700, 4.800, 4.900, 5.000, 5.100, 5.200, 5.300, 5.400, 5.500, 5.600,
 5.700, 5.800, 5.900, 6.000, 6.100, 6.200, 6.300, 6.400, 6.500, 6.600, 6.700, 6.800, 6.900, 7.000, 7.100, 7.200, 7.300, 7.400, 7.500, 7.600, 7.700, 7.800, 7.900, 8.000, 8.100, 8.200, 8.300, 8.400,
 8.500, 8.600, 8.700, 8.800, 8.900, 9.000, 9.100, 9.200, 9.300, 9.400, 9.500, 9.600, 9.700, 9.800, 9.900, 10.000, 10.100, 10.200, 10.300, 10.400, 10.500, 10.600, 10.700, 10.800, 10.900, 11.000,
 11.100, 11.200, 11.300, 11.400, 11.500, 11.600, 11.700, 11.800, 11.900, 12.000, 12.100, 12.200, 12.300, 12.400, 12.500, 12.600, 12.700, 12.800, 12.900, 13.000, 13.100, 13.200, 13.300, 13.400, 13.500,
 13.600, 13.700, 13.800, 13.900, 14.000, 14.100, 14.200, 14.300, 14.400, 14.500, 14.600, 14.700, 14.800, 14.900, 15.000, 15.100, 15.200, 15.300, 15.400, 15.500, 15.600, 15.700, 15.800, 15.900, 16.000,
 16.100, 16.200, 16.300, 16.400, 16.500, 16.600, 16.700, 16.800, 16.900, 17.000, 17.100, 17.200, 17.300, 17.400, 17.500, 17.600, 17.700, 17.800, 17.900, 18.000, 18.100, 18.200, 18.300, 18.400, 18.500,
 18.600, 18.700, 18.800, 18.900, 19.000, 19.100, 19.200, 19.300, 19.400, 19.500, 19.600, 19.700, 19.800, 19.900, 20.000, 20.100, 20.200, 20.300, 20.400, 20.500, 20.600, 20.700, 20.800, 20.900, 21.000,
 21.100, 21.200, 21.300, 21.400, 21.500, 21.600, 21.700, 21.800, 21.900, 22.000, 22.100, 22.200, 22.300, 22.400, 22.500, 22.600, 22.700, 22.800, 22.900, 23.000, 23.100, 23.200, 23.300, 23.400, 23.500,
 23.600, 23.700, 23.800, 23.900, 24.000, 24.100, 24.200, 24.300, 24.400, 24.500, 24.600, 24.700, 24.800, 24.900, 25.000, 25.100, 25.200, 25.300, 25.400, 25.500, 25.600, 25.700, 25.800, 25.900, 26.000,
 26.100, 26.200, 26.300, 26.400, 26.500, 26.600, 26.700, 26.800, 26.900, 27.000, 27.100, 27.200, 27.300, 27.400, 27.500, 27.600, 27.700, 27.800, 27.900, 28.000, 28.100, 28.200, 28.300, 28.400, 28.500,
 28.600, 28.700, 28.800, 28.900, 29.000, 29.100, 29.200, 29.300, 29.400, 29.500, 29.600, 29.700, 29.800, 29.900, 30.000])

albedo_model = np.array([0.018, 0.019, 0.027, 0.063, 0.125, 0.209, 0.241, 0.240, 0.220, 0.220, 0.220, 0.220, 0.222, 0.225, 0.228, 0.228, 0.228, 0.228, 0.228, 0.228, 0.227, 0.226, 0.225, 0.224, 0.223, 0.220, 0.216, 0.211,
 0.207, 0.203, 0.199, 0.195, 0.191, 0.187, 0.183, 0.178, 0.174, 0.170, 0.167, 0.163, 0.160, 0.157, 0.153, 0.150, 0.146, 0.143, 0.139, 0.136, 0.132, 0.129, 0.125, 0.120, 0.116, 0.112, 0.107, 0.103,
 0.099, 0.094, 0.090, 0.085, 0.081, 0.077, 0.072, 0.068, 0.063, 0.059, 0.055, 0.050, 0.046, 0.041, 0.037, 0.032, 0.028, 0.023, 0.019, 0.014, 0.010, 0.010, 0.010, 0.010, 0.013, 0.017, 0.020, 0.023,
 0.027, 0.029, 0.031, 0.033, 0.035, 0.037, 0.037, 0.037, 0.037, 0.037, 0.037, 0.037, 0.036, 0.036, 0.036, 0.035, 0.035, 0.034, 0.034, 0.033, 0.033, 0.032, 0.032, 0.031, 0.031, 0.030, 0.029, 0.028,
 0.028, 0.027, 0.026, 0.025, 0.024, 0.023, 0.022, 0.021, 0.021, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.023, 0.023, 0.023, 0.023, 0.023, 0.023, 0.023, 0.024, 0.024, 0.024, 0.024, 0.024,
 0.024, 0.024, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026, 0.026, 0.027, 0.027, 0.027, 0.027, 0.027, 0.027, 0.027, 0.028, 0.028, 0.028, 0.028,
 0.028, 0.028, 0.028, 0.029, 0.029, 0.029, 0.029, 0.029, 0.029, 0.029, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.031, 0.031, 0.031, 0.031, 0.031, 0.031, 0.031, 0.032, 0.032, 0.032, 0.032,
 0.032, 0.032, 0.032, 0.033, 0.033, 0.033, 0.033, 0.033, 0.033, 0.033, 0.033, 0.033, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.035, 0.035, 0.035, 0.035, 0.035, 0.035,
 0.035, 0.035, 0.035, 0.035, 0.036, 0.036, 0.036, 0.036, 0.036, 0.036, 0.036, 0.036, 0.036, 0.036, 0.037, 0.037, 0.037, 0.037, 0.037, 0.037, 0.037, 0.037, 0.037, 0.037, 0.038, 0.038, 0.038, 0.038,
 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.041, 0.041, 0.041,
 0.041, 0.041, 0.041, 0.041, 0.041, 0.041, 0.041, 0.042, 0.042, 0.042, 0.042, 0.042, 0.042, 0.042, 0.042, 0.042, 0.042, 0.043, 0.043, 0.043])

