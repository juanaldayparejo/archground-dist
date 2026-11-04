import numpy as np
from struct import *
import sys,os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import archground
import archnemesis as ans


##############################################################################################
##############################################################################################
#                                CREATE ARCHNEMESIS CLASSES
##############################################################################################
##############################################################################################


def create_retrieval_class(niter,philimit=0.1,ncores=1):
    """
    FUNCTION NAME : create_retrieval_class()

    DESCRIPTION : Create an instance of the `Retrieval_0` class containing the retrieval parameters for the Venus case

    INPUTS : 

        None

    OPTIONAL INPUTS:
    
        None
            
    OUTPUTS : 
 
        Retrieval :: An instance of the `Retrieval_0` class containing the retrieval parameters for the Venus case

    CALLING SEQUENCE:

        Retrieval = create_retrieval_class()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    Retrieval = ans.OptimalEstimation_0(IRET=0)

    Retrieval.NITER = niter       #Number of iterations
    Retrieval.PHILIMIT = philimit   #Convergence criterion
    Retrieval.NCORES = ncores      #Number of available cores

    Retrieval.assess_input()

    return Retrieval

###########################################################################################################################

def create_layer_class(Atmosphere,nlay=101):
    """
    FUNCTION NAME : create_layer_class()

    DESCRIPTION : Create an instance of the `Layer_0` class containing the layer parameters

    INPUTS : 

        Atmosphere :: Instance of the Atmosphere class

    OPTIONAL INPUTS:
    
        None
            
    OUTPUTS : 
 
        Layer :: An instance of the `Layer_0` class containing the layer parameters for the Venus case

    CALLING SEQUENCE:

        Layer = create_layer_class()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    Layer = ans.Layer_0(Atmosphere.RADIUS)

    Layer.NLAY = nlay        #Number of layers
    Layer.LAYTYP = 1         #Layering performed with equal changes in log pressure
    Layer.LAYINT = 1       

    mask = (Atmosphere.T > 0.0)
    hvalid = Atmosphere.H[mask]
    layht = hvalid[0]
    Layer.LAYHT = layht #m
    Layer.assess()

    return Layer

###########################################################################################################################

def create_measurement(wavenumber,radiance,error,lat,lon,emiss_ang,sol_ang,azi_ang,
                        lat_subobs=None,lon_subobs=None,v_doppler=0.0,resolving_power=90000.):
    """
    FUNCTION NAME : create_measurement()

    DESCRIPTION : Create a measurement instance with the given parameters

    INPUTS : 

        wavenumber :: Wavenumber grid (cm-1)
        radiance :: Radiance (erg s-1 cm-2 sr-1 (cm-1)-1)
        error :: Radiance error (erg s-1 cm-2 sr-1 (cm-1)-1)
        lat :: Latitude (degrees)
        lon :: Longitude (degrees)
        emiss_ang :: Emission angle (degrees)
        sol_ang :: Solar zenith angle (degrees)
        azi_ang :: Azimuth angle (degrees)

    OPTIONAL INPUTS:
    
        lat_subobs :: Latitude of the sub-observer point (degrees)
        lon_subobs :: Longitude of the sub-observer point (degrees)
            
    OUTPUTS : 
 
        Measurement :: An instance of the `Measurement` class containing the reference measurement from TEXES

    CALLING SEQUENCE:

        Measurement = create_measurement(wavenumber,radiance,error,lat,lon,emiss_ang,sol_ang,azi_ang)

    MODIFICATION HISTORY : Juan Alday (23/08/2026)
    """

    Measurement = ans.Measurement_0()
    Measurement.NGEOM = 1
    Measurement.NCONV = np.zeros(Measurement.NGEOM,dtype='int') + len(wavenumber)
    Measurement.edit_VCONV(wavenumber[:,None])
    Measurement.edit_MEAS(radiance[:,None] * 1.0e-7) #Converting from erg s-1 cm-2 sr-1 (cm-1)-1 to W cm-2 sr-1 (cm-1)-1
    Measurement.edit_ERRMEAS(error[:,None] * 1.0e-7) #Converting from erg s-1 cm-2 sr-1 (cm-1)-1 to W cm-2 sr-1 (cm-1)-1

    Measurement.NAV = np.ones(Measurement.NGEOM,dtype='int')
    Measurement.edit_FLAT(np.ones((Measurement.NGEOM,Measurement.NAV.max())) * lat )
    Measurement.edit_FLON(np.ones((Measurement.NGEOM,Measurement.NAV.max())) * lon )
    Measurement.edit_EMISS_ANG(np.ones((Measurement.NGEOM,Measurement.NAV.max())) * emiss_ang )
    Measurement.edit_SOL_ANG(np.ones((Measurement.NGEOM,Measurement.NAV.max())) * sol_ang )
    Measurement.edit_AZI_ANG(np.ones((Measurement.NGEOM,Measurement.NAV.max())) * azi_ang )
    Measurement.edit_WGEOM(np.ones((Measurement.NGEOM,Measurement.NAV.max())) )

    Measurement.V_DOPPLER = v_doppler

    if lat_subobs is not None:
        Measurement.SUBOBS_LAT = lat_subobs
    else:
        Measurement.SUBOBS_LAT = lat

    if lon_subobs is not None:
        Measurement.SUBOBS_LON = lon_subobs
    else:
        Measurement.SUBOBS_LON = lon

    #Building the ILS
    midwave = np.mean(Measurement.VCONV[:,0])
    fwhm = midwave / resolving_power
    Measurement.ISHAPE = 2 #Gaussian ILS
    Measurement.FWHM = fwhm
    Measurement.build_ils()

    Measurement.assess()

    return Measurement

#########################################################################################################

def create_telluric_online(observatory,waven,date_obs,time_obs,elevation,
                            delv=0.001,hitran_file=archground.paths.archnemesis_hitran24,tips_file=archground.paths.archnemesis_tips,
                            era5=True,
                            resolving_power=90000.,id_gas_act=[1,2,3,4,5,6,7]):
    """
    FUNCTION NAME : create_telluric_online()

    DESCRIPTION : Create an instance of the `Telluric_0` class containing the telluric atmosphere
                   and the spectroscopy calculated from the HITRAN24 and TIPS files

    INPUTS : 

        observatory :: Name of the observatory (see dictionary in geometry.py)
        waven :: Wavenumber array of the measurement (cm-1)
        date_obs :: Observation date
        time_obs :: Observation time (UTC)
        elevation :: Elevation of the target above the local horizon (degrees)

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

    #Extracting the observatory parameters
    lon_obs = archground.geometry.observatory_dict[observatory]["lon"]
    lat_obs = archground.geometry.observatory_dict[observatory]["lat"]
    alt_obs = archground.geometry.observatory_dict[observatory]["elevation"]

    #Defining the Telluric class
    Telluric = ans.Telluric_0()

    #Atmosphere
    ##################################################################################

    #Defining the inputs
    Telluric.DATE=date_obs               #UTC date of the observation
    Telluric.TIME=time_obs               #UTC time of the observation
    Telluric.LATITUDE=lat_obs            #Latitude of the observatory
    Telluric.LONGITUDE=lon_obs           #Longitude of the observatory
    Telluric.ALTITUDE=alt_obs * 1000.    #Altitude of the observatory
    Telluric.EMISS_ANG=180. - (90. - elevation)    #Observing angle, looking at the zenith is 180 degrees, and at the horizon is 90 degrees

    if era5 is True:
        #Extracting the atmosphere from the ERA5 model
        Telluric.extract_atmosphere_era5()
    else:
        Telluric.extract_atmosphere_circ()

    Telluric.Atmosphere.NDUST = 0

    #Spectroscopy
    ##################################################################################

    Telluric.Spectroscopy = ans.Spectroscopy_0(ILBL=1)

    Telluric.Spectroscopy.NGAS = len(id_gas_act)
    Telluric.Spectroscopy.ID = np.array(id_gas_act,dtype='int32')
    Telluric.Spectroscopy.ISO = np.zeros_like(id_gas_act,dtype="int32")
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

#########################################################################################################

def create_stellar_class(waven,dist_sun,v_sun_planet,t_sun=5772.0):
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

    #Reading Toon's linelist
    toon_file = os.path.join(archground.paths.archground_path,'../','data','solar','solar_merged_20240731_600_33300_100.out')
    data = np.loadtxt(toon_file, skiprows=3)
    waven_solar = data[:,0] 
    trans_solar = data[:,1]

    if( ( (waven.min()-5.) < waven_solar.min()) or ( (waven.max()+5.) > waven_solar.max() ) ):
        print("warning :: Toon's solar linelist outside of range. Modelling solar spectrum as blackbody without lines")

        waven_solarx = np.arange(waven.min() - 5. , waven.max() + 5., 0.1)
        trans_solarx = np.ones_like(waven_solarx)

    else:
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
    Stellar.DIST = dist_sun
    Stellar.calc_luminosity_blackbody(t_sun)

    Stellar.SOLSPEC *= trans_solarx

    #Applying Doppler shift between star and planet
    c_kms = 299792.458
    doppler_factor = 1.0 + v_sun_planet / c_kms
    Stellar.WAVE /= doppler_factor
    Stellar.assess()

    return Stellar

###########################################################################################################################

def create_spectroscopy_class_online(waven,id_gases,iso_gases,
                                    delv=0.001,
                                    hitran_file=archground.paths.archnemesis_hitran24,
                                    tips_file=archground.paths.archnemesis_tips,
                                    resolving_power=90000.):
    """
    FUNCTION NAME : create_spectroscopy_class()

    DESCRIPTION : Create an instance of the `Spectroscopy_0` class containing the spectroscopy parameters for a calculation of 
                   the absorption cross sections at runtime

    INPUTS : 

        waven :: Wavenumber array of the measurement (cm-1)
        delv :: Wavenumber step for the line-by-line opacity calculation in cm-1

    OPTIONAL INPUTS:
    
        id_gases :: Array containing the gas IDs of the gases to be included in the spectroscopy (default: None, which means that only the main gases in the Venus atmosphere will be included)
        iso_gases :: Array containing the isotopologue IDs of the gases to be included in the spectroscopy (default: None, which means that only the main isotopologues of the main gases
        hitran_file :: Path to the HITRAN file to be used for the line-by-line opacity calculation (default: texes.paths.archnemesis_hitran24, which is a custom HITRAN24 file containing the main gases in the Venus atmosphere)
        tips_file :: Path to the TIPS file to be used for the partition function calculation (default: texes.paths.archnemesis_tips, which is a custom TIPS file containing the partition function data for the main gases in the Venus atmosphere)
    
    OUTPUTS : 
 
        Spectroscopy :: An instance of the `Spectroscopy_0` class containing the spectroscopy parameters

    CALLING SEQUENCE:

        Spectroscopy = create_spectroscopy_class_online()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    Spectroscopy = ans.Spectroscopy_0()

    Spectroscopy.ILBL = 1  #runtime calculation of the line-by-line opacity

    if id_gases is not None:
        ids =  []
        iso = []
        for i in range(len(id_gases)):
            ids.append(id_gases[i])
            iso.append(iso_gases[i])
    else:
        raise ValueError("error :: no gases are specified when creating spectroscopy class")

    Spectroscopy.NGAS = len(ids)
    Spectroscopy.ID = ids
    Spectroscopy.ISO = iso
    Spectroscopy.ISPACE = 0
    Spectroscopy.IPROC = np.zeros(Spectroscopy.NGAS,dtype='int32')
    Spectroscopy.LOCATION_LD = [hitran_file] * Spectroscopy.NGAS
    Spectroscopy.LOCATION_PF = [tips_file] * Spectroscopy.NGAS
    Spectroscopy.LOCATION_CD = [hitran_file] * Spectroscopy.NGAS

    Spectroscopy.LINE_DATA_PARAMS = [ans.MolLineDataParams()] * Spectroscopy.NGAS

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

    Spectroscopy.NWAVE = nwave
    Spectroscopy.WAVE = wave

    return Spectroscopy


##############################################################################################
##############################################################################################
#                                  ATMOSPHERIC CLASS
##############################################################################################
##############################################################################################

##############################################################################################

def split_H2O_isotopes(Atmosphere,dhratio=5.,o18ratio=1.,o17ratio=1.):

    """
    FUNCTION NAME : split_water_isotopes()

    DESCRIPTION : Function to split H2O into its 4 main isotopes

    INPUTS : 

        Atmosphere :: archNEMESIS atmosphere class

    OPTIONAL INPUTS:
    
        dhratio :: Value of D/H (VSMOW)
        o18ratio :: Value of 18O/16O (VSMOW)
        o17ratio :: Value of 17O/16O (VSMOW)

    OUTPUTS : 
 
        Atmosphere :: archNEMESIS Atmosphere class

    CALLING SEQUENCE:

        Atmosphere = split_H2O_isotopes(Atmosphere,dhratio=5.,o18ratio=1.,o17ratio=1.)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Splitting H2O into 4 isotopes
    ih2o = np.where( (Atmosphere.ID==1) & (Atmosphere.ISO==0) )[0][0]

    vmr_h2o = Atmosphere.VMR[:,ih2o]

    Atmosphere.remove_gas(1,0)

    #Adding the isotopes of H2O
    Atmosphere.add_gas(1,1,vmr_h2o)
    Atmosphere.add_gas(1,2,vmr_h2o*2005.2e-6*o18ratio)      #VSMOW
    Atmosphere.add_gas(1,3,vmr_h2o*379.9e-6*o17ratio)       #VSMOW
    Atmosphere.add_gas(1,4,vmr_h2o*(155.76e-6*2.)*dhratio)  #VSMOW

    return Atmosphere

###############################################################################################

def split_CO2_isotopes(Atmosphere,c13ratio=1.,o18ratio=1.,o17ratio=1.):

    """
    FUNCTION NAME : split_water_isotopes()

    DESCRIPTION : Function to split H2O into its 4 main isotopes

    INPUTS : 

        Atmosphere :: archNEMESIS atmosphere class

    OPTIONAL INPUTS:
    
        c13ratio :: Value of 13C/12C (VPDB)
        o18ratio :: Value of 18O/16O (VSMOW)
        o17ratio :: Value of 17O/16O (VSMOW)

    OUTPUTS : 
 
        Atmosphere :: archNEMESIS Atmosphere class

    CALLING SEQUENCE:

        Atmosphere = split_CO2_isotopes(Atmosphere,c13ratio=1.,o18ratio=1.,o17ratio=1.)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Splitting CO2 into 4 isotopes
    ico2 = np.where( (Atmosphere.ID==2) & (Atmosphere.ISO==0) )[0][0]

    vmr_co2 = Atmosphere.VMR[:,ico2]

    Atmosphere.remove_gas(2,0)

    #Adding the isotopes of CO2 isotopes (following HITRAN fractionation)
    Atmosphere.add_gas(2,1,vmr_co2*0.984204)
    Atmosphere.add_gas(2,2,vmr_co2*0.011057*c13ratio)       #VPDB
    Atmosphere.add_gas(2,3,vmr_co2*0.003947*o18ratio)       #VSMOW
    Atmosphere.add_gas(2,4,vmr_co2*7.339890e-4*o17ratio)    #VSMOW

    return Atmosphere

###############################################################################################

def split_CO_isotopes(Atmosphere,c13ratio=1.,o18ratio=1.,o17ratio=1.):

    """
    FUNCTION NAME : split_CO_isotopes()

    DESCRIPTION : Function to split H2O into its 4 main isotopes

    INPUTS : 

        Atmosphere :: archNEMESIS atmosphere class

    OPTIONAL INPUTS:
    
        c13ratio :: Value of 13C/12C (VPDB)
        o18ratio :: Value of 18O/16O (VSMOW)
        o17ratio :: Value of 17O/16O (VSMOW)

    OUTPUTS : 
 
        Atmosphere :: archNEMESIS Atmosphere class

    CALLING SEQUENCE:

        Atmosphere = split_CO2_isotopes(Atmosphere,c13ratio=1.,o18ratio=1.,o17ratio=1.)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Splitting CO into 4 isotopes
    ico = np.where( (Atmosphere.ID==5) & (Atmosphere.ISO==0) )[0][0]

    vmr_co = Atmosphere.VMR[:,ico]

    Atmosphere.remove_gas(5,0)

    #Adding the isotopes of CO
    Atmosphere.add_gas(5,1,vmr_co)
    Atmosphere.add_gas(5,2,vmr_co*0.0112372*c13ratio)   #VPDB
    Atmosphere.add_gas(5,3,vmr_co*2005.2e-6*o18ratio)   #VSMOW
    Atmosphere.add_gas(5,4,vmr_co*379.9e-6*o17ratio)    #VSMOW

    return Atmosphere

###############################################################################################

def split_O2_isotopes(Atmosphere,o18ratio=1.,o17ratio=1.):

    """
    FUNCTION NAME : split_o2_isotopes()

    DESCRIPTION : Function to split O2 into its 3 main isotopes

    INPUTS : 

        Atmosphere :: archNEMESIS atmosphere class

    OPTIONAL INPUTS:
    
        o18ratio :: Value of 18O/16O (VSMOW)
        o17ratio :: Value of 17O/16O (VSMOW)

    OUTPUTS : 
 
        Atmosphere :: archNEMESIS Atmosphere class

    CALLING SEQUENCE:

        Atmosphere = split_O2_isotopes(Atmosphere,o18ratio=1.,o17ratio=1.)

    MODIFICATION HISTORY : Juan Alday (18/06/2026)

    """

    #Splitting CO2 into 4 isotopes
    io2 = np.where( (Atmosphere.ID==7) & (Atmosphere.ISO==0) )[0][0]

    vmr_o2 = Atmosphere.VMR[:,io2]

    Atmosphere.remove_gas(7,0)

    #Adding the isotopes of CO2 isotopes (following HITRAN fractionation)
    Atmosphere.add_gas(7,1,vmr_o2*0.995262)
    Atmosphere.add_gas(7,2,vmr_o2*0.003991*o18ratio)       #VSMOW
    Atmosphere.add_gas(7,3,vmr_o2*7.422350e-4*o17ratio)    #VSMOW

    return Atmosphere

###########################################################################################################################

def read_solar_transmission_toon():
    """
    FUNCTION NAME : read_solar_transmission_toon()

    DESCRIPTION : Read the solar transmission data from:
                    Toon, G. C., Solar line list for GGG2014, TCCON data archive, 
                    hosted by the Carbon Dioxide Information Analysis Center, 
                    Oak Ridge National Laboratory, Oak Ridge, Tennessee, U.S.A., 
                    doi:10. 14291/tccon.ggg2014.solar.R0/1221658, 2014.

    INPUTS : 

        None

    OPTIONAL INPUTS:
    
        None

    OUTPUTS : 
 
        wave :: Wavenumber array (cm-1)
        trans :: Dimensionless solar pseudo-transmittance

    CALLING SEQUENCE:

        wave,trans = read_solar_transmission_toon()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    #Reading the Toon solar transmission data
    toon_file = os.path.join(archground.paths.texespy_path,'data','solar','solar_merged_20240731_600_33300_100.out')
    data = np.loadtxt(toon_file, skiprows=3)
    wave = data[:,0] 
    trans = data[:,1]

    return wave,trans