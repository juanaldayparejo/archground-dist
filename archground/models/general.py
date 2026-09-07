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

def create_stellar_class(dist_sun):
    """
    FUNCTION NAME : create_stellar_class()

    DESCRIPTION : Create an instance of the `Stellar_0` class containing the stellar parameters for the any case

    INPUTS : 

        dist_sun :: Distance of planet from the Sun in AU

    OPTIONAL INPUTS:
    
        None
            
    OUTPUTS : 
 
        Stellar :: An instance of the `Stellar_0` class containing the stellar parameters for the any case

    CALLING SEQUENCE:

        Stellar = create_stellar_class()

    MODIFICATION HISTORY : Juan Alday (13/03/2025)
    """

    Stellar = ans.Stellar_0()

    #Defining the Sun-planet distance
    Stellar.DIST = dist_sun

    #Defining the file containing the solar spectrum
    solfile = 'houghtonsolarwn.dat'

    return Stellar, solfile

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
    Layer.LAYHT = Atmosphere.H[0] #m
    Layer.assess()

    return Layer


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
    Atmosphere.add_gas(2,4,vmr_co2*7.339890e-4*o18ratio)    #VSMOW

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