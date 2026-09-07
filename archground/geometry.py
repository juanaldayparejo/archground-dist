#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
# archGROUND - Python package for ground-based planetary spectroscopy and radiative transfer with archNEMESIS. 
# geometry.py - Set of functions to compute geometry for ground-based observations.
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

from astropy.time import Time
from astroquery.jplhorizons import Horizons
import archground
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import patches

#############################################################################################################
#############################################################################################################
# DICTIONARIES
#############################################################################################################
#############################################################################################################

planet_dict = {
    'MERCURY': {
        'horizons_id': "199",
        'equatorial_radius_km': 2440.53,
        'polar_radius_km': 2438.26,
    },
    'VENUS': {
        'horizons_id': "299",
        'equatorial_radius_km': 6051.80,
        'polar_radius_km': 6051.80,
    },
    'EARTH': {
        'horizons_id': "399",
        'equatorial_radius_km': 6378.1366,
        'polar_radius_km': 6356.7519,
    },
    'MARS': {
        'horizons_id': "499",
        'equatorial_radius_km': 3396.19,
        'polar_radius_km': 3376.20,
    },
    'JUPITER': {
        'horizons_id': "599",
        'equatorial_radius_km': 71492.0,
        'polar_radius_km': 66854.0,
    },
    'SATURN': {
        'horizons_id': "699",
        'equatorial_radius_km': 60268.0,
        'polar_radius_km': 54364.0,
    },
    'URANUS': {
        'horizons_id': "799",
        'equatorial_radius_km': 25559.0,
        'polar_radius_km': 24973.0,
    },
    'NEPTUNE': {
        'horizons_id': "899",
        'equatorial_radius_km': 24764.0,
        'polar_radius_km': 24341.0,
    },
}

observatory_dict = {
    'IRTF': {
        "name": "NASA Infrared Telescope Facility",
        "horizons_id": "568",
        "lon": -155.472000,
        "lat": 19.8262575,
        "elevation": 4.20124,
    },
    'CAHA': {
        "name": "Calar Alto Observatory",
        "horizons_id": "493",
        "lon": -2.545800,
        "lat": 37.2234714,
        "elevation": 2.17561,
    },
}


#############################################################################################################
#############################################################################################################
# JPL HORIZONS QUERY FUNCTIONS
#############################################################################################################
#############################################################################################################

def get_ephemerides_epoch(date,time,planet,observatory):
    """
    Get ephemerides for a given date, time, planet and observatory.

    Parameters
    ----------
    date : str
        Date in the format 'YYYY-MM-DD'.
    time : str
        Time in the format 'HH:MM:SS'.
    planet : str
        Name of the planet (see archground dictionaries).
    observatory : str
        Name of the observatory (see archground dictionaries).

    Returns
    -------
    dict
        Dictionary containing ephemerides information.
    """
    
    # Combine date and time into a single string
    datetime_str = f"{date} {time}"
    
    # Convert to astropy Time object
    epoch = Time(datetime_str, scale="utc").jd
    
    # Get the Horizons ID for the planet and observatory
    planet_id = planet_dict[planet.upper()]["horizons_id"]
    observatory_id = observatory_dict[observatory.upper()]["horizons_id"]

    #Query the JPL Horizons system for ephemerides
    query = Horizons(
        id=planet_id,    
        location=observatory_id, 
        epochs=epoch,
    )

    eph = query.ephemerides(quantities="1,2,3,4,8,10,13,14,15,16,17,19,20,21,23,24,44")

    return eph

def planetodetic_to_planetocentric(
    latitude_deg,
    equatorial_radius=3396.19,
    polar_radius=3376.20,
):
    """Convert latitude in degrees for a point on the reference ellipsoid.

    Radii must use the same units. Accepts scalars or arrays.
    """
    lat = np.deg2rad(latitude_deg)
    ratio_squared = (polar_radius / equatorial_radius)**2

    return np.rad2deg(
        np.arctan2(ratio_squared * np.sin(lat), np.cos(lat))
    )

##############################################################################################
##############################################################################################
#                                     MAPPING GEOMETRY
##############################################################################################
##############################################################################################


def build_maps_epoch(ephemerides,res_bin=0.1):
    '''
    FUNCTION NAME : build_maps_epoch()

    DESCRIPTION :   Function to create an angular grid of pixels (in arcsec) and calculate the associated
                    latitude and longitude of the body disk observed in the centre of the grid
                    
                    Formulas were taken from https://en.wikipedia.org/wiki/Orthographic_map_projection

    INPUTS : 

        ephemerides :: dict
            Dictionary containing the ephemerides from the JPL Horizons System
        
    OPTIONAL INPUTS:
    
        res_bin :: Resolution of the map (arcsec/bin)
            
    OUTPUTS : 
 
        maps :: dict
            Dictionary containing the maps with the geometry

    CALLING SEQUENCE:

        maps = build_maps_epch(ephemerides,res_bin=0.1)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    '''

    #Getting the relevant parameters from the ephemerides
    angsize = ephemerides["ang_width"][0]          #Angular size of the planet's disk (arcsec)
    pa = ephemerides["NPole_ang"][0]               #Position angle of the North pole wrt the astronomical north
    lat_subobs = ephemerides["PDObsLat"][0]        #Planetodetic sub-observer latitude (degrees)
    lon_subobs = (-ephemerides["PDObsLon"][0] + 180.0) % 360.0 - 180.0        #Planetodetic sub-observer longitude (degrees) - East positive
    lat_subsol = ephemerides["PDSunLat"][0]        #Planetodetic sub-solar latitude (degrees)
    lon_subsol = (-ephemerides["PDSunLon"][0] + 180.0) % 360.0 - 180.0        #Planetodetic sub-solar longitude (degrees) - East positive
    phase_angle = ephemerides["alpha"][0]

    #Changing from planetodetic to planetocentric
    target = str(ephemerides["targetname"][0]).split("(")[0].strip().upper()   #Name of planet
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

    #Calculating the maps
    angmax = np.ceil(angsize*0.75)
    res = res_bin
    rotan = -pa
    hx,hy,lat,lon = archground.geometry.xy2latlon(angmax,res,angsize,lat_subobs,lon_subobs,rot_angle=rotan)
    emiss_ang = archground.geometry.latlon2emiss(lat,lon,lat_subobs,lon_subobs)
    sol_ang = archground.geometry.latlon2sza(lat,lon,lat_subsol,lon_subsol)
    azi_ang = archground.geometry.phase2azi(emiss_ang,sol_ang,phase_angle)
    lst = archground.geometry.lon2lst(lon, lon_subsol)

    #Creating dictionary
    maps = {
        "hx" : hx, #(ny,nx) Position across x-axis (arcsec)
        "hy" : hy, #(ny,nx) Position across y-axis (arcsec)
        "lat" : lat, #(ny,nx) Latitude of each bin (degrees)
        "lon" : lon, #(ny,nx) Longitude of each bin (degrees)
        "emiss_ang" : emiss_ang, #(ny,nx) Emission angle of each bin (degrees)
        "sol_ang" : sol_ang, #(ny,nx) Solar zenith angle of each bin (degrees)
        "azi_ang" : azi_ang, #(ny,nx) Azimuth angle of each bin (degrees)
        "lst" : lst, #(ny,nx) Local time of each bin (degrees)
    }

    return maps

############################################################################################################

def plot_maps(maps,
              add_fibre_plot=False,x_fibre=0.,y_fibre=0.,d_fibre=1.5):
    '''
    FUNCTION NAME : plot_maps()

    DESCRIPTION :   Function to plot the geometry maps computed with build_maps_epoch()

    INPUTS : 

        maps :: dict
            Dictionary containing the geometry maps for the selected epoch
        
    OPTIONAL INPUTS:

        add_fibre_plot :: If True, it adds a circular fibre to the plot, according to x, y and d_fibre
            
    OUTPUTS : 
 
        plots

    CALLING SEQUENCE:

        plot_maps(maps)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    '''

    #Plotting some parameters
    fig,([ax1,ax2],[ax3,ax4]) = plt.subplots(2,2,figsize=(8,6))

    cmap = 'turbo'

    nlevels = 51

    im1 = ax1.contourf(maps["hx"],maps["hy"],maps["lat"],cmap=cmap,levels=np.linspace(-90.,90.,nlevels),vmin=-90.,vmax=90.)
    add_compass(ax1)
    if add_fibre_plot is True:
        add_fibre(
            ax1,
            x_fibre,
            y_fibre,
            d_fibre,
            edgecolor="black",
            facecolor="none",
            linewidth=1.0,
            alpha=1.0,
            label=None,
            zorder=10,
        )
    ax1.grid()

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im1, cax=cax)
    cbar.set_label('Latitude ($^\circ$)')

    im2 = ax2.contourf(maps["hx"],maps["hy"],maps["lst"],cmap=cmap,levels=np.linspace(0.,24.,nlevels),vmin=0.,vmax=24.)
    add_compass(ax2)
    if add_fibre_plot is True:
        add_fibre(
            ax2,
            x_fibre,
            y_fibre,
            d_fibre,
            edgecolor="black",
            facecolor="none",
            linewidth=1.0,
            alpha=1.0,
            label=None,
            zorder=10,
        )
    ax2.grid()

    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im2, cax=cax)
    cbar.set_label('LST (hour)')


    im3 = ax3.contourf(maps["hx"],maps["hy"],maps["emiss_ang"],cmap=cmap+'_r',levels=np.linspace(0.,90.,nlevels),vmin=0.,vmax=90.)
    add_compass(ax3)
    if add_fibre_plot is True:
        add_fibre(
            ax3,
            x_fibre,
            y_fibre,
            d_fibre,
            edgecolor="black",
            facecolor="none",
            linewidth=1.0,
            alpha=1.0,
            label=None,
            zorder=10,
        )
    ax3.grid()

    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im3, cax=cax)
    cbar.set_label('Emission angle ($^\circ$)')

    im4 = ax4.contourf(maps["hx"],maps["hy"],maps["sol_ang"],cmap=cmap+'_r',levels=np.linspace(0.,180.,nlevels),vmin=0.,vmax=180.)
    add_compass(ax4)
    if add_fibre_plot is True:
        add_fibre(
            ax4,
            x_fibre,
            y_fibre,
            d_fibre,
            edgecolor="black",
            facecolor="none",
            linewidth=1.0,
            alpha=1.0,
            label=None,
            zorder=10,
        )
    ax4.grid()

    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im4, cax=cax)
    cbar.set_label('Solar zenith angle ($^\circ$)')

    plt.tight_layout()
            


###############################################################################################

def xy2latlon(angmax,res,angsize_mars,lat_subobs,lon_subobs,rot_angle=0.):
    
    '''
    FUNCTION NAME : xy2latlon()

    DESCRIPTION :   Function to create an angular grid of pixels (in arcsec) and calculate the associated
                    latitude and longitude of the body disk observed in the centre of the grid
                    
                    Formulas were taken from https://en.wikipedia.org/wiki/Orthographic_map_projection

    INPUTS : 

        angmax :: Maximum size to include in the image (image goes from -angmax to angmax) (arcsec)
        res :: Resolution of the image (arcsec/pix)
        angsize :: Angular size of the observed body (arcsec)
        lat_subobs :: Sub-observer latitude (degrees)
        lon_subobs :: Sub-observer lontitude (degrees)
        
    OPTIONAL INPUTS:
    
        rot_angle :: Angle to rotate the disk to align it with the observation (0 angles in North and positive angle is to the East)
            
    OUTPUTS : 
 
        hx :: Angular position of each pixel in x direction (arcsec)
        hy :: Angular position of each pixel in y direction (arcsec)
        lat :: Latitude of each point (degrees)
        lon :: Longitude of each point (degrees)

    CALLING SEQUENCE:

        hx,hy,lat,lon = xy2latlon(angmax,res,angsize,lat_subobs,lon_subobs)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    
    '''

    lat0 = np.radians(lat_subobs)  # Convert degrees to radians
    lon0 = np.radians(lon_subobs)

    # Image pixel coordinates
    px = np.arange(-angmax, angmax + res, res)
    py = np.arange(-angmax, angmax + res, res)
    hx, hy = np.meshgrid(px, py)

    # If rotation is needed, apply a proper 2D rotation matrix
    if rot_angle != 0.0:
        theta = np.radians(rot_angle)  # Convert degrees to radians
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)

        # Apply rotation matrix
        hxrot = hx * cos_theta - hy * sin_theta
        hyrot = hx * sin_theta + hy * cos_theta
    else:
        hxrot = hx
        hyrot = hy

    # Compute z of sphere hit position, if pixel's ray hits
    rho = np.sqrt(hxrot**2 + hyrot**2)
    c = np.ones(rho.shape)
    c[rho<=angsize_mars/2.] = np.arcsin(rho[rho<=angsize_mars/2.] / (angsize_mars / 2.0))

    # Compute latitude and longitude
    lat = np.arcsin(np.cos(c) * np.sin(lat0) + (hyrot * np.sin(c) * np.cos(lat0)) / rho)
    lon = lon0 + np.arctan2(hxrot * np.sin(c), (rho * np.cos(c) * np.cos(lat0) - hyrot * np.sin(c) * np.sin(lat0)))

    # Handle singularities (when rho == 0)
    ic = np.where(rho == 0.)
    if len(ic[0]) > 0:
        lat[ic] = lat0
        lon[ic] = lon0

    # Convert back to degrees
    lat = np.degrees(lat)
    lon = np.degrees(lon)

    # Adjust longitudes to stay within [-180,180]
    lon[lon <= -180] += 360
    lon[lon > 180] -= 360

    lat[hx**2. + hy**2. > (angsize_mars/2.)**2.] = np.nan
    lon[hx**2. + hy**2. > (angsize_mars/2.)**2.] = np.nan

    return hx, hy, lat, lon


###############################################################################################

def latlon2emiss(lat,lon,lat_subobs,lon_subobs):
    
    '''
    FUNCTION NAME : latlon2emiss()

    DESCRIPTION :   Given the latitude and longitude of a grid of pixels, and the sub-observer latitude and longitude
                    it calculates the emission angle at each pixel
                    
                    Formulas were taken from https://math.stackexchange.com/questions/2688803/angle-between-two-points-on-a-sphere

    INPUTS : 

        lat :: Latitude of each pixel in the image (deg)
        lon :: Longitude of each pixel in the image (deg)
        lat_subobs :: Sub-observer latitude (deg)
        lon_subobs :: Sub-observer longitude (deg)
    
    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        emiss_ang :: Emission angle of each point (degrees)

    CALLING SEQUENCE:

        emiss_ang = latlon2emiss(lat,lon,lat_subobs,lon_subobs)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    
    '''
    
    #Changing from degrees to radians
    lat = lat / 180. * np.pi
    lon = lon / 180. * np.pi
    lat_subobs = lat_subobs / 180. * np.pi
    lon_subobs = lon_subobs / 180. * np.pi
    
    #Calculating the angle between each point and the sub-observer point
    #emiss_ang = np.arccos( np.sin(lon)*np.sin(lon_subobs) + np.cos(lon)*np.cos(lon_subobs)*np.cos(lat-lat_subobs) ) / np.pi * 180.
    emiss_ang = np.arccos( np.sin(lat)*np.sin(lat_subobs) + np.cos(lat)*np.cos(lat_subobs)*np.cos(lon-lon_subobs) ) / np.pi * 180.
    
    return emiss_ang

###############################################################################################

def latlon2sza(lat,lon,lat_subsol,lon_subsol):
    
    '''
    FUNCTION NAME : latlon2sza()

    DESCRIPTION :   Given the latitude and longitude of a grid of pixels, and the sub-solar latitude and longitude
                    it calculates the incident angle at each pixel
                    
                    Formulas were taken from https://math.stackexchange.com/questions/2688803/angle-between-two-points-on-a-sphere

    INPUTS : 

        lat :: Latitude of each pixel in the image (deg)
        lon :: Longitude of each pixel in the image (deg)
        lat_subsol :: Sub-solar latitude (deg)
        lon_subsol :: Sub-solar longitude (deg)
    
    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        sol_ang :: Incident solar angle of each point (degrees)

    CALLING SEQUENCE:

        sol_ang = latlon2sza(lat,lon,lat_subsol,lon_subsol)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    
    '''
    
    #Changing from degrees to radians
    lat = lat / 180. * np.pi
    lon = lon / 180. * np.pi
    lat_subsol = lat_subsol / 180. * np.pi
    lon_subsol = lon_subsol / 180. * np.pi
    
    #Calculating the angle between each point and the sub-observer point
    #emiss_ang = np.arccos( np.sin(lon)*np.sin(lon_subobs) + np.cos(lon)*np.cos(lon_subobs)*np.cos(lat-lat_subobs) ) / np.pi * 180.
    sza = np.arccos( np.sin(lat)*np.sin(lat_subsol) + np.cos(lat)*np.cos(lat_subsol)*np.cos(lon-lon_subsol) ) / np.pi * 180.
    
    return sza


###############################################################################################

def phase2azi(emiss_ang,sza,phase):
    
    '''
    FUNCTION NAME : phase2azi()

    DESCRIPTION :   Given the emission and incident angles of a grid of pixels, together with the 
                    phase angle it calculates the azimuth angle at each pixel. The azimuth angle here
                    follows the convention used in NEMESIS (phi=0 for forward scattering)

    INPUTS : 

        emiss_ang :: Emission angle (deg)
        sza :: Solar zenith angle (deg)
        phase :: Phase angle (deg)
    
    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        azi_ang :: Azimuth angle of each point (degrees)

    CALLING SEQUENCE:

        azi_ang = phase2azi(emiss_ang,sza,phase)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    
    '''
    
    #First of all let's calculate the scattering phase angle
    mu = np.cos(emiss_ang/180.*np.pi)   #Cosine of the reflection angle
    mu0 = np.cos(sza/180.*np.pi)    #Coside of the incidence angle
    cg = np.cos(phase/180.*np.pi)    #Cosine of the phase angle
    
    cazi = ((mu * mu0 - cg)/(np.sqrt(1. - mu**2.) * np.sqrt(1.-mu0**2.)))
    
    iin = np.where( (np.isnan(emiss_ang)==False) & (cazi<=-1.0) )
    cazi[iin] = -1.0
    
    iin = np.where( (np.isnan(emiss_ang)==False) & (cazi>=1.0) )
    cazi[iin] = 1.0
    
    iin = np.where( (np.isnan(cazi)==True) )
    cazi[iin] = 0.0
    
    azi_ang = np.arccos(cazi) / np.pi * 180.
    
    return azi_ang

###############################################################################################

def lon2lst(lon, lon_subsol):
    '''
    FUNCTION NAME : lon2lst()

    DESCRIPTION : Given the longitude in a grid of pixels and the sub-solar longitude, it 
                  calculates the Local Solar Time in each pixel.

    INPUTS : 
        lon :: Longitude of each pixel (deg)
        lon_subsol :: Sub-solar longitude (deg)
    
    OUTPUTS : 
        LST :: Local solar time (h)

    CALLING SEQUENCE:
        lst = lon2lst(lon, lon_subsol)

    MODIFICATION HISTORY : Juan Alday (24/04/2024)
    '''
    
    # Calculating the Universal Local Time (LST at lon=0)
    ULT = 12. - lon_subsol / 180. * 12.
    
    # Compute Local Solar Time
    lst = ULT + lon / 180. * 12.
    
    # Ensure LST remains in the range [0, 24]
    lst = np.mod(lst, 24.)

    return lst

###############################################################################################

# Convert RA from degrees to hours, minutes, seconds
def convert_ra_to_hms(ra_deg):
    '''
    Convert the right ascensions to hours, minutes, second
    '''
    ra_hours = ra_deg / 15
    h = int(ra_hours)
    m = int((ra_hours - h) * 60)
    s = ((ra_hours - h) * 60 - m) * 60
    return f"{h}h {m}m {s:.2f}s"

###############################################################################################

# Convert Dec to degrees, arcminutes, arcseconds
def convert_dec_to_dms(dec_deg):
    '''
    Convert declination for float to degrees,arcminutes,arcseconds
    '''
    d = int(dec_deg)
    m = int(abs(dec_deg - d) * 60)
    s = (abs(dec_deg - d) * 60 - m) * 60
    sign = "+" if dec_deg >= 0 else "-"
    return f"{sign}{abs(d)}° {m}′ {s:.2f}″"

###############################################################################################

def add_compass(ax,compass_x=0.15,compass_y=0.15,size=0.075):
    '''
    Add a compass with sky cardinal points to a given axis
    '''

    from matplotlib.patches import FancyArrow
    
    ax.annotate("N", xy=(compass_x, compass_y + size), xycoords='axes fraction',
                fontsize=10, ha="center", va="bottom")
    
    ax.annotate("S", xy=(compass_x, compass_y - size), xycoords='axes fraction',
                fontsize=10, ha="center", va="top")
    
    ax.annotate("W", xy=(compass_x + size, compass_y), xycoords='axes fraction',
                fontsize=10, ha="left", va="center")
    
    ax.annotate("E", xy=(compass_x - size, compass_y), xycoords='axes fraction',
                fontsize=10, ha="right", va="center")
    
    # Adding arrows for the compass
    ax.add_patch(FancyArrow(compass_x, compass_y, 0, size*0.7, transform=ax.transAxes, width=0.005))
    ax.add_patch(FancyArrow(compass_x, compass_y, 0, -size*0.7, transform=ax.transAxes, width=0.005))
    ax.add_patch(FancyArrow(compass_x, compass_y, size*0.7, 0, transform=ax.transAxes, width=0.005))
    ax.add_patch(FancyArrow(compass_x, compass_y, -size*0.7, 0, transform=ax.transAxes, width=0.005))

###############################################################################################

def add_slit(
    ax,
    x_center,
    y_center,
    length,
    width,
    angle=0.0,
    edgecolor="black",
    facecolor="none",
    linewidth=1.0,
    alpha=1.0,
    label=None,
    zorder=10,
):
    """
    Add a rectangular slit footprint to a Matplotlib axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis on which the slit will be drawn.
    x_center, y_center : float
        Center of the slit in the coordinate system of the plot.
    length : float
        Slit dimension along the y-axis
    width : float
        Slit dimension along the x-axis
    angle : float, optional
        Counterclockwise rotation in degrees. An angle of 0° produces
        a horizontal slit; 90° produces a vertical slit.

    Returns
    -------
    slit : matplotlib.patches.Rectangle
        The created rectangle.
    """
    slit = patches.Rectangle(
        (x_center - width / 2, y_center - length / 2),
        width=width,
        height=length,
        angle=angle,
        rotation_point="center",
        edgecolor=edgecolor,
        facecolor=facecolor,
        linewidth=linewidth,
        alpha=alpha,
        label=label,
        zorder=zorder,
    )

    ax.add_patch(slit)
    return slit

###############################################################################################

def add_fibre(
    ax,
    x_center,
    y_center,
    diameter,
    edgecolor="black",
    facecolor="none",
    linewidth=1.0,
    alpha=1.0,
    label=None,
    zorder=10,
):
    """
    Add a circular fibre footprint to a Matplotlib axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis on which the fibre will be drawn.
    x_center, y_center : float
        Centre of the fibre in plot coordinates.
    diameter : float
        Fibre diameter in the same units as the plot coordinates.

    Returns
    -------
    fibre : matplotlib.patches.Circle
        The created circle.
    """
    fibre = patches.Circle(
        (x_center, y_center),
        radius=diameter / 2,
        edgecolor=edgecolor,
        facecolor=facecolor,
        linewidth=linewidth,
        alpha=alpha,
        label=label,
        zorder=zorder,
    )

    ax.add_patch(fibre)
    return fibre

###############################################################################################
