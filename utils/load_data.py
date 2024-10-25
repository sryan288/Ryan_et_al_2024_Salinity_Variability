'''
 Routines to load datasets
'''
'''
    1) load_pioneer()    Loads all pioneer datasets downloaded with python/python_poseidon/Salinity_Intrusion/download_pioneer_data.ipynb
    
    2) load_gulf_stream()    Reads mat file from Magdalena Andres and creates xarray
'''

import xarray as xr


#
#------------------------------------------------------------------------
# 1) load pioneer data

def load_pioneer():
    datapath = '/mnt/data/OOI_Pioneer/'
    inshore7m = xr.open_dataset(datapath + 'pioneer_inshore7m.nc')
    inshore_prof = xr.open_dataset(datapath + 'pioneer_inshore_prof.nc')
    inshore_surf = xr.open_dataset(datapath + 'pioneer_inshore_surf.nc')
    inshore_surf = inshore_surf.rename({'sea_surface_temperature (degree_Celsius)':'sea_water_temperature (degree_Celsius)'})
    ups_inshore_prof = xr.open_dataset(datapath + 'pioneer_ups_inshore_prof.nc')
    central_inshore_prof = xr.open_dataset(datapath + 'pioneer_central_inshore_prof.nc')
    offshore7m = xr.open_dataset(datapath + 'pioneer_offshore7m.nc')
    central7m = xr.open_dataset(datapath + 'pioneer_central7m.nc')
    central_surf = xr.open_dataset(datapath + 'pioneer_central_surf.nc')
    central_surf = central_surf.rename({'sea_surface_temperature (degree_Celsius)':'sea_water_temperature (degree_Celsius)'})
    offshore_surf = xr.open_dataset(datapath + 'pioneer_offshore_surf.nc')
    offshore_surf = offshore_surf.rename({'sea_surface_temperature (degree_Celsius)':'sea_water_temperature (degree_Celsius)'})
    central_offshore_prof = xr.open_dataset(datapath + 'pioneer_central_offshore_prof.nc')
    offshore_prof = xr.open_dataset(datapath + 'pioneer_offshore_prof.nc')
    
    ###########################################################################
    # create dictionary for profilers and surface moorings
    prof = {}
    prof['inshore'] = inshore_prof
    prof['upstream_inshore'] = ups_inshore_prof
    prof['central_inshore'] = central_inshore_prof
    prof['central_offshore'] = central_offshore_prof
    prof['offshore'] = offshore_prof

    surf7m = {}
    surf7m['inshore'] = inshore7m
    surf7m['central'] = central7m
    surf7m['offshore'] = offshore7m
    # surf['central_offshore_prof'] = central_offshore_prof
    # surf['offshore_prof'] = offshore_prof

    surf = {}
    surf['inshore'] = inshore_surf
    surf['central'] = central_surf
    surf['offshore'] = offshore_surf
    
    
    ## put all into one dictionary
    pioneer = {}
    pioneer['surf'] = surf
    pioneer['surf7m'] = surf7m
    pioneer['prof'] = prof
    return pioneer