#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 09:40:35 2024

@author: sryan

----------------------
FIGURE 8
Example of additive model at Pioneer
----------------------

"""

import os
# go to working directory
workingpath = '/home/sryan/python/NASA_salinity/p2_salinity_variability/'
os.chdir(workingpath)
# filepath = workingpath+'/scripts/'
# printpath = filepath+'woa_box_analysis.py'
STARTUP_2023_coldfronts_mhws = "./utils/2023_salinity_variability_startup.py"
exec(open(STARTUP_2023_coldfronts_mhws).read())
# my own packages
import sys
sys.path.append("./") # go to parent dir
from utils.plot_utils import finished_plot

#%% Data processing
# load save data further below for plotting

#=============================================================================
# #%% Decomposition: Pioneer array data
# # =============================================================================
# #

# ## load pioneer data 
# pioneer_surf = np.load('/mnt/data/OOI_Pioneer/moor_surf_xr.npy',allow_pickle=True)


# #%%% 1) processes pioneer data, i.e. fill gaps and apply low-pass filter, derive 
# #seasonal cycle by group by
# # Data loaded at beginning of script
# surf_lp30 = {}
# surf_lp30_mm = {}
# surf_lp30_seas = {}

# # inshore and offshore mooring only have small gaps than can be interpolates, however
# # central mooring has large gaps across seasonal maximum, so we handle it separately
# for name in ['issm','ossm']:
#     ds = pioneer_surf[name]['surf']['sal'].resample(time='1D').mean(skipna=True)
#     ds = ds.interpolate_na(dim="time", method="linear").to_dataset()
#     for var in ['sal','temp','dens']:
#         dummy = pioneer_surf[name]['surf'][var].resample(time='1D').mean(skipna=True)
#         dummy = dummy.interpolate_na(dim="time", method="linear").to_dataset()
        
#         wo = np.where(np.isfinite(dummy[var]))[0]
#         ds_filt = xr.DataArray(data=butterworth_lowpass_filter(dummy[var][wo], order=2, cutoff_freq=1.0/30, axis=0),
#                                 dims=["time"],
#                                 coords=dict(time=dummy['time'][wo].values))
#                                 #30 day lp filter
#         ds[var] =  ds_filt.interp(time=ds.time)

#     surf_lp30[name] = ds
#     surf_lp30_mm[name] = surf_lp30[name].resample(time='MS').mean() 
#     ds_seas = cut2fullyears(surf_lp30_mm[name])
#     surf_lp30_seas[name] = ds_seas.groupby('time.month').mean('time')
    
# # I am now filling the large gaps at cnsm with the average of issm and ossm but
# # for filtering and am then setting the original values back to NaN
# for name in ['cnsm']:
#     ds = pioneer_surf[name]['surf']['sal'].resample(time='1D').mean(skipna=True)
#     wonan = np.where(np.isnan(ds))
#     ds = ds.interpolate_na(dim="time", method="linear").to_dataset()
#     for var in ['sal','temp','dens']:
#         dummy = pioneer_surf[name]['surf'][var].resample(time='1D').mean(skipna=True).to_dataset()
#         dummy_fill = (pioneer_surf['issm']['surf'][var].resample(time='1D').mean(skipna=True)+
#                     pioneer_surf['ossm']['surf'][var].resample(time='1D').mean(skipna=True))/2
#         dummy[var][wonan[0]] = dummy_fill[wonan[0]]
#         # dummy = dummy.interpolate_na(dim="time", method="cubic")
#         wo = np.where(np.isfinite(dummy[var]))[0]
#         ds_filt = xr.DataArray(data=butterworth_lowpass_filter(dummy[var][wo], order=2, cutoff_freq=1.0/30, axis=0),
#                                 dims=["time"],
#                                 coords=dict(time=dummy['time'][wo].values))
#                                 #30 day lp filter
#         ds[var] =  ds_filt.interp(time=ds.time)
#         # ds[var][wonan] = np.nan  # uncomment if you want to add original gaps back in

#     surf_lp30[name] = ds
#     surf_lp30_mm[name] = surf_lp30[name].resample(time='MS').mean() 
#     ds_seas = cut2fullyears(surf_lp30_mm[name])
#     surf_lp30_seas[name] = ds_seas.groupby('time.month').mean('time')
# #  
# # #%%% load saved data
# # ## save to pickle 
# # surf_lp30 = np.load('./data/pioneer_TSrho_surf_lp30.npy',allow_pickle=True)
# # surf_lp30_mm = np.load('./data/pioneer_TSrho_surf_lp30_mm.npy',allow_pickle=True)
# # surf_lp30_seas = np.load('./data/pioneer_TSrho_surf_lp30_seas_fullyears.npy',allow_pickle=True)


# =============================================================================
# %% Load saved data data
# =============================================================================
surf_lp30 = np.load('./data/pioneer_TSrho_surf_lp30.npy',allow_pickle=True)
 
#%% 1) decompose timeseries into trend, seasonal and residual components using
#    statsmodels Python package
import statsmodels.api as sm  # https://www.statsmodels.org --> code by chat GPT
plt.close('all')
font_for_pres()
cols = cmap = plt.get_cmap('tab10')
fig,ax = plt.subplots(nrows=4,figsize=(6, 6),sharex=False,gridspec_kw={'height_ratios': [1.2, 0.9,0.8,1]})
for name,i in zip(['issm','ossm'],range(2)):
    
    labelname = ['Pioneer Inshore','Pioneer Offshore']
    
    # define your data and change to pandas dataframe
    x = surf_lp30[name].sel(time=slice('2015-01','2022-12')).time
    y = surf_lp30[name].sel(time=slice('2015-01','2022-12'))['sal']
    wo = np.where(np.isfinite(y.values))[0]
    data = pd.DataFrame({'Date': x[wo], 'Value': y[wo]})
    # Decompose the time series data
    result = sm.tsa.seasonal_decompose(data['Value'], model='additive', period=365)  # Assuming a seasonal cycle of 12 months (365 for daily values)
    
    # plot
    # fig,ax = plt.subplots(nrows=4,figsize=(8, 6),sharex=True)
    ax[0].plot(data['Date'], data['Value'], label=labelname[i],linewidth=1)
    ax[0].legend(loc='lower left',ncol=3,bbox_to_anchor=(0.02, -0.15))
    ax[0].xaxis.set_ticks_position('top')
    ax[0].xaxis.set_label_position('top')
    ax[0].spines['bottom'].set_visible(False)
    ax[0].minorticks_off()
    ax[0].set_ylim(31,36.5)
    ax[0].set_ylabel('full salinity',labelpad=18)
    
    # ax[0].set_yticks([32,34,36])
    ax[1].plot(data['Date'], result.trend-np.nanmean(result.trend),linewidth=1)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].axhline(0,color='k',linestyle='dashed',linewidth=0.5,zorder=0)
    ax[1].set_xticks([])
    ax[1].set_ylabel('low-frequency \n component',labelpad=0)
    
    # ax[1].set_title('Trend')
    ax[2].plot(data['Date'], result.seasonal,linewidth=1)
    ax[2].spines['top'].set_visible(False)
    ax[2].spines['bottom'].set_visible(False)
    ax[2].set_xticks([])
    ax[2].axhline(0,color='k',linestyle='dashed',linewidth=0.5,zorder=0)
    ax[2].set_ylabel('seasonal \n component',labelpad=18)


    # ax[2].set_title('Seasonal')
    ax[3].plot(data['Date'], result.resid,linewidth=1)
    # ax[3].set_title('Residual')
    ax[3].axhline(0,color='k',linestyle='dashed',linewidth=0.5,zorder=0)
    ax[3].spines['top'].set_visible(False)
    ax[3].set_ylabel('residuals',labelpad=16)
    
    [axh.grid(False) for axh in ax]
    [axh.set_xlim(np.datetime64('2015-04-01'),np.datetime64('2022-06-30')) for axh in ax]
    
# finished_plot(fig,'./manuscript/plots/Figure8_pioneer_additive_model.png')
