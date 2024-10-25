#!/usr/bin/env python
# coding: utf-8

#
import os
# go to working directory
workingpath = '/home/sryan/python/NASA_salinity/p2_salinity_variability/'
os.chdir(workingpath)
filepath = workingpath+'/scripts/'
printpath = filepath+'salinity_analysis.py'
# 
# An overview - https://salinity.oceansciences.org/data-salinity-01.htm
STARTUP_2023_coldfronts_mhws = "./utils/2023_salinity_variability_startup.py"
exec(open(STARTUP_2023_coldfronts_mhws).read())

import sys
sys.path.append("~/home/sryan/python/") # go to parent dir
from utils.plot_utils import finished_plot,latlon_label,datetime_tick,calculate_ticks
# from utils.orca_utilities import deseason
from utils.xarray_utilities import xr_lagged_pearson
from utils.datafun import butterworth_lowpass_filter,matlab2datetime
from utils.general import pickle_save

# Some definitions
## overlapping period
time_overlap = ['2016-01-01','2022-12-31']
## region to plot on maps
plot_region = [-77,-63,36,45]#[-78,-68,35,43]


# =============================================================================
#%% Load data processed in prepare_satellite_surface_salinity.py
# =============================================================================

## sea surface salinity data preprocesses for three products
sss_monthly = np.load('./data/sss_monthly.npy',allow_pickle=True)
# sss_daily = np.load('./data/sss_daily.npy',allow_pickle=True)


## load pioneer data prepared by prepare_pioneer.ipynb
pioneer_surf = np.load('/mnt/data/OOI_Pioneer/moor_surf_xr.npy',allow_pickle=True)


## wcr ring trajectories
wcr = xr.open_dataset('/mnt/data/WCRs/wcrs_2011_2020_silver_et_al_2022.nc')


# =============================================================================
#%% Maps Seasonal cycle SSS
# =============================================================================
#
#%%% Figure 4 + SI: Maps mean seasonal cycle of overlapping period
##########################################

font_for_pres()
plt.close('all')
for name,i in zip(['oisss','smap','smos'],range(3)):
#     plt.close('all')
    ds = sss_monthly[name].sel(time=slice(*time_overlap))
    fs=11
    plt.rcParams.update({'font.size': fs,'xtick.labelsize': fs, 'ytick.labelsize': fs})
    fig,axh = plt.subplots(figsize=(8,8),nrows=4,ncols=3,subplot_kw = dict(projection=proj),
                           sharex=True,sharey=True,constrained_layout=False)
    plt.subplots_adjust(wspace=0.01, hspace=0.01) 
    ax = axh.flatten()
    cmap='Spectral_r'
    vmax=36
    vmin=31.6
    levels=np.arange(31.6,36.2,0.2)

    ####### oisss
    ds = cut2fullyears(ds).groupby('time.month').mean('time')
#     ds = cut2fullyears(ds).sel(time=slice(*time_overlap)).groupby('time.month').mean()
    ds = ds.transpose('month','lat','lon')
    for i in range(12):
        gl=fig_map(ax[i],plot_region)#[-74,-67,37,43]
        cc1 = ds[i,::].plot.contourf(ax=ax[i],levels=levels,transform=ccrs.PlateCarree(),add_colorbar=False,
                             vmin=vmin,vmax=vmax,cmap=cmap,
                             extend='both')
        ds[i,::].plot.contour(ax=ax[i],levels=[34],transform=ccrs.PlateCarree(),colors='mediumpurple',linewidths=0.5)
        cc=ax[i].contour(bathy.x,bathy.y,bathy.z*(-1),levels=[100,1000],
                         transform=ccrs.PlateCarree(),colors='k',alpha=0.5,linewidths=0.5,
                         linestyles='dashed')
        t=ax[i].text(0.05, 0.8, month_converter(ds[i,::].month.values), va='center',
                   fontsize = 10,transform = ax[i].transAxes,backgroundcolor='w')
        t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='white'))
        manual_locations = [(-74, 38.5)]
        ax[i].clabel(cc,fontsize=8,fmt='%1dm',inline=True, manual=manual_locations)
        # add Pioneer positions (inshore moorings)
        ax[i].plot(-70.8,40.4,marker='.',markersize=10, color='orchid',transform=ccrs.PlateCarree(),markeredgecolor='k',
                 label='Pioneer Inshore')
        if i in np.arange(0,9):
            gl.xlabels_bottom = False
        if i in [1,2,4,5,7,8,10,11]:
            gl.ylabels_left = False
            
        ## Gulf Stream time-mean path
        gs_lat_seas = cut2fullyears(da.sel(time=slice('2010-01','2017-12'))).mean('time') # data until 2018-05
        gs_lat_seas.plot(color='k',ax=ax[i],transform=ccrs.PlateCarree(),linewidth=0.5);
        ax[i].set_title('')

    # misc and save figure   
    cb = fig.colorbar(cc1,ax=axh[:],shrink=0.5,label='sea surface salinity [psu]',pad=0.02) 
    cb.ax.minorticks_off()
    # finished_plot(fig,f'./plots/manuscript/{name}_mean_seasonal_cycle_overlap_time_2016_2022_V2_incl_GoM_mean_GS_path_rev1.png',dpi=300)

##########################################
#%%% SI: Maps yearly seasonal cycle SMAP of overlapping perio
##########################################
font_for_pres()
plt.close('all')
for name,i in zip(['smap'],range(3)):#'smap',
#     plt.close('all')
    ## loop through years
    for year in np.arange(2015,2023):
        ds = sss_monthly[name].sel(time=slice(f'{year}-01',f'{year}-12'))
        fs=11
        plt.rcParams.update({'font.size': fs,'xtick.labelsize': fs, 'ytick.labelsize': fs})
        fig,axh = plt.subplots(figsize=(8,8),nrows=4,ncols=3,subplot_kw = dict(projection=proj),
                               sharex=True,sharey=True,constrained_layout=False)
        plt.subplots_adjust(wspace=0.01, hspace=0.01) 
        ax = axh.flatten()
        cmap='Spectral_r'
        vmax=36
        vmin=31.6
        levels=np.arange(31.6,36.2,0.2)
    
        ## loop over months
        ds = ds.transpose('time','lat','lon')
        for i in range(12):
            gl=fig_map(ax[i],plot_region)#[-74,-67,37,43]
            cc1 = ds[i,::].plot.contourf(ax=ax[i],levels=levels,transform=ccrs.PlateCarree(),add_colorbar=False,
                                 vmin=vmin,vmax=vmax,cmap=cmap,
                                 extend='both')
            ds[i,::].plot.contour(ax=ax[i],levels=[34],transform=ccrs.PlateCarree(),colors='mediumpurple',linewidths=0.5)
            cc=ax[i].contour(bathy.x,bathy.y,bathy.z*(-1),levels=[100],
                             transform=ccrs.PlateCarree(),colors='k',alpha=0.5,linewidths=0.5,
                             linestyles='dashed')
            t=ax[i].text(0.05, 0.8, month_converter(i+1), va='center',
                       fontsize = 10,transform = ax[i].transAxes,backgroundcolor='w')
            t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='white'))
            manual_locations = [(-74, 38.5)]
            ax[i].clabel(cc,fontsize=8,fmt='%1dm',inline=True, manual=manual_locations)
            # add Pioneer positions (inshore moorings)
            ax[i].plot(-70.8,40.4,marker='.',markersize=10, color='orchid',transform=ccrs.PlateCarree(),markeredgecolor='k',
                     label='Pioneer Inshore')
            if i in np.arange(0,9):
                gl.xlabels_bottom = False
            if i in [1,2,4,5,7,8,10,11]:
                gl.ylabels_left = False
                
            ## Gulf Stream path
            if year<2018:
                gs_lat_seas = da.sel(time=slice(f'{year}-01',f'{year}-12'))
        #         gs_lat_seas = cut2fullyears(da).groupby('time.month').mean('time')
                gs_lat_seas.isel(time=i).plot(color='k',ax=ax[i],transform=ccrs.PlateCarree(),linewidth=1);
            ax[i].set_title('')
            
            ## add wcr trajectories
            time = f'{year}-{i+1}'
            ax[i].plot(wcr.sel(time=time).Lon,wcr.sel(time=time).Lat,transform=ccrs.PlateCarree(),color='k',marker='.',linestyle='None')
    
        ## add colorbar   
        cb = fig.colorbar(cc1,ax=axh[:],shrink=0.5,label='sea surface salinity',pad=0.02) 
        cb.ax.minorticks_off()
        
        ## save figure
        # finished_plot(fig,f'./plots/satellite/SSS/{name}_{year}_seasonal_cycle_incl_GoM_and_WCRs.jpg',dpi=300)


