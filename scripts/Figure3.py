#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 13:25:46 2024

@author: sryan

----------------------
FIGURE 2
Mean seasonal cycle and std of SSS and CP data at pioneer locations
Plus temperature and density at pioner
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

# =============================================================================
#%% Load data
# =============================================================================

## sea surface salinity data preprocesses for three products
sss_monthly = np.load('./data/sss_monthly.npy',allow_pickle=True)

## load pioneer data prepared by prepare_pioneer.ipynb
pioneer_surf = np.load('/mnt/data/OOI_Pioneer/moor_surf_xr.npy',allow_pickle=True)
pioneer_prof = np.load('/mnt/data/OOI_Pioneer/moor_prof_xr.npy',allow_pickle=True)
surf_lp30_seas = np.load('./data/pioneer_TSrho_surf_lp30_seas_fullyears.npy',allow_pickle=True)

## overlaopping time of all datasets
time_overlap = ['2016-01-01','2022-12-31']

# =============================================================================
#%% Plot salinity
# =============================================================================
plt.close('all')
font_medium()
# plt.rc('figure', titlesize=8, titleweight='bold', autolayout=False) #,
cols = plt.get_cmap('tab10')
fig,axf = plt.subplots(figsize=(4,3.5),constrained_layout=True,nrows=2,sharex=True,
                       gridspec_kw={'height_ratios': [1.7, 1]}) # full time period

ax=axf[0]
#
# loop over SSS datasets and plot mean seasonal cycle at inshore position
for name,i,label in zip(['oisss','smap','smos'],range(3),['OISSS','SMAP','SMOS']):
    ds = sss_monthly[name].copy(deep=True)
    ds = ds.sel(lon=pioneer_loc['inshore'][1],
                lat=pioneer_loc['inshore'][0],
                method='nearest')
    ds_seas = cut2fullyears(ds).groupby('time.month').mean()
    ds_seas_std = cut2fullyears(ds).groupby('time.month').std()
    ds_seas_overlap = cut2fullyears(ds.sel(time=slice(*time_overlap))).groupby('time.month').mean()
    ds_seas_overlap_std = cut2fullyears(ds.sel(time=slice(*time_overlap))).groupby('time.month').std()
    
    # plot
    ds_seas_overlap.plot(ax=ax,color=cols(i),label=label,linewidth=1,zorder=10)
    ax.fill_between(ds_seas_overlap_std.month,
                     ds_seas_overlap-ds_seas_overlap_std,
                     ds_seas_overlap+ds_seas_overlap_std,
                     alpha=0.1,color=cols(i),zorder=0)
 

#
# add pioneer data    
surf_lp30_seas['issm']['sal'].plot(ax=ax,label='Pioneer Inshore',color='k',linewidth=0.5)
# surf_lp30_seas['cnsm']['sal'].plot(ax=ax,label='Pioneer Central',color='gray',linewidth=0.5)
surf_lp30_seas['ossm']['sal'].plot(ax=ax,label='Pioneer Offshore',color='lightgray',linewidth=0.5)
# plt.legend(['Pioneer Inshore','Pioneer Central','Pioneer Offshore'],
#           loc='lower right',fontsize=6,bbox_to_anchor=(0.9, 0))
# plt.gca().add_artist(legend1);
legend1 = ax.legend(loc='upper left',fontsize=6,bbox_to_anchor=(0.1, 1),ncol=2);

#
# misc
ax.set_ylabel('salinity [psu]')
ax.grid(False)
ax.minorticks_off()
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
ax.tick_params(axis='x', which='major',length=3)
ax.tick_params(axis='y', which='major',length=3)
ticks_month(ax)
ax.spines['bottom'].set_visible(False)


# finished_plot(fig,'./manuscript/sss_seasonal_cycle_pioneer_loc_full.png')


#################tempe and dens
# plt.close('all')
# fig,ax = plt.subplots(figsize=(4,2.3),constrained_layout=False) # full time period
# add pioneer data 
ax1=axf[1]   
surf_lp30_seas['issm']['temp'].plot(ax=ax1,label='Pioneer Inshore',color='k',linewidth=0.5)
# surf_lp30_seas['cnsm']['temp'].plot(ax=ax1,label='Pioneer Central',color='gray',linewidth=0.5)
surf_lp30_seas['ossm']['temp'].plot(ax=ax1,label='Pioneer Offshore',color='lightgray',linewidth=0.5)
#
# misc
# ax.legend(ncol=2,loc='upper left',fontsize=6);
ax1.set_ylabel('temperature [$^\circ$C]',labelpad=8)
ax1.grid(False)
ticks_month(ax1)
ax1.minorticks_off()
ax1.tick_params(axis='x', which='major',length=3)
ax1.tick_params(axis='y', which='major', direction='out', length=3)

ax2 = ax1.twinx()
# fig,ax = plt.subplots(figsize=(3.5,2),constrained_layout=True) # full time period
# add pioneer data    
surf_lp30_seas['issm']['dens'].plot(ax=ax2,label='Pioneer Inshore',color='k',linewidth=0.5,linestyle='dashed')
# surf_lp30_seas['cnsm']['dens'].plot(ax=ax2,label='Pioneer Central',color='gray',linewidth=0.5,linestyle='dashed')
surf_lp30_seas['ossm']['dens'].plot(ax=ax2,label='Pioneer Offshore',color='lightgray',linewidth=0.5,linestyle='dashed')
#
# misc
# ax.legend(ncol=2,loc='upper left',fontsize=6);
ax2.set_ylabel('density [kg/$m^3$]')
ax2.grid(False)
ax2.set_yticks(np.arange(1022,1027))
ax2.minorticks_off()
ax2.tick_params(axis='y', which='both', direction='out', length=3)
ax2.spines['top'].set_visible(False)
ax1.spines['top'].set_visible(False)


finished_plot(fig,'./plots/manuscript/Fig3_seasonal_cycle_pioneer_loc_sal_temp_dens_rev1.png')

# =============================================================================
#%% Plot seasonal data in alpha beta space 
# =============================================================================

# derive thermal expansion and haline contraction coefficient for TS plots
T = np.arange(-2,30,0.5)
S = np.arange(10,38)
Tm,Sm = np.meshgrid(T,S)
rho,alpha,beta = gsw.rho_alpha_beta(Sm,Tm,0)
pden = sw.pden(Sm,Tm,0,0)


#%%% plot alpha
# plot TS
plt.close('all')
font_medium()
fig,ax = plt.subplots(figsize=(2.5,2),constrained_layout=True)
cc=ax.contourf(Sm,Tm,alpha/1e-04,levels=np.arange(0,4,0.25),cmap='fire')
pp=plt.colorbar(cc,shrink=0.8,label=r'thermal expansion coefficient $\alpha$')
pp.ax.minorticks_off()
ax.grid(False)
ax.set_ylabel('temperature [$^\circ$C]')
ax.set_xlabel('salinity [psu]')
# ax.set_title(r'Thermal expansion coefficient $\alpha$')
ax.set_xlim(29,36.5)
ax.set_ylim(2.5,28)
# cc=ax.contour(Sm,Tm,alpha/1e-04,levels=[0],colors='w',linewidths=0.5)
cc=ax.contour(Sm,Tm,rho,levels=np.arange(1015,1030,1),colors='k',linestyles='dashed',linewidths=0.2)
plt.clabel(cc,levels=[1022,1024,1026],manual=True)

# add pioneer data
size = 8
cols = plt.get_cmap('gray',12)

for i in range(12):
    plt.plot(surf_lp30_seas['issm']['sal'].isel(month=i),
             surf_lp30_seas['issm']['temp'].isel(month=i),
             marker='*',markersize=size,linewidth=1,color=cols(i),zorder=100)

# finished_plot(fig,'./manuscript/alpha_TS_pioneer_inshore_seas_rev1.jpg')

#%%% beta

fig,ax = plt.subplots(figsize=(2.5,2),constrained_layout=True)
cc=ax.contourf(Sm,Tm,beta/1e-04,levels=np.arange(7,7.85,0.05),cmap='haline')
pp=plt.colorbar(cc,shrink=0.8,label=r'haline contraction coefficient $\beta $')
# pp.set_title('*$10^{-4}$')
ax.set_xlim(29,36.5)
ax.set_ylim(2.5,28)
cc=ax.contour(Sm,Tm,rho,levels=np.arange(1015,1030,1),colors='k',linestyles='dashed',linewidths=0.2)
plt.clabel(cc,levels=[1022,1024,1026],manual=True)
pp.ax.minorticks_off()
ax.grid(False)
ax.set_ylabel(r'temperature [$^\circ$C]')
ax.set_xlabel(r'salinity [psu]')
# ax.set_title(r'Haline contraction coefficient $\beta $')


# add pioneer data
size = 8
cols = plt.get_cmap('gray',12)
for i in range(12):
    plt.plot(surf_lp30_seas['issm']['sal'].isel(month=i),
             surf_lp30_seas['issm']['temp'].isel(month=i),
             marker='*',markersize=size,linewidth=1,color=cols(i))

# 
# finished_plot(fig,'./manuscript/beta_TS_pioneer_inshore_seas_rev1.jpg')

