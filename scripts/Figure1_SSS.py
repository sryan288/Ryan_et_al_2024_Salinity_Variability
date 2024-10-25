#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 09:40:05 2024

@author: sryan

----------------------
FIGURE 1
Maps showing standard deviation of SSS and trend over full timeperiod
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
from utils.plot_utils import finished_plot,plot_map_NWA
from utils.xarray_utilities import cut2fullyears


def plot_map_NWA(ax,extent,c=np.arange(0,6000,500),plotbathy=None,gulfstream=False):
    """
    Produces map with coastlines using cartopy.
    
    INPUT:
    ax       : axis handle 
    extent   : [lon1,lon2,lat1,lat2]
    
    OUTPUT:
    gl       : grid handle
    
    """
    import scipy.io as sc
    import matplotlib.ticker as mticker

    
    ## load bathymetry
    bathy = xr.open_dataset('/vast/clidex/data/bathymetry/ETOPO2v2/ETOPO2v2c_f4.nc').sel(x=slice(-85,-40),y=slice(25,55)).load()
    # smooth bathy
    bathy = bathy.rolling(x=5,center=True).mean().rolling(y=5,center=True).mean()
    
    # projection
    proj = ccrs.PlateCarree()
    
    ## plot bathymetry
    if plotbathy=='contourf':
        ax.contourf(bathy.x,bathy.y,(bathy.z*(-1)),levels=c,cmap=plt.get_cmap('Blues',len(c)),transform=ccrs.PlateCarree())
    elif plotbathy=='contour':
        ax.contour(bathy.x,bathy.y,(bathy.z*(-1)),levels=c,colors='k',
                   transform=ccrs.PlateCarree(),
                   linewidths=0.5,alpha=0.5)
        
    ## add mean Gulf Stream path
    if gulfstream:
        # gs = sc.loadmat('/mnt/data/GS_monthly_yearly_two_sat_contours_1993_to_2018.mat')
        # gs_time = gs['time_monthly']
        # gs_lon = gs['lon_cell'][0,0][0,:]
        # gs_lat = gs['lat_cell'][0,0][0,:]
        # plt.plot(gs_lon,gs_lat,transform=ccrs.PlateCarree(),color='k',linewidth=1)
        gs_lat_seas = cut2fullyears(da.sel(time=slice('2010-01','2021-12'))).mean('time')
        gs_lat_seas.plot(color='k',ax=ax,transform=ccrs.PlateCarree(),linewidth=1);

    ## formatting
    ax.set_extent(extent)
#     ax.coastlines(resolution='50m',color='gray')
    gl = ax.gridlines(crs=proj,draw_labels=True,linestyle='-',alpha=0.1,
                      xlocs=range(-120,40,2),ylocs=range(0,80,1),x_inline=False)
    ax.add_feature(cartopy.feature.GSHHSFeature(scale='intermediate'), facecolor='lightgray',edgecolor='gray')
    gl.xlocator = mticker.FixedLocator(np.arange(-75,-50,5))
    gl.ylocator = mticker.FixedLocator(np.arange(30,50,5))
    gl.rotate_labels = False
    
    # gl.ylabels_left= False
    # gl.xlabels_bottom = False
    return gl

######
# function to derive trend and significance
#######

def trend_significance(data,alpha,timevar='time_counter'):
    import matplotlib.dates as dates
    import scipy.stats as stats
    # data has to have time,lon,lat dimensions, i.e. 3D.
    # Won't work with additional depth dimension
    # need to stack data again to have 1D array, which is required by linregress function
    stacked = data.fillna(0).stack(loc=('lon','lat'))
    timenum=dates.date2num(stacked[timevar].to_index().to_pydatetime()) #matlab time

    # initialize fields
    p_value=[]
    lintrend=[]
    trend = []
    for i in np.arange(0,stacked.shape[1]):
        # do linear regression
        slope,intercept,rval,pval,std_err = stats.linregress(timenum,stacked[:,i])
        p_value.append(pval)
        # derive linear trend
        dummy = slope*timenum+intercept
        trend.append(dummy)
        lintrend.append(dummy[-1]-dummy[0])
    # add values back into xarray & unstack
    foo = xr.DataArray(np.array(trend).transpose(),
                   coords=[stacked[timevar], stacked['loc']],
                   dims=['time', 'loc'],name=stacked.name)
    trend = foo.unstack('loc')
    foo = xr.DataArray(np.array(lintrend),
                  coords=[stacked['loc']],
                  dims=['loc'],name=stacked.name)
    lintrend = foo.unstack('loc')
    foo = xr.DataArray(np.array(p_value),
                  coords=[stacked['loc']],
                  dims=['loc'],name=stacked.name)
    pval = foo.unstack('loc')
    # now derive mask for desired significance value alpha
    mask = pval.values-alpha
    mask[mask>0]=np.nan
    mask[mask<0]=1
    
    # trend is full size trend, that can be subtracted from original field
    # lintrend is an integer, which gives trend over whole time period
    return mask,trend,lintrend



#%% Load data 
# =============================================================================

### salinity
product = 'smos'
## sea surface salinity data preprocesses for three products
sss_monthly = np.load('./data/sss_monthly.npy',allow_pickle=True)[product]
sss_daily = np.load('./data/sss_daily.npy',allow_pickle=True)[product]

## cut to time period where data exists 
sss_monthly = sss_monthly.where(np.isfinite(sss_monthly),drop=True)
sss_daily = sss_daily.where(np.isfinite(sss_daily),drop=True)

## if year at beginning or end not full, it will be cut off to not spur data
# sss_monthly = cut2fullyears(sss_monthly)
# sss_daily = cut2fullyears(sss_daily)

# derive annuaml mean for trend
sss_am = sss_daily.groupby('time.year').mean(skipna=True)
sss_am['year'] = (('year'),pd.date_range(start=f'{sss_daily.time.dt.year[0].values}-01-01',
                                         end=f'{sss_daily.time.dt.year[-1].values}-01-01',
                                         freq='YS').values)

### derive trend and std
#1) STD
monthly_std = sss_monthly.std('time',skipna=True)
daily_std = sss_daily.std('time',skipna=True)

# 2) trend
mask,trend,lintrend = trend_significance(sss_am,alpha=0.1,timevar='year')
lintrend = lintrend.where(np.isfinite(sss_am.isel(year=1)))

#%% plot Map

#### Salinity standard deviation
plt.close('all')
font_medium()
fig = plt.figure(figsize=(10*cm,6*cm),constrained_layout=True)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],projection=proj)
gl=plot_map_NWA(ax, [-78,-57,36,46],c=[1000,4000],plotbathy='contour',gulfstream=False)
gl.ylabels_right = False
gl.xlabels_bottom = False
# add Gulf Stream
gs_lat_seas = cut2fullyears(da.sel(time=slice('2010-01','2021-12'))).mean('time')
gs_lat_seas.plot(color='black',ax=ax,transform=ccrs.PlateCarree(),linewidth=1);
cc = monthly_std.plot(cmap='spectral_r',add_colorbar=False,vmin=0.3,vmax=1.2,
                      transform=ccrs.PlateCarree())
pp = plt.colorbar(cc,shrink=0.6,ticks=np.arange(0,2,0.2),extend='both',
             label='salinity standard deviation [psu]',orientation='horizontal')
pp.ax.minorticks_off()
ax.add_feature(cartopy.feature.RIVERS)
finished_plot(fig,f'./manuscript/{product}_monthly_std_gulfstream_mean_rev1.png')


###  Trend
plt.close('all')
font_medium()
fig = plt.figure(figsize=(10*cm,6*cm),constrained_layout=True)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8],projection=proj)
gl=plot_map_NWA(ax, [-78,-57,36,46],c=[1000,4000],plotbathy='contour',gulfstream=False)
gl.ylabels_left= False
gl.xlabels_bottom = False
# add Gulf Stream
gs_lat_seas = cut2fullyears(da.sel(time=slice('2010-01','2021-12'))).mean('time')
gs_lat_seas.plot(color='black',ax=ax,transform=ccrs.PlateCarree(),linewidth=1);

import proplot
cmap1 = 'BrBG_r'#proplot.Colormap('DryWet_r')
cc=lintrend.T.plot(transform=ccrs.PlateCarree(),
                            vmin=-0.8,vmax=0.8,cmap=cmap1,zorder=0,add_colorbar=False)
pp=plt.colorbar(cc,shrink=0.6,ticks=np.arange(-0.8,1,0.4),extend='both',
                label=f'salinity change {sss_daily.time.dt.year[0].values}-{sss_daily.time.dt.year[-1].values} [psu]',orientation='horizontal')
pp.ax.minorticks_off()
ax.add_feature(cartopy.feature.RIVERS)
switched_data = np.where(np.isnan(mask), 1, np.nan)
ax.contourf(lintrend.lon,lintrend.lat,switched_data.T,colors='None',level=[1],hatches=['....'],transform=ccrs.PlateCarree())
# ax.contourf(lintrend.lon,lintrend.lat,mask.T,
#             levels=[0,1],colors='None',hatches=['....'],transform=ccrs.PlateCarree())
ax.set_title('')
finished_plot(fig,f'./manuscript/{product}_am_trend_gulfstream_mean_rev1_alpha0.1.png')


# cc = daily_std.plot(cmap='spectral_r',add_colorbar=False,vmin=0.3,vmax=1.2,
#                       transform=ccrs.PlateCarree())
# pp = plt.colorbar(cc,shrink=0.8,ticks=np.arange(0,2,0.2),extend='both',
#              label='salinity standard deviation')
# pp.ax.minorticks_off()
# ax.add_feature(cartopy.feature.RIVERS)
# finished_plot(fig,'./manuscript/smos_daily_std_2010_2022.jpg')

