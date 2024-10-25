#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 13:33:24 2023

@author: sryan
"""

import scipy.io
import os
# go to working directory
workingpath = '/home/sryan/python/NASA_salinity/p2_salinity_variability/'
os.chdir(workingpath)
filepath = workingpath+'./'
printpath = filepath+'oisst_trend.py'

STARTUP_2023_coldfronts_mhws = "./utils/2023_salinity_variability_startup.py"
exec(open(STARTUP_2023_coldfronts_mhws).read())


# my own packages
import sys
# sys.path.append("/home/sryan/python/") # go to parent dir
from utils.plot_utils import finished_plot,datetime_tick,plot_map_NWA
from utils.xarray_utilities import deseason
plt.close('all')

#%% Functions
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.dates as dates
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry as sgeom
import proplot

def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.
    
    """
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
              'right': [(maxx, miny), (maxx, maxy)],
              'bottom': [(minx, miny), (maxx, miny)],
              'top': [(minx, maxy), (maxx, maxy)],}
    return sgeom.LineString(points[side])


def lambert_xticks(ax, ticks):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
    te = lambda xy: xy[0]
    lc = lambda t, n, b: np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T
    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels])
    

def lambert_yticks(ax, ticks):
    """Draw ricks on the left y-axis of a Lamber Conformal projection."""
    te = lambda xy: xy[1]
    lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels])

def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """Get the tick locations and labels for an axis of a Lambert Conformal projection."""
    outline_patch = sgeom.LineString(ax.outline_patch.get_path().vertices.tolist())
    axis = find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(ccrs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(ccrs.Geodetic(), xy[:, 0], xy[:, 1])
        xyt = proj_xyz[..., :2]
        ls = sgeom.LineString(xyt.tolist())
        locs = axis.intersection(ls)
        if not locs:
            tick = [None]
        else:
            tick = tick_extractor(locs.xy)
        _ticks.append(tick[0])
    # Remove ticks that aren't visible:    
    ticklabels = copy(ticks)
    while True:
        try:
            index = _ticks.index(None)
        except ValueError:
            break
        _ticks.pop(index)
        ticklabels.pop(index)
    return _ticks, ticklabels

def detrend_xarray(da):
    dt = xr.apply_ufunc(
                linear_detrend,
                da,
                input_core_dims=[['time']],
                output_core_dims=[['time']],
                vectorize=True,
                output_dtypes=[da.dtype],
                dask="parallelized",
            )
    dt = dt.transpose('time','lat','lon')
    return dt

def trend_xarray(da,timevar='time'):
    dt = xr.apply_ufunc(
                linear_trend,
                da,
                input_core_dims=[[timevar]],
                output_core_dims=[[timevar]],
                vectorize=True,
                output_dtypes=[da.dtype],
                dask="parallelized",
            )
    dt = dt.transpose(timevar,'lat','lon')
    return dt

def linear_detrend(da):
    ds = da.copy()
    mask = ~np.isnan(ds)
    if mask.sum() == 0:
        return ds
    else:
        ds_masked = ds[mask]
        time_masked = np.arange(0,len(ds))[mask]
        coeff = np.polyfit(time_masked, ds_masked, 1)
        trend = np.polyval(coeff, time_masked)
        detrended = ds_masked - trend
        ds[mask] = detrended
        return ds

def linear_trend(da):
    ds = da.copy()
    mask = ~np.isnan(ds)
    if mask.sum() == 0:
        return ds
    else:
        ds_masked = ds[mask]
        time_masked = np.arange(0,len(ds))[mask]
        coeff = np.polyfit(time_masked, ds_masked, 1)
        trend = np.polyval(coeff, time_masked)
        # detrended = ds_masked - trend
        ds[mask] = trend
        return ds

def cut2fullyears(ds,dt='day'):
    years = ds['time.year']
    year1 = np.min(years)
    yearend = np.max(years)
    if dt=='day': lim = 365
    elif dt=='month':lim = 12
    if len(ds['time'].where(years==year1,drop=True))<lim:
        year1 = year1+1
    if len(ds['time'].where(years==yearend,drop=True))<lim:
        yearend=yearend-1
#         print(f'yearend: {yearend}')
    return ds.sel(time=slice(f'{year1.values}-01',f'{yearend.values}-12'))

import scipy.stats as stats
def trend_significance(data,alpha,timevar='time'):
    # data has to have time,lon,lat dimensions, i.e. 3D.
    # Won't work with additional depth dimension
    # need to stack data again to have 1D array, which is required by linregress function
    stacked = data.fillna(0).stack(loc=('lon','lat'))
    # need numerical time for regression
    if timevar=='time':
        timenum=dates.date2num(stacked[timevar].to_index().to_pydatetime()) #matlab time
    elif timevar=='year':
        timenum=stacked[timevar].values #matlab time

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


# =============================================================================
#%% Load data, manipulate
# =============================================================================

# OISST data
sst = xr.open_dataset('/vast/clidex/data/obs/SST/OISST/OISSTv2.1/oisst.monmean_1982.01-2023.07_90W40W20N60N.nc').sel(time=slice('1982-01','2022-12'))

# convert longitude values
sst['lon']=sst['lon']-360
# cut data to only include full years - important for trend analysis
sst_cut = cut2fullyears(sst,dt='month')
# remove seasonal cycle to derive anomalies
sst_deseason = deseason(sst_cut,timevar='time',refperiod=['1983','2012']).drop('month')
# create annual mean
sst_am = sst_cut.groupby('time.year').mean()


# =============================================================================
#%%% Derive trend and significance
# =============================================================================
mask,trend,lintrend = trend_significance(sst_deseason['sst'],alpha=0.01,timevar='time')
mask_am,trend_am,lintrend_am = trend_significance(sst_am['sst'],alpha=0.01,timevar='year')

# =============================================================================
#%% Plotting
# =============================================================================

# define colormap for plot
cmap1 = proplot.Colormap('fire',name='fraction')
plt.rcParams.update({'font.size': 10})

# =============================================================================
#%%% Timeseries temperature anomaly
# =============================================================================

fig,ax2 = plt.subplots(figsize=(2.5,1.5),constrained_layout=True)
# plt.rc('axes.spines', top= True, right=True )
dummy=sst_deseason.sel(lon=slice(-77,-50),lat=slice(35,45)).mean(('lat','lon'),skipna=True)['sst'].rolling(time=24,center=True).mean()
ax2.plot(dummy.time,dummy,color=cmap1(150),linewidth=0.5)
plt.minorticks_off()
ax2.xaxis.set_label_position('top') 
ax2.set_xticklabels(ax2.get_xticks(), rotation = 45)
ax2.xaxis.set_major_locator(mdates.YearLocator(10))
ax2.axhline(0,color='k',linestyle='solid',linewidth=0.5,zorder=0)

ax2.xaxis.set_major_formatter(dates.DateFormatter('%Y')) 
ax2.xaxis.tick_top()
ax2.set_ylabel('SST anomaly [\N{DEGREE SIGN}C]')

# finished_plot(fig,'temp_anomaly_74W50W_35N45N_large_font_V2.jpg')

# =============================================================================
#%%% Temp trend plot final
# =============================================================================
# font_for_pres()
# plt.rcParams.update({'font.size': 10})
font_medium()

plt.close('all')
fig,ax = plt.subplots(figsize=(7,4),sharex=True,sharey=True,
                          subplot_kw = dict(projection=proj),constrained_layout=False)
gl=plot_map_NWA(ax,[-77,-47,34,47],c=[1000,4000],plotbathy='contour',gulfstream=False)
cc=lintrend.T.plot.contourf(levels=np.arange(0,4.25,0.25),transform=ccrs.PlateCarree(),
                                                          cmap=cmap1,vmin=0,vmax=4,
                                                          zorder=0,add_colorbar=False)
plt.colorbar(cc,shrink=0.6,ticks=np.arange(0,5.25,1),label='temperature change \n 1983-2022 [\N{DEGREE SIGN}C]')
# add Gulf Stream
gs_lat_seas = cut2fullyears(da.sel(time=slice('2010-01','2021-12'))).mean('time')
gs_lat_seas.plot(color='black',ax=ax,transform=ccrs.PlateCarree(),linewidth=1);

#mask significance
switched_data = np.where(np.isnan(mask), 1, np.nan)
ax.contourf(lintrend.lon,lintrend.lat,switched_data.T,colors='None',level=[1],hatches=['....'],transform=ccrs.PlateCarree())

# add rectangle for time series
ax.plot(np.linspace(-75,-50,20),np.ones([20])*35,transform=ccrs.PlateCarree(),color='k',linestyle='dashed',linewidth=0.5)
ax.plot(np.linspace(-75,-50,20),np.ones([20])*45,transform=ccrs.PlateCarree(),color='k',linestyle='dashed',linewidth=0.5)
ax.plot(np.ones([20])*-50,np.linspace(35,45,20),transform=ccrs.PlateCarree(),color='k',linestyle='dashed',linewidth=0.5)
ax.plot(np.ones([20])*-75,np.linspace(35,45,20),transform=ccrs.PlateCarree(),color='k',linestyle='dashed',linewidth=0.5)

finished_plot(fig,'./plots/manuscript/map_sst_trend_am_significance_rev1_smaller_font.jpg')
# 