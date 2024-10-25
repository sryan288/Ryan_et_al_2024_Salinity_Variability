#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create overview maps for publication


Created on Mon Dec 18 15:15:05 2023

@author: sryan
"""


import os
# go to working directory
workingpath = '/home/sryan/python/NASA_salinity/p2_salinity_variability/'
os.chdir(workingpath)
filepath = workingpath+'/scripts/'
printpath = filepath+'woa_box_analysis.py'

STARTUP_2023_coldfronts_mhws = "./utils/2023_salinity_variability_startup.py"
exec(open(STARTUP_2023_coldfronts_mhws).read())

# my own packages
import sys
sys.path.append("./") # go to parent dir
from utils.plot_utils import finished_plot,plot_map_NWA,datetime_tick,print_path,ts_append


#%% Load data
# woa = xr.open_dataset('/mnt/data/WOA/woa_xarray.nc')
woa_box = np.load('/mnt/data/WOA/woa_boxes.npy',allow_pickle=True)
pioneer = np.load('/mnt/data/OOI_Pioneer/pioneer_Lukas_dict.npy',allow_pickle=True)

#### Polygons
# extract positions of 20m and 200m isobath
plt.figure()
cc=bathy.sel(x=slice(-76,-68),y=slice(36,42)).z.plot.contour(levels=[-200])
x200 = np.squeeze(cc.allsegs[0][0][:,0])
y200 = cc.allsegs[0][0][:,1]
cc=bathy.sel(x=slice(-76,-68),y=slice(36,42)).z.plot.contour(levels=[-20])
x0 = cc.allsegs[0][0][:,0]
y0 = cc.allsegs[0][0][:,1]
plt.plot(x0,y0,color='r')
plt.plot(x200,y200,color='k')
plt.close('all') 

# =============================================================================
# 1) create Polygon
# =============================================================================
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
# construct polygon
def construct_polygon(pos,posvar):  
    """
    Function to create a polygon on the shelf bound by the 20m and 200m isobath
    
    
    Parameters
    ----------
    pos : vector with coordinates for 4 corner points [SE,NE,SW,NW], [p1,p2,p3,p4]
    var : variables for 4 corners, e.g. x200 vs y200

    Returns
    -------
    xbox : Vector with polygon longitudes
    ybox : Vector with polygon latitudes

    """
    # offshore edge
    wo_p1p2 = np.squeeze(np.where((posvar[0]>pos[0]) & (posvar[1]<=pos[1])))
    yp1 = y200[wo_p1p2[0]]
    xp1 = x200[wo_p1p2[0]]
    xp1p2 = x200[wo_p1p2]
    yp1p2 = y200[wo_p1p2]
    yp2 = y200[wo_p1p2[-1]]
    xp2 = x200[wo_p1p2[-1]]
    # inshore edge
    wo_p3p4 = np.squeeze(np.where((posvar[2]<=pos[2]) & (posvar[3]>pos[3])))
    yp3 = y0[wo_p3p4[-1]]
    xp3 = x0[wo_p3p4[-1]]
    xp3p4 = x0[wo_p3p4]
    yp3p4 = y0[wo_p3p4]
    yp4 = y0[wo_p3p4[0]]
    xp4 = x0[wo_p3p4[0]]
    # final point
    xp5 = xp1
    yp5 = yp1
    
    
    xbox = []
    ybox = []
    for x,y in  zip([xp1p2,xp3,np.flipud(xp3p4),xp5],
                    [yp1p2,yp3,np.flipud(yp3p4),yp5]):
        xbox = np.append(xbox,x)
        ybox = np.append(ybox,y)
        
    return xbox,ybox

###### call function to create four boxes
xbox1,ybox1 = construct_polygon([36,37.2,37.6,36],[y200,y200,y0,y0])
xbox2,ybox2 = construct_polygon([37.2,38.5,39.3,37.6],[y200,y200,y0,y0])
xbox3,ybox3 = construct_polygon([38.5,39.7,40.8,39.3],[y200,y200,y0,y0])
xbox4,ybox4 = construct_polygon([39.7,-70,-70,40.8],[y200,x200,x0,y0])


#%% Data map
plt.close('all')

## simple bathymetry map
font_for_pres()
extent = [-79,-69,37,42]
fig,ax = plt.subplots(figsize=(4,3),
                      subplot_kw = dict(projection=proj_rot),
                      constrained_layout=True)
gl = plot_map_NWA(ax,extent,
                  gulfstream=False)
ax.add_feature(cartopy.feature.RIVERS)
gl.xlocator = mticker.FixedLocator(np.arange(-80,60,2.5))
gl.ylocator = mticker.FixedLocator(np.arange(35,50,2.5))

## add bathymetry
cc=ax.contourf(bathy.x,bathy.y,bathy.z/(-1000),levels=np.arange(0,6,0.5),cmap='deep',transform=ccrs.PlateCarree())
axc = fig.add_axes([0.17, 0.9, 0.25, 0.03])
cb=fig.colorbar(cc,cax=axc,orientation='horizontal',ticks=np.arange(0,6,1))
cb.ax.tick_params(labelsize=8)
cb.minorticks_off()
cb.set_label(label='depth [km]',size=8)

## add river stations
riverid = [1184000,1358000,1446500,1570500,1638500,2037500,2080500]
lat = [41.987319,42.752222,40.826389,40.254812,39.273583,37.563202,36.460000]
lon = np.array([72.605367,73.688889,75.082500,76.886085,77.543111,77.546931,77.633611])*[-1]
ax.plot(lon,lat,marker='.',linestyle='None',transform=ccrs.PlateCarree(),color='royalblue',markersize=10)


cols = plt.get_cmap('tab10')
## add woa datapoints
for box in woa_box.keys():
    ax.plot(woa_box[box].lon,woa_box[box].lat,transform=ccrs.PlateCarree(),
              marker='.',
              linestyle='None',
              markersize=0.5,
              zorder=1,color=cols(0))

## add boxes
ax.plot(xbox1,ybox1,color='k',transform=ccrs.PlateCarree(),linewidth=1,zorder=10)
ax.plot(xbox2,ybox2,color='k',transform=ccrs.PlateCarree(),linewidth=1,zorder=10)
ax.plot(xbox3,ybox3,color='k',transform=ccrs.PlateCarree(),linewidth=1,zorder=10)
ax.plot(xbox4,ybox4,color='k',transform=ccrs.PlateCarree(),linewidth=1,zorder=10)

# ## add Oleander
# lat_oleander = [40.55,32+17/60+40/60]
# lon_oleander = [-74.05,-(64+47/60+9/60)]
# ax.plot(lon_oleander,lat_oleander,color='orchid',transform=ccrs.PlateCarree(),linewidth=0.5,zorder=10)

## add Pioneer location --> do I need insert for separate moorings?
ax.plot(-70.8,40.3,marker='s',color='darkred',markersize=5,transform=ccrs.PlateCarree())
ax_insert = fig.add_axes([0.69, 0.15, 0.2, 0.2],projection=proj_rot)
gl = plot_map_NWA(ax_insert,[-71.5,-69.8,40,40.5],
                  gulfstream=False)
gl.xlocator = mticker.FixedLocator([])
gl.ylocator = mticker.FixedLocator([])
cc=ax_insert.contourf(bathy.x,bathy.y,bathy.z*(-1),levels=np.arange(0,1000,50),cmap='gray_r',transform=ccrs.PlateCarree(),extend='both')

for key in pioneer_loc.keys():
    ax_insert.plot(pioneer_loc[key][1],pioneer_loc[key][0],transform=ccrs.PlateCarree(),marker='.',markersize=3,color='k')
for key in ['inshore','centr_inshore','offshore']:#pioneer_loc.keys():
    ax_insert.plot(pioneer_loc[key][1],pioneer_loc[key][0],transform=ccrs.PlateCarree(),marker='.',markersize=5,color='darkred')


finished_plot(fig,'/home/sryan/python/NASA_salinity/p2_salinity_variability/plots/manuscript/data_map_rev1.png')