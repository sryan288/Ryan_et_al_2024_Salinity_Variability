#!/usr/bin/env python
# coding: utf-8

# # Do salinity variability analysis
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
# pioneer_prof = np.load('/mnt/data/OOI_Pioneer/moor_prof_xr.npy',allow_pickle=True)

# ## oisst
# oisst = xr.open_dataset('/vast/clidex/data/obs/SST/OISST/OISSTv2.1/oisst.day.mean_1982.01-2023.6_90W0W20N70N.nc')
# oisst['lon'] = oisst['lon']-360
#
## pioneer (data originally created further down in this script)
surf_lp30 = np.load('./data/pioneer_TSrho_surf_lp30.npy',allow_pickle=True)
surf_lp30_mm = np.load('./data/pioneer_TSrho_surf_lp30_mm.npy',allow_pickle=True)
surf_lp30_seas = np.load('./data/pioneer_TSrho_surf_lp30_seas_fullyears.npy',allow_pickle=True)

## SSS in ecoboxes (derived further down)
smap_eco = np.load('./data/smap_ecoboxes.npy',allow_pickle=True)
smos_eco = np.load('./data/smos_ecoboxes.npy',allow_pickle=True)
oisss_eco = np.load('./data/oisss_ecoboxes.npy',allow_pickle=True)

## wcr ring trajectories
# wcr = xr.open_dataset('/mnt/data/WCRs/wcrs_2011_2020_silver_et_al_2022.nc')

## river discharge data
# load discharge data
conv = 0.0283168 # convertion coefficient from ft^3 to m^3
path = '../data/compiled_discharge/'
discharge_all = xr.open_dataset(path + 'river_discharge_long.nc')*conv
discharge_MAB = xr.open_dataset(path + 'river_discharge_MAB_long.nc')*conv
discharge_GoM = xr.open_dataset(path + 'river_discharge_GoM_long.nc')*conv

discharge_cumsum = xr.open_dataset(path +'discharge_MAB_annual_cumsum_wateryear.nc')['__xarray_dataarray_variable__']*conv
discharge_cumsum_GoM = xr.open_dataset(path +'discharge_GoM_annual_cumsum_wateryear.nc')['__xarray_dataarray_variable__']*conv

        
# =============================================================================       
#%% Figure7a: Timeseries SSS and CP at Pioneer location
moor = 'inshore'

plt.close('all')
font_for_pres()
fig,ax = plt.subplots(figsize=(8,3))

# satellite data
sss_monthly['oisss'].sel(lon=pioneer_loc[moor][1],lat=pioneer_loc[moor][0],method='nearest').plot(label='OISSS',lw=1)
sss_monthly['smap'].sel(lon=pioneer_loc[moor][1],lat=pioneer_loc[moor][0],method='nearest').plot(lw=1,color='gold',label='SMAP')
sss_monthly['smos'].sel(lon=pioneer_loc[moor][1],lat=pioneer_loc[moor][0],method='nearest').plot(lw=1,color='green',label='SMOS')

# Pioneer data
surf_lp30_mm['issm']['sal'].plot(lw=0.5,color='k',label='Pioneer Inshore',alpha=0.5)
# surf_lp30_mm['cnsm']['sal'].plot(lw=0.5,color='gray',label='Central',alpha=0.5)
surf_lp30_mm['ossm']['sal'].plot(lw=0.5,color='lightgray',label='Pioneer Offshore',alpha=0.5)

## make pretty
plt.xlim(np.datetime64('2010-01'),np.datetime64('2022-12'))
plt.ylabel('sea surface salinity [psu]')
plt.xlabel('')
plt.title('')
plt.legend(ncol=2,fontsize=8,loc='best')
datetime_tick(ax)
plt.xticks(rotation=0,ha='center')
plt.grid(False)
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
# ax.spines['bottom'].set_visible(False)
ax.minorticks_off()

## save figure
finished_plot(fig,f'./manuscript/plots/Fig7a_sal_pioneer_all_oisss_smap5_smos.jpg',dpi=300)

# =============================================================================
#%% Figure7b: Discharge
# =============================================================================

## plot setup
plt.close('all')
font_for_pres()
plt.rcParams.update({'font.size': 10,'xtick.labelsize': 10, 'ytick.labelsize': 10})
fig,ax = plt.subplots(figsize=(8,2.5))

## discharge 
plt.plot(discharge_MAB.datetime,discharge_MAB['sum_discharge']/1e03,label='MAB')
# plt.plot(discharge_MAB.datetime,discharge_MAB['CONNECTICUT RIVER']*5,color='gray',linewidth=0.5)
plt.plot(discharge_GoM.datetime,discharge_GoM['sum_discharge']/1e03,label='GoM',color='gray',linewidth=0.5)
plt.xlim([np.datetime64('2010-01'),np.datetime64('2022-12')])
plt.ylim(0,11000/1e03)
datetime_tick(ax)
ax.grid(False)
ax.set_ylabel('discharge [m$^3$/s]')
ax.legend(ncol=2,fontsize=10)
ax.minorticks_off()

## cumulative discharge twin axis
ax2=ax.twinx()
cols = plt.get_cmap('tab10')
ax2.fill_between(discharge_cumsum['time'].values+np.timedelta64(30*9,'D'),discharge_cumsum/1e03,0,color=cols(0),zorder=0,alpha=0.1)
ax2.plot(discharge_cumsum['time'].values+np.timedelta64(30*9,'D'),discharge_cumsum/1e03,marker='*',linestyle='None',color=cols(0),zorder=0,alpha=0.3)
ax2.set_ylim(2.5e04/1e03,5.5e04/1e03)
ax2.set_ylabel('cumulative discharge [m$^3$/s]')
datetime_tick(ax2)
ax2.grid(False)
ax2.minorticks_off()

## fix time axis
import matplotlib.ticker
import matplotlib.dates as dt
dates = discharge_cumsum['time'].values+np.timedelta64(30*9,'D')
ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(dt.date2num(dates)))

## make pretty
ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

## save figure
finished_plot(fig,'./manuscript/plots/Fig7b_discharge_m3_2010_2022_rev1.png')


# =============================================================================
#%% Figure 7c: Hovmoeller plot SSS
# =============================================================================

# define lat & lon
lats = [38,41.5]
lons = [-71.3,-70.5]

#---------------------------------
# plot Hovmoeller
#---------------------------------
plt.close('all')
font_for_print()
def plot_hovmoeller(ds):
    plt.rcParams.update({'font.size': 10,'xtick.labelsize': 10, 'ytick.labelsize': 10})
    fig,ax = plt.subplots(figsize=(8,4))
    cc=ds.sel(lon=slice(*lons),
                   lat=slice(*lats)).mean('lon').T.plot.contourf(levels=np.arange(32,36,0.25),
                                                                 cmap = 'Spectral_r',
                                                                   add_colorbar=False)
    cb=plt.colorbar(cc,orientation='horizontal',shrink=0.4,pad=0.09)   
    cb.ax.tick_params(labelsize=9)
    cb.set_label('sea surface salinity [pau]',size=10)
    cb.ax.minorticks_off()                               
    cc=ds.sel(lon=slice(*lons),
                   lat=slice(*lats)).mean('lon').T.plot.contour(levels=[34,35],colors='k',linewidths=0.5,alpha=0.3)                                                            
    # cc=ds.sel(lon=slice(*lons),
    #                lat=slice(*lats)).mean('lon').T.plot.contour(levels=[34,35],colors='k')
    # plt.clabel(cc,fmt='%1d',labelspacing=100)
    ax.axhline(pioneer_loc['inshore'][0],color='k',linewidth=0.5,linestyle='dashed')
    ax.set_ylim(38.2,41.1)
    # ax.set_title(f'sss [lon={lons}]')
    ax.set_yticks([39,40,41])
    latlon_label(ax,coord='latitude',axis='y',fmt='%1d')
    plt.xticks(rotation=0,ha='center')
    ax.set_ylabel('')
    ax.set_xlim([np.datetime64('2010-01'),np.datetime64('2022-12')])
    datetime_tick(ax)
    ax.set_xlabel('')
    # plt.xticks(position='both')
    # plt.grid(False)
    return fig


# loop through datasets and plot hovmoeller
for name in ['smos','oisss','smap'][:1]:
    ds = sss_monthly[name].squeeze().sel(time=slice('2010-01','2022-12'))
    fig = plot_hovmoeller(ds)
    finished_plot(fig,f'./manuscript/plots/Fig7c_hovmoeller_71W_{name}_2010_2022.jpg')
# plot_hovmoeller(ds1)



# =============================================================================
#%%  Figure9: Extract SSS in ecoregions
# =============================================================================

# - OISST and OISSS are on the same grid, so the same mask can be applied  
# - have to interpolate on to SMAP mask

mask_names = ['Southern MAB','Northern MAB','Georges Bank','Western GoM',
         'Eastern GoM','Western SS','Eastern SS','Southern GSL','Northern GSL',
         'NFL Shelf','Northern NFL Shelf','Labrador Shelf1', 'Labrador Shelf2',
         'Labrador Shelf3','A','B','C','D','E','F','G','H']

################################################
#%%%interpolate masks on grid and extract regions
###############################################

def int_ecomask(lon,lat,mask):
    mask_eco = xr.open_dataset('/vast/clidex/data_old/ORCA/ecoregions_NEMO_masks.nc')
    mask_names = ['Southern MAB','Northern MAB','Georges Bank','Western GoM',
             'Eastern GoM','Western SS','Eastern SS','Southern GSL','Northern GSL',
             'NFL Shelf','Northern NFL Shelf','Labrador Shelf1', 'Labrador Shelf2',
             'Labrador Shelf3','A','B','C','D','E','F','G','H']
    [xm,ym] = np.meshgrid(lon,lat,indexing='ij')
    for i in range(22):
        dummy = griddata((mask_eco['nav_lon'].values.ravel(),mask_eco['nav_lat'].values.ravel()),
                        mask_eco['mask'][i,::].values.ravel(), (xm, ym),'nearest')
        dummy[dummy==1]=np.nan
        dummy[dummy==0]=1
        mask[mask_names[i]] = (('lon','lat'),dummy)
    return mask

# extract regions
def extract_regions(ds,mask):
    """
    ds has to have dimensions (time,lon,lat)
    """
    ds_eco = ds*mask
    for key in list(mask.keys())[:1]:
        dummy = mask[key].values
        dummy[dummy==0]=np.nan
        ds_eco[key] = dummy*ds
        
    return ds_eco 

#################################################################
def extract_eco(ds,lon,lat,mask):
    ##################
    # 1) interpolate mask onto grid
    mask_eco = xr.open_dataset('/vast/clidex/data_old/ORCA/ecoregions_NEMO_masks.nc')
    mask_names = ['Southern MAB','Northern MAB','Georges Bank','Western GoM',
             'Eastern GoM','Western SS','Eastern SS','Southern GSL','Northern GSL',
             'NFL Shelf','Northern NFL Shelf','Labrador Shelf1', 'Labrador Shelf2',
             'Labrador Shelf3','A','B','C','D','E','F','G','H']
    [xm,ym] = np.meshgrid(lon,lat,indexing='ij')
    for i in range(22):
        dummy = griddata((mask_eco['nav_lon'].values.ravel(),mask_eco['nav_lat'].values.ravel()),
                        mask_eco['mask'][i,::].values.ravel(), (xm, ym),'nearest')
        dummy[dummy==1]=np.nan
        dummy[dummy==0]=1
        mask[mask_names[i]] = (('lon','lat'),dummy)
        
    ##################
    # 2) extract regions
    ds_eco = ds*mask
    for key in list(mask.keys())[:1]:
        dummy = mask[key].values
        dummy[dummy==0]=np.nan
        ds_eco[key] = dummy*ds
       
    # return output
    return ds_eco 

################################################################################
# extract ecoregions in datasets
################################################################################
mask = xr.ones_like(sss_monthly['smap'][0,::]).to_dataset().drop('smap')
smap_eco = extract_eco(sss_monthly['smap'].transpose('time','lon','lat'),mask.lon,mask.lat,mask)
# smap_unc_eco = extract_eco(smap5['sss_smap_unc'].transpose('time','lon','lat'),smap5.lon,smap5.lat,mask)
smos_eco = extract_eco(sss_monthly['smos'].transpose('time','lon','lat'),mask.lon,mask.lat,mask)
# smap5_unc_eco = extract_eco(smap5['sss_smap_unc'].transpose('time','lon','lat'),smap5.lon,smap5.lat,mask)
oisss_eco = extract_eco(sss_monthly['oisss'].transpose('time','lon','lat'),mask.lon,mask.lat,mask)


# =============================================================================
#%%% Decompose using stats model
# =============================================================================
import statsmodels.api as sm  # https://www.statsmodels.org --> code by chat GPT

smos_eco_decompose = {}
smos_eco_decompose_trend = []
smos_eco_decompose_seas = []
smos_eco_decompose_resid = []


mask_names = ['Southern MAB','Northern MAB','Georges Bank','Western GoM',
         'Eastern GoM','Western SS','Eastern SS','Southern GSL','Northern GSL',
         'NFL Shelf','Northern NFL Shelf','Labrador Shelf1', 'Labrador Shelf2',
         'Labrador Shelf3','A','B','C','D','E','F','G','H']

for box,i in zip(mask_names,range(len(mask_names))):
    print(box)
    # define your data and change to pandas dataframe
    x = smos_eco[box].mean(('lat','lon')).time
    y = smos_eco[box].mean(('lat','lon'))
    wo = np.where(np.isfinite(y.values))[0]
    if len(wo)>24:
        data = pd.DataFrame({'Date': x[wo], 'Value': y[wo]})
        print(wo)
        # Decompose the time series data
        result = sm.tsa.seasonal_decompose(data['Value'], model='additive', period=12)  # Assuming a seasonal cycle of 12 months
        smos_eco_decompose[box] = result
        # extract trend
        if i ==0: smos_eco_decompose_trend = np.array(result.trend).T
        else: smos_eco_decompose_trend = np.vstack([smos_eco_decompose_trend,np.array(result.trend).T])
        # extract seasonal cycle
        if i ==0: smos_eco_decompose_seas = np.array(result.seasonal).T
        else: smos_eco_decompose_seas = np.vstack([smos_eco_decompose_seas,np.array(result.seasonal).T])
        # extract residuals
        if i ==0: smos_eco_decompose_res = np.array(result.resid).T
        else: smos_eco_decompose_res = np.vstack([smos_eco_decompose_res,np.array(result.resid).T])
        
        #plot
        # fig,ax = plt.subplots(nrows=4,figsize=(8, 6),sharex=True)
        # ax[0].plot(data['Date'], data['Value'], label='Original Data')
        # ax[0].legend(loc='best')
        # ax[1].plot(data['Date'], result.trend, label='Trend')
        # ax[1].legend(loc='best')
        # ax[2].plot(data['Date'], result.seasonal, label='Seasonal')
        # ax[2].legend(loc='best')
        # ax[3].plot(data['Date'], result.resid, label='Residuals')
        # ax[3].legend(loc='best')
        # ax[3].axhline(0,color='k',linestyle='dashed',linewidth=0.5)
        # fig.suptitle(box)
        # ax[3].set_xticks(pd.date_range(start='2010-01-01',end='2023-12-31',freq='YS').strftime('%Y'))
        # ax[3].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        # ax.xaxis.set_major_locator(locator) ## calling the locator for the x-axis
        # finished_plot(fig,f'./manuscript/plots/satellite/SSS/additive_model/smos_additive_model_box_{box}.png')
        
    else: print(box + 'not enough values')
        

        
#%%%% create xarray

da_trend = xr.Dataset(
    data_vars=dict(smos_eco_trend=(["time", "index"],smos_eco_decompose_trend.T)),
    coords=dict(
        time=(["time"], smos_eco['A'].time.values),
        box=(["index"], list(smos_eco_decompose.keys())),
    ),
    attrs=dict(
        description="Applied additive model from statsmodel https://www.statsmodels.org. This is the trend component in ecoboxes"
    ),
)

da_seas = xr.Dataset(
    data_vars=dict(smos_eco_seas=(["time", "index"],smos_eco_decompose_seas.T)),
    coords=dict(
        time=(["time"], smos_eco['A'].time.values),
        box=(["index"], list(smos_eco_decompose.keys())),
    ),
    attrs=dict(
        description="Applied additive model from statsmodel https://www.statsmodels.org. This is the seasonal component in ecoboxes"
    ),
)

da_res = xr.Dataset(
    data_vars=dict(smos_eco_res=(["time", "index"],smos_eco_decompose_res.T)),
    coords=dict(
        time=(["time"], smos_eco['A'].time.values),
        box=(["index"], list(smos_eco_decompose.keys())),
    ),
    attrs=dict(
        description="Applied additive model from statsmodel https://www.statsmodels.org. This is the residual component in ecoboxes"
    ),
)


#%%%% Figure 9: Decomposed SMOS LF signal

########################################
# Trend
#########
# convert to pandas datafram to plot heatmap
test = (da_trend['smos_eco_trend']-np.nanmean(da_trend['smos_eco_trend'],axis=0)).to_pandas()
test.columns = list(smos_eco_decompose.keys())

font_for_pres()
import seaborn as sns

fig = plt.figure(figsize=(8,4))

ax1 = fig.add_axes([0.15, 0.4, 0.75, 0.5])
ax2 = fig.add_axes([0.15, 0.1, 0.75, 0.25])
cbar_ax = fig.add_axes([.9, .3, .02, .4])
# shelf boxes
ax=sns.heatmap(test.iloc[:,0:7].T,ax=ax1, annot=False,cmap=plt.get_cmap('curl',len(np.arange(-0.7,0.75,0.05))),
                 fmt='1.0f',cbar_ax=cbar_ax,cbar_kws={'label':'salinity anomaly [psu]','ticks':np.arange(-0.6,0.8,0.2)},vmin=-0.7,vmax=0.7)
ax.set_xticks(np.arange(0,154,12))
ticklabels = [test.index[int(tick)].strftime('%Y') for tick in ax.get_xticks()]
ax.set_xticklabels('',rotation=0);
ax.invert_yaxis()
ax.set_xlabel('')
# slope boxes
sns.heatmap(test.iloc[:,9:].T, ax=ax2,annot=False,cmap=plt.get_cmap('curl',len(np.arange(-0.7,0.75,0.05))),
                 fmt='1.0f',cbar=False,cbar_kws={'label':'salinity anomaly [psu]','ticks':np.arange(-0.6,0.8,0.2)},vmin=-0.7,vmax=0.7)
ax2.set_xticks(np.arange(0,154,12))
ticklabels = [test.index[int(tick)].strftime('%Y') for tick in ax.get_xticks()]
ax2.set_xticklabels(ticklabels,rotation=0);
ax2.invert_yaxis()
ax2.set_xlabel('')
# save figure
# finished_plot(fig,'./manuscript/plots/manuscript/Figure9_smos_additive_model_trend anomaly_ecoboxes.png')


