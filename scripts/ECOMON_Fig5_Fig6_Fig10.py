# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 21:30:57 2023

@author: svenj
"""

"""
This file reads the WOA matlab file provided by Adrienne (Paula Fratantoni)
and visualizes the data
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
from utils.datafun import matlab2datetime
from utils.xarray_utilities import deseason
import dask
import gsw
import calendar 



# file create by woa_mat2xarray.py
woa = xr.open_dataset('/mnt/data/WOA/woa_xarray.nc')
woa_box = np.load('/mnt/data/WOA/woa_boxes.npy',allow_pickle=True)



##################################################
#%% Figure 5 and Figure 6: Mean cross-shelf structure in each box
##################################################

#%%% groupby distance bins
# time average
mean_sec = {}
mean_sec_count = {}
interv = 10 # km
for name in ['box1','box2','box3','box4']:
    mean_sec[name] = woa_box[name].groupby_bins('dist2sb',np.arange(0,160,interv)).median(skipna=True)
    mean_sec_count[name] = woa_box[name]['SA'].groupby_bins('dist2sb',np.arange(0,160,interv)).count()
    
#############    
# seasonal averages    
mean_sec_seas = {}
mean_sec_seas_count = {}

for name in ['box1','box2','box3','box4']:
    mean_sec_seas[name]={}
    mean_sec_seas_count[name]={}
    for seas in ['DJF','MAM','JJA','SON']:
        mean_sec_seas[name][seas] = woa_box[name].where(woa_box[name].time.dt.season==seas,drop=True).groupby_bins('dist2sb',np.arange(0,160,interv)).mean(skipna=True)
        mean_sec_seas_count[name][seas] = woa_box[name]['SA'].where(woa_box[name].time.dt.season==seas,drop=True).groupby_bins('dist2sb',np.arange(0,160,interv)).count()
#


#%%% functions to plot sections
#    
# Time-mean sections (subplot)
#-----------------------------------------------------------------------------
def plot_box_sections(var,levels,cmap,label,factor=1):
    fig,ax = plt.subplots(figsize=(4,6),nrows=4,sharex=True,sharey=True,constrained_layout=True)
    for name,i in zip(['box1','box2','box3','box4'],
                      range(4)):
        cc=(mean_sec[name][var]*factor).plot.contourf(ax=ax,levels=levels,
                                            cmap=cmap,ylim=(150,0),xlim=(155,0),
                                            y='pres',add_colorbar=False)
        (mean_sec[name]['sal']).plot.contour(ax=ax,levels=[34],
                                            colors='r',linewidths=0.5,ylim=(150,0),xlim=(155,0),
                                            y='pres')
        ccc=(mean_sec[name]['pden']).plot.contour(ax=ax,levels=np.arange(24,27,0.5),
                                            colors='k',linewidths=0.2,ylim=(150,0),xlim=(155,0),
                                            y='pres',linestyles='dashed',alpha=0.5)
        plt.clabel(ccc,fontsize=6,levels=[24,25,26],fmt='%1d')
        
        if var=='Tu':
            (mean_sec[name][var]*factor).plot.contourf(ax=ax,levels=[0,45,90,-45,-90],
                                                cmap=cmap,ylim=(150,0),xlim=(155,0),
                                                y='pres',add_colorbar=False)
        
        # add count to bins
        count = mean_sec_count[name].max(dim='pres')
        for x,j in zip(np.arange(5,155,interv),range(len(np.arange(5,155,interv)))):
            if np.isfinite(count[j]):
                ax.text(x,-2,f'{count[j].values.astype(int)}',ha='center',fontsize=6)
                ax.plot(x,-5,marker='v',markersize=10,color='k',zorder=50)
        # fill in bathymetry
        ax.fill_between(np.arange(5,155,interv),
                           mean_sec[name]['zetopo'],
                           np.ones(len(np.arange(5,155,interv)))*149,
                           color='gray',zorder=50)
        ax.set_xlabel('')
        ax.text(120,130,name,fontweight='bold',zorder=100)
        ax.grid(False)
        ax.set_ylabel('pressure [db]')
        ax.set_xlim(150,0)
        plt.rcParams["axes.axisbelow"] = False
        # ax.minorticks_off()
    
    ax[-1].set_xlabel('distance to shelfbreak [km]')
    cb = plt.colorbar(cc,ax=ax[:],shrink=0.8,label=label,aspect=30,pad=0.03,location='bottom')
    cb.ax.minorticks_off()
    # finished_plot(fig,f'./plots/ECOMON/sections/{interv}km/woa_boxes_subplot_median_sections_{var}.eps')
#-----------------------------------------------------------------------------   
#
# Time-mean sections (individual) 
#-----------------------------------------------------------------------------
def plot_box_sections_individual(var,levels,cmap,label,factor=1):
    plt.rc('axes.spines', top= False, left=False, right=False,bottom=False)
    plt.rc('ytick.major', size= 3, width=0.5 )
    plt.rc('xtick.major', size= 3, width=0.5 )
    plt.close('all')
    for name,i in zip(['box1','box2','box3','box4'],
                      range(4)):
        fig,ax = plt.subplots(figsize=(3.5,2.3),constrained_layout=True)
        cc=(mean_sec[name][var]*factor).plot.contourf(ax=ax,levels=levels,
                                            cmap=cmap,ylim=(150,0),xlim=(155,0),
                                            y='pres',add_colorbar=False)
        (mean_sec[name]['sal']).plot.contour(ax=ax,levels=[34],
                                            colors='r',linewidths=0.5,ylim=(150,0),xlim=(155,0),
                                            y='pres')
        (mean_sec[name]['temp']).plot.contour(ax=ax,levels=[10],
                                            colors='white',linewidths=0.5,ylim=(150,0),xlim=(155,0),
                                            y='pres',linestyle='dotted')
        ccc=(mean_sec[name]['pden']).plot.contour(ax=ax,levels=np.arange(24,27,0.5),
                                            colors='k',linewidths=0.2,ylim=(150,0),xlim=(155,0),
                                            y='pres',linestyles='dashed',alpha=0.5)
        plt.clabel(ccc,fontsize=6,levels=[24,25,26],fmt='%1d')
        
        if var=='Tu':
            (mean_sec[name][var]*factor).plot.contourf(ax=ax,levels=[0,45,90,-45,-90],
                                                cmap=cmap,ylim=(150,0),xlim=(155,0),
                                                y='pres',add_colorbar=False)
        
        # add count to bins
        count = mean_sec_count[name].max(dim='pres')
        for x,j in zip(np.arange(5,155,interv),range(len(np.arange(5,155,interv)))):
            if np.isfinite(count[j]):
                ax.text(x,-2,f'{count[j].values.astype(int)}',ha='center',fontsize=6)
                ax.plot(x,-5,marker='v',markersize=10,color='k',zorder=50)
        # fill in bathymetry
        ax.fill_between(np.arange(5,155,interv),
                           mean_sec[name]['zetopo'],
                           np.ones(len(np.arange(5,155,interv)))*149,
                           color='gray',zorder=50)
        ax.set_xlabel('')
        ax.text(40,140,name,fontweight='bold',zorder=100)
        ax.grid(False)
        ax.set_ylabel('pressure [db]')#,rotation=-90,labelpad=10)
        ax.set_xlim(150,3)
        plt.rcParams["axes.axisbelow"] = False
        # ax.minorticks_off()
        ax.set_yticks(np.arange(0,175,25))
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.minorticks_off()
        ax.set_xlabel('distance to shelfbreak [km]')
        cb = plt.colorbar(cc,ax=ax,shrink=0.8,label=label,aspect=30,pad=0.03,location='bottom')
        cb.ax.minorticks_off()
        # finished_plot(fig,f'./plots/ECOMON/sections/{interv}km/individual4manuscript/woa_{name}_median_section_{var}.eps')
#-----------------------------------------------------------------------------      
#
# Seasonal mean sections (subplot)
#-----------------------------------------------------------------------------   
def plot_box_sections_seas(box,var,levels,cmap,label,factor=1):
    fig,ax = plt.subplots(figsize=(4,6),nrows=4,sharex=True,sharey=True,constrained_layout=True)
    for name,i in zip(['DJF','MAM','JJA','SON'],
                      range(4)):
        cc=(mean_sec_seas[box][name][var]*factor).plot.contourf(ax=ax,levels=levels,
                                            cmap=cmap,ylim=(150,0),xlim=(155,0),
                                            y='pres',add_colorbar=False)
        (mean_sec_seas[box][name]['sal']).plot.contour(ax=ax,levels=[34],
                                            colors='r',linewidths=0.5,ylim=(150,0),xlim=(155,0),
                                            y='pres')
        ccc=(mean_sec_seas[box][name]['pden']).plot.contour(ax=ax,levels=np.arange(24,27,0.5),
                                            colors='k',linewidths=0.2,ylim=(150,0),xlim=(155,0),
                                            y='pres',linestyles='dashed',alpha=0.5)
        plt.clabel(ccc,fontsize=6,levels=[24,25,26],fmt='%1d')
        
        # if var=='Tu':
        #     (mean_sec_seas[box][name][var]*factor).plot.contour(ax=ax,levels=[0,45,90,-45,-90],
        #                                         cmap=cmap,ylim=(150,0),xlim=(155,0),
        #                                         y='pres',add_colorbar=False)
        
        # add count to bins
        count = mean_sec_seas_count[box][name].max(dim='pres')
        for x,j in zip(np.arange(5,155,interv),range(len(np.arange(5,155,interv)))):
            if np.isfinite(count[j]):
                ax.text(x,-2,f'{count[j].values.astype(int)}',ha='center',fontsize=6)
                ax.plot(x,-5,marker='v',markersize=10,color='k',zorder=50)
        # fill in bathymetry
        ax.fill_between(np.arange(5,155,interv),
                           mean_sec_seas[box][name]['zetopo'],
                           np.ones(len(np.arange(5,155,interv)))*149,
                           color='gray',zorder=50)
        ax.set_xlabel('')
        ax.text(120,130,name,fontweight='bold',zorder=100)
        ax.grid(False)
        ax.set_ylabel('pressure [db]')
        ax.set_xlim(150,0)
        plt.rcParams["axes.axisbelow"] = False
        # ax.minorticks_off()
    
    ax[-1].set_xlabel('distance to shelfbreak [km]')
    cb = plt.colorbar(cc,ax=ax[:],shrink=0.8,label=label,aspect=30,pad=0.03,location='bottom')
    cb.ax.minorticks_off()
    # finished_plot(fig,f'./plots/ECOMON/sections/{interv}km/woa_{box}_subplot_median_sections_{var}.jpg')
#-----------------------------------------------------------------------------   
#
# Seasonal mean sections (indidvidual)
#-----------------------------------------------------------------------------
def plot_box_sections_seas_individual(box,var,levels,cmap,label,factor=1):
    # plt.rc('axes.spines', top= False, left=False, right=False,bottom=False)
    plt.rc('ytick.major', size= 3, width=0.5 )
    plt.rc('xtick.major', size= 3, width=0.5 )
    plt.close('all')
    for name,i in zip(['DJF','MAM','JJA','SON'],
                      range(4)):
        fig,ax = plt.subplots(figsize=(3.5,2.3),constrained_layout=True)
        ax.spines[['right', 'top','left','bottom']].set_visible(False)
        cc=(mean_sec_seas[box][name][var]*factor).plot.contourf(ax=ax,levels=levels,
                                            cmap=cmap,ylim=(150,0),xlim=(155,0),
                                            y='pres',add_colorbar=False)
        (mean_sec_seas[box][name]['sal']).plot.contour(ax=ax,levels=[34],
                                            colors='r',linewidths=0.5,ylim=(150,0),xlim=(155,0),
                                            y='pres')
        (mean_sec_seas[box][name]['temp']).plot.contour(ax=ax,levels=[10],
                                            colors='white',linewidths=0.5,ylim=(150,0),xlim=(155,0),
                                            y='pres',linestyle='dotted')
        ccc=(mean_sec_seas[box][name]['pden']).plot.contour(ax=ax,levels=np.arange(24,27,0.5),
                                            colors='k',linewidths=0.2,ylim=(150,0),xlim=(155,0),
                                            y='pres',linestyles='dashed',alpha=0.5)
        plt.clabel(ccc,fontsize=6,levels=[24,25,26],fmt='%1d')
        
        # if var=='Tu':
        #     (mean_sec_seas[box][name][var]*factor).plot.contour(ax=ax,levels=[0,45,90,-45,-90],
        #                                         cmap=cmap,ylim=(150,0),xlim=(155,0),
        #                                         y='pres',add_colorbar=False)
        
        # add count to bins
        count = mean_sec_seas_count[box][name].max(dim='pres')
        for x,j in zip(np.arange(5,155,interv),range(len(np.arange(5,155,interv)))):
            if np.isfinite(count[j]):
                ax.text(x,-2,f'{count[j].values.astype(int)}',ha='center',fontsize=6)
                ax.plot(x,-5,marker='v',markersize=10,color='k',zorder=50)
        # fill in bathymetry
        ax.fill_between(np.arange(5,155,interv),
                           mean_sec_seas[box][name]['zetopo'],
                           np.ones(len(np.arange(5,155,interv)))*149,
                           color='gray',zorder=50)
        ax.set_xlabel('')
        ax.text(40,140,name,fontweight='bold',zorder=100)
        ax.grid(False)
        ax.set_ylabel('pressure [db]')#,rotation=-90,labelpad=10)
        ax.set_xlim(150,3)
        plt.rcParams["axes.axisbelow"] = False
        # ax.minorticks_off()
        ax.set_yticks(np.arange(0,175,25))
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.minorticks_off()
        ax.set_xlabel('distance to shelfbreak [km]')
        cb = plt.colorbar(cc,ax=ax,shrink=0.8,label=label,aspect=30,pad=0.03,location='bottom')
        cb.ax.minorticks_off()
        # finished_plot(fig,f'./plots/ECOMON/sections/{interv}km/individual4manuscript/woa_{box}_median_sections_{var}_{name}.eps')

#%%% plot time-mean median sections
plt.close('all')
font_for_print()
plot_box_sections_individual('temp',np.arange(8,17,0.25),'Thermal','pot. temperature [\N{DEGREE SIGN}C]')
plot_box_sections_individual('sal',np.arange(32,36,0.2),'Spectral_r','salinity [psu]')
# plot_box_sections_individual('pden',np.arange(24,26.8,0.2),'Haline','pot. density [kg/m$^{-3}$]')
# plot_box_sections_individual('CT_contr',np.arange(-3,3,0.1),'RdBu_r','CT_contr*$10^{-4}$',factor=1e04)
# plot_box_sections_individual('SA_contr',np.arange(-3,3,0.1),'RdBu_r','SA_contr*$10^{-4}$',factor=1e04)
# plot_box_sections_individual('bfrq',np.arange(0,3,0.1),'Spectral_r','N2*$10^{-4}$',factor=1e04)
# plot_box_sections_individual('spice_strat',np.arange(-3,3,0.2),'rdbu_r','spice stratification [s$^{-2}$]*$10^{-4}$',factor=1e04)
# plot_box_sections_individual('Tu',np.arange(-60,65,5),'RdBu_r','Tu')
# plot_box_sections_individual('SCI',np.arange(-3,3,0.2),'Rdbu_r','SCI')


#%%% plot seasonal-mean median sections
plt.rc('axes.spines', top= True, right=True )

for name in ['box1','box2','box3','box4'][:1]:
    plt.close('all')    
    plot_box_sections_seas_individual(name,'temp',np.arange(5,19,0.25),'Thermal','pot. temperature [\N{DEGREE SIGN}C]')
    plot_box_sections_seas_individual(name,'sal',np.arange(32,36,0.2),'Spectral_r','salinity [psu]')
    # plot_box_sections_seas_individual(name,'pden',np.arange(24,26.8,0.2),'Haline','pot. density [kg/m^{-3}]')
    # plot_box_sections_seas_individual(name,'CT_contr',np.arange(-5,5.1,0.1),'Rdbu_r','CT_contr*$10^{-4}$',factor=1e04)
    # plot_box_sections_seas_individual(name,'SA_contr',np.arange(-5,5.1,0.1),'Rdbu_r','SA_contr*$10^{-4}$',factor=1e04)
    # plot_box_sections_seas_individual(name,'bfrq',np.arange(0,10.1,0.1),'spectral_r','N2 [s$^{-2}$]*$10^{-4}$',factor=1e04)
    # plot_box_sections_seas_individual(name,'SCI',np.arange(-4,4.2,0.2),'Rdbu_r','SCI')
    # # plot_box_sections_seas_individual(name,'Tu',30,'Rdbu_r','Tu')
    # plot_box_sections_seas_individual(name,'spice_strat',np.arange(-3,3.2,0.2),'rdbu_r','spice stratification [s$^{-2}$]*$10^{-4}$',factor=1e04)

## =============================================================================
#
#
#
# =============================================================================
#%% FUNCTION for seasonal histogram plots
# =============================================================================

def plot_hist(ds,var,locvar,locrange,botdepthrange,presrange,bins=np.arange(29,36,0.2)):
    """
    
    Parameters
    ----------
    ds :                xarray data array (i.e. variable already chosen)
    sal :               string, variable name  
    lonrange :          min and max of longitude to cut data
    botdepthrange :     range of bottom depth to limit to e.g. inner or outer shelf
    presrange :         min max of pressure range to focus on specific levels in water column

    Returns
    -------
    None.

    """
    
    # inititalize hist
    fig,ax = plt.subplots(figsize=(10,8),ncols=4,nrows=3,sharex=True,sharey=True,constrained_layout=True)
    axh = ax.ravel()
    
    # inititalize map
    ax2 = fig.add_axes([0.23, 0.8, 0.1, 0.1])
    ax2.plot(ds['lon'],ds['lat'],marker='.',linestyle='',color='gray')
    
    # loop through months
    for i in np.arange(1,13):
        # woa['sal'].where((woa.lon>=-72) & (woa.lon<=-70) & (woa.botdepth<50),drop=True).where(woa.time.dt.month==i,drop=True).plot.line(ax=axh[i-1],hue='profile',add_legend=False,y='pres')
        dummy = ds.where((ds[locvar]>=locrange[0]) & 
                        (ds[locvar]<=locrange[1]) &
                        (ds.botdepth>botdepthrange[0]) &
                        (ds.botdepth<botdepthrange[1]) &
                        (ds.time.dt.month==i),
                        drop=True).sel(pres=slice(*presrange))[var]
        if np.shape(dummy)[1]>0:
            dummy.plot.hist(ax=axh[i-1],bins=bins)
            axh[i-1].text(0.05,0.9,f'# profiles: {np.shape(dummy)[1]}',fontsize=8,transform=axh[i-1].transAxes)
            axh[i-1].text(0.8,0.9,month_converter(i),transform=axh[i-1].transAxes,fontweight='bold',fontsize=8)
            axh[i-1].axvline(dummy.median().values,color='r',ymax=100)
            axh[i-1].grid(False)
            
            # prettify
            if i<9: axh[i-1].set_xlabel('')
            else: axh[i-1].set_xlabel(var)
        
            # plot on map
            ax2.plot(dummy['lon'],dummy['lat'],marker='.',linestyle='',color=plt.get_cmap('tab10')(0))
        ax2.set_xticklabels('')
        ax2.set_yticklabels('')
        ax2.text(0.3,0.2,f'{botdepthrange[0]}-{botdepthrange[1]}m depth',
                 transform=ax2.transAxes,fontsize=6)

    plt.suptitle(f"Pressure range: {presrange[0]}-{presrange[1]}")  
    return fig
        

#%% FUNCTIONS: Box and whisker plots
# =============================================================================
def plot_boxplot(ax,ds,label,var='sal'):
    import matplotlib.patches as mpatches
    cols = plt.get_cmap('tab10')
    font_for_pres()
    
    for i in np.arange(1,13):
        dummy = ds.where(ds.time.dt.month==i).mean(dim='pres',skipna=True).stack(z=['profile'])

        if i==1:
            pdbox = pd.DataFrame({f'{i}': dummy})
            # pdbox2 = pd.DataFrame({f'{i}': dummy})
        else:
            pdbox[f'{i}'] = dummy
        
        
    if var=='sal':
        ylim=[29,35.8]
        ylabel='salinity'
    elif var=='temp':
        ylim=[3,28] 
        ylabel='temperature [\N{DEGREE SIGN}C]'
        
    # plotting
    bp_dict = pdbox.plot(kind='box',ylim=ylim,ax=ax,xticks=np.arange(1,13),
                         patch_artist=True,return_type='both',sym='.')
    months = [month_converter(i+1) for i in range(12)]
    count1 = pdbox.count().values
    ticks = np.arange(1,13)

    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{months[i]}' for i in range(12)])#\n ({count1[i]}) \n ({count2[i]})
    for i in np.arange(0,12):
        ax.text(ticks[i], 28.2, f"({count1[i]})", size=7, ha='center',color=cols(0))
    fig.subplots_adjust(bottom=0.15)     # Add space at bottom
    
    # for row_key, (ax,row) in bp_dict.iteritems():
    for i in np.arange(0,len(bp_dict[1]['boxes'])):
        bp_dict[1]['boxes'][i].set_facecolor('lightblue')
        # bp_dict[1]['fliers'][i].set_color('blue')
        bp_dict[1]['medians'][i].set_color(cols(0))
    # for i in np.arange(4,len(bp_dict[1]['whiskers']),6):
    #     bp_dict[1]['fliers'][i].set_color('lightgray')
    # bp_dict[1]['fliers'][0].set_color('gray')
    ax.set_ylabel(ylabel)
    
    # #add legend
    # patch1 = mpatches.Patch(color=cols(0),label=label[0])
    # # patch2 = mpatches.Patch(color=cols(1),label=label[1])
    # plt.legend(handles=patch1,loc='lower left',fontsize=6)

    #add map
    ax2 = ax.inset_axes([0.8, 0.1, 0.15, 0.15])
    ax2.plot(woa['lon'],woa['lat'],marker='.',linestyle='',color='gray')
    ax2.plot(ds.dropna(dim='profile')['lon'],ds.dropna(dim='profile')['lat'],
              marker='.',linestyle='',color='lightblue',markersize=1)
    ax2.set_xticks([-75,-70])
    ax2.set_xticklabels(['75W','70W'],fontsize=7)
    ax2.set_yticklabels(['35N','40N'],fontsize=7)
    ax2.minorticks_off()
    ax2.set_xlim(-76.5,-65.8)
    ax2.set_ylim(35,42.1)
    return fig,ax,pdbox




def plot_boxplot2(ax,ds,ds2,label,var='sal'):
    import matplotlib.patches as mpatches
    cols = plt.get_cmap('tab10')
    font_for_pres()
    
    for i,j in zip(np.arange(1,13),np.arange(1,26,2)):
        dummy = ds.where(ds.time.dt.month==i).mean(dim='pres',skipna=True).stack(z=['profile'])
        dummy2 = ds2.where(ds2.time.dt.month==i).mean(dim='pres',skipna=True).stack(z=['profile'])

        if i==1:
            pdbox = pd.DataFrame({f'{j}': dummy})
            # pdbox2 = pd.DataFrame({f'{i}': dummy})
        else:
            pdbox[f'{j}'] = dummy
            # pdbox2[f'{i}'] = dummy2
        pdbox[f'{j+0.5}'] = dummy2 
        pdbox[f'{j+1}'] = np.ones(np.size(dummy2)) *np.nan
        
        
    if var=='sal':
        ylim=[30,36.3]
        ylabel='salinity'
    elif var=='temp':
        ylim=[3,28] 
        ylabel='temperature [\N{DEGREE SIGN}C]'
        
    # plotting
    # fig,ax = plt.subplots(figsize=(8,4))
    bp_dict = pdbox.plot(kind='box',ylim=ylim,ax=ax,xticks=np.arange(1.5,35,3),
                         patch_artist=True,return_type='both',sym='.')
    months = [month_converter(i+1) for i in range(12)]
    count1 = pdbox.iloc[:,0:-1:3].count().values
    count2 = pdbox.iloc[:,1:-1:3].count().values
    ticks = np.arange(1.5,35,3)

    ax.set_xticklabels([f'{months[i]}' for i in range(12)] )#\n ({count1[i]}) \n ({count2[i]})
    for i, x in enumerate(ax.get_xticks()):
        ax.text(ticks[i], 29.1, f"({count1[i]})", size=7, ha='center',color=cols(0))
        ax.text(ticks[i], 28.8, f"({count2[i]})", size=7, ha='center',color=cols(1))
    fig.subplots_adjust(bottom=0.15)     # Add space at bottom
    # pdbox2.plot(kind='box',ylim=[29,36.2],ax=ax[1])
    # ax.grid()
    # for row_key, (ax,row) in bp_dict.iteritems():
    for i in np.arange(1,len(bp_dict[1]['boxes']),3):
        bp_dict[1]['boxes'][i].set_facecolor(cols(1))
        # bp_dict[1]['fliers'][i].set_color('blue')
        bp_dict[1]['medians'][i].set_color(cols(0))
    # for i in np.arange(4,len(bp_dict[1]['whiskers']),6):
    #     bp_dict[1]['fliers'][i].set_color('lightgray')
    bp_dict[1]['fliers'][0].set_color('gray')
    ax.set_ylabel(ylabel)
    ax.minorticks_off()
    
    #add legend
    patch1 = mpatches.Patch(color=cols(0),label=label[0])
    patch2 = mpatches.Patch(color=cols(1),label=label[1])
    plt.legend(handles=[patch1,patch2],loc='lower right',fontsize=8)
    
    # #add map
    # ax2 = ax.inset_axes([0.82, 0.15, 0.15, 0.15])
    # ax2.plot(woa['lon'],woa['lat'],marker='.',linestyle='',color='gray')
    # ax2.plot(ds.dropna(dim='profile')['lon'],ds.dropna(dim='profile')['lat'],
    #          marker='.',linestyle='',color=cols(3),markersize=1)
    # ax2.plot(ds2.dropna(dim='profile')['lon'],ds2.dropna(dim='profile')['lat'],
    #          marker='.',linestyle='',color=cols(3),markersize=1)
    # ax2.set_xticks([-75,-70])
    # ax2.set_xticklabels(['75W','70W'],fontsize=7)
    # ax2.set_yticklabels(['35N','40N'],fontsize=7)
    # ax2.set_xlim(-76.5,-65.8)
    # ax2.set_ylim(35,42.1)
    # return fig,ax,pdbox


def plot_boxplot3(ax,ds,label,var='sal'):
    """
    Does not average across specified depth range like functions above

    """
    import matplotlib.patches as mpatches
    cols = plt.get_cmap('tab10')
    font_for_pres()
    
    for i in np.arange(1,13):
        dummy = ds.where(ds.time.dt.month==i).stack(z=['profile','pres'])

        if i==1:
            pdbox = pd.DataFrame({f'{i}': dummy})
            # pdbox2 = pd.DataFrame({f'{i}': dummy})
        else:
            pdbox[f'{i}'] = dummy
        
        
    if var=='sal':
        ylim=[29,35.8]
        ylabel='salinity'
    elif var=='temp':
        ylim=[3,28] 
        ylabel='temperature [\N{DEGREE SIGN}C]'
        
    # plotting
    bp_dict = pdbox.plot(kind='box',ylim=ylim,ax=ax,xticks=np.arange(1,13),
                         patch_artist=True,return_type='both',sym='.')
    months = [month_converter(i+1) for i in range(12)]
    count1 = pdbox.count().values
    ticks = np.arange(1,13)

    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{months[i]}' for i in range(12)])#\n ({count1[i]}) \n ({count2[i]})
    for i in np.arange(0,12):
        ax.text(ticks[i], 28.2, f"({count1[i]})", size=7, ha='center',color=cols(0))
    fig.subplots_adjust(bottom=0.15)     # Add space at bottom
    
    # for row_key, (ax,row) in bp_dict.iteritems():
    for i in np.arange(0,len(bp_dict[1]['boxes'])):
        bp_dict[1]['boxes'][i].set_facecolor(cols(1))
        # bp_dict[1]['fliers'][i].set_color('blue')
        bp_dict[1]['medians'][i].set_color(cols(0))
    # for i in np.arange(4,len(bp_dict[1]['whiskers']),6):
    #     bp_dict[1]['fliers'][i].set_color('lightgray')
    # bp_dict[1]['fliers'][0].set_color('gray')
    ax.set_ylabel(ylabel)
    
    # #add legend
    # patch1 = mpatches.Patch(color=cols(0),label=label[0])
    # # patch2 = mpatches.Patch(color=cols(1),label=label[1])
    # plt.legend(handles=patch1,loc='lower left',fontsize=6)

    #add map
    ax2 = ax.inset_axes([0.8, 0.05, 0.15, 0.15])
    ax2.plot(woa['lon'],woa['lat'],marker='.',linestyle='',color='gray')
    ax2.plot(ds.dropna(dim='profile')['lon'],ds.dropna(dim='profile')['lat'],
              marker='.',linestyle='',color=cols(3),markersize=1)
    ax2.set_xticks([-75,-70])
    ax2.set_xticklabels(['75W','70W'],fontsize=7)
    ax2.set_yticklabels(['35N','40N'],fontsize=7)
    ax2.minorticks_off()
    ax2.set_xlim(-76.5,-65.8)
    ax2.set_ylim(35,42.1)
    return fig,ax,pdbox
    
    
def subset_ds(ds,botdepth,presrange,loc,varloc,timerange=[1990,2023],var='sal'):
    ds = ds[var].where((woa[varloc]>=loc[0]) &
                          (woa[varloc]<=loc[1]) & 
                           (woa.time.dt.year>=timerange[0]) &
                           (woa.time.dt.year<timerange[1]) & 
                          (woa.botdepth>botdepth[0]) &
                          (woa.botdepth<botdepth[1])).sel(pres=slice(*presrange))
    return ds

# =============================================================================   
#%% Box plots
# =============================================================================
plt.close('all')

presrange = [0,50]
loc = [36,42] # first number has to be smaller than second
varloc = 'lat'
# fig,ax = plt.subplots(figsize=(8,8),sharex=True,nrows=4,sharey=True)
for i,name in zip(range(1),['box3']):#zip(range(4),['box1','box2','box3','box4']):
    # fig,ax = plt.subplots(figsize=(8,4))
    # plot_boxplot(ax,
    #               subset_ds(woa_box[name],botdepth=[0,200],
    #                             presrange=presrange,loc=loc,varloc=varloc,
    #                             timerange=[1990,2020]),
    #               label=['1990-2023'])
    # ax.set_title(f"{name}:{presrange[0]}-{presrange[1]}m")
    # finished_plot(fig,f'./plots/ECOMON/box_plot_{name}_{presrange[0]}-{presrange[1]}m.jpg')
    
    fig,ax = plt.subplots(figsize=(5,3))
    plot_boxplot2(ax,
                 subset_ds(ds=woa_box[name],botdepth=[0,200],
                                presrange=presrange,loc=loc,varloc=varloc,
                                timerange=[1990,2000]),
                 subset_ds(ds=woa_box[name],botdepth=[0,200],
                                presrange=presrange,loc=loc,varloc=varloc,
                                timerange=[2010,2020]),
                 label=['1990-2000','2010-2020'])
    ax.set_title(f"{name}:{presrange[0]}-{presrange[1]}m")
    ax.grid(False)
    finished_plot(fig,f'./plots/manuscript/box_plot_{name}_{presrange[0]}-{presrange[1]}m_P1vsP2.jpg')
    # 
    
# =============================================================================
#%% TS diagram for 5 annual averages over boxes
# =============================================================================
plt.close('all')
fig,ax = plt.subplots(figsize=(5,5))
ds = woa
cols = plt.get_cmap('tab10')
i=0
for year in np.arange(1990,2020,5):
    ds_mean = ds.where((ds.time.dt.year>=year) & (ds.time.dt.year<year+5),drop=True).mean('profile',skipna=True)
    plt.plot(ds_mean['sal'],ds_mean['temp'],label=str(year),color=cols(i))
    ax.set_xlim(32,36)
    ax.set_ylim(7,15)
    i=i+1
ts_append(np.arange(22,28,0.5),ax)
plt.legend()
# finished_plot(fig,'ts_diagram_woa_5y_interval.jpg',printpath=printpath)


#=============================================================================
#%% derive annual mean values
#
# Want to apply a monthly weighted average 
# =============================================================================

var='sal'
def annual_mean_timeseries(ds):
    var_am = np.ones(200)*np.nan
    count_am = np.ones(200)*np.nan
    for year in range(1990,2020):
        dummy = ds.where(ds.time.dt.year==year,drop=True)
        if np.shape(dummy.profile)[0]>0:
            count = dummy.groupby('time.month').count()
            wght = count/count.sum(dim='month')
            var_am = np.column_stack([var_am,
                                      (dummy.groupby('time.month').mean('profile')*wght).sum('month').values])
            count_am = np.column_stack([count_am,count.reindex({"month": np.arange(1,13)}, fill_value=0)])
        else:
            var_am = np.column_stack([var_am,np.ones(200)*np.nan])
            count_am = np.column_stack([count_am,np.zeros([200,12])])

    return var_am[:,1::],count_am[:,1::]

# derive timeseries for different seasons
def seasonal_mean_timeseries(ds):
    var_am = np.ones(200)*np.nan
    count_am = np.ones(200)*np.nan
    for year in range(1990,2020):
        dummy = ds.where(ds.time.dt.year==year,drop=True)
        if np.shape(dummy.profile)[0]>0:
            count = dummy.groupby('time.month').count()
            wght = count/count.sum(dim='month')
            var_am = np.column_stack([var_am,
                                      (dummy.groupby('time.month').mean('profile')*wght).sum('month').values])
            count_am = np.column_stack([count_am,count.reindex({"month": np.arange(1,13)}, fill_value=0)])
        else:
            var_am = np.column_stack([var_am,np.ones(200)*np.nan])
            count_am = np.column_stack([count_am,np.zeros([200,12])])

    return var_am[:,1::],count_am[:,1::]


woa_box_am = {}
woa_box_am_count = {}
for name in list(woa_box.keys()):
    print(name)
    sal_am,sal_count = annual_mean_timeseries(woa_box[name]['sal'])
    temp_am,temp_count = annual_mean_timeseries(woa_box[name]['temp'])
    #
    # create data set
    woa_box_am[name] = xr.Dataset(data_vars=dict(sal=(["pres","time"],sal_am),
                                                 temp=(["pres","time"],temp_am)),
                            coords=dict(pres=(["pres"],np.arange(1,201)),
                                        time=(["time"],pd.date_range(start='1990',end='2019',freq='YS'))
                                    )
                        )
    woa_box_am_count[name] = xr.Dataset(data_vars=dict(sal_count=(["pres","time"],sal_count),
                                                 temp_count=(["pres","time"],temp_count)),
                            coords=dict(pres=(["pres"],np.arange(1,201)),
                                        time=(["time"],pd.date_range(start='1990',end='2019-12',freq='MS'))
                                    )
                        )
    
#%%
def datetime_tick(ax,axis='xaxis',interval='year'):
    eval(f"ax.{axis}.set_major_locator(mdates.YearLocator(5))")
    eval(f"ax.{axis}.set_major_formatter(mdates.DateFormatter('%Y'))")  

plt.close('all')
cols = plt.get_cmap('spectral',10)
norm = colors.Normalize(vmin=5,vmax=50)
sm = plt.cm.ScalarMappable(cmap=cols, norm=norm)
sm.set_array([])
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.Spectral(np.linspace(0,1,10)))


#############3
# Salinity
fig,ax = plt.subplots(figsize=(5,8),nrows=4,sharex=True,sharey=False)
for name,i in zip(list(woa_box.keys()),range(4)):
    dummy = woa_box_am[name]['sal'].groupby_bins('pres',np.arange(0,55,5)).mean(skipna=True)
    (dummy).plot(ax=ax,hue='pres_bins',
                                        add_legend=False,
                                        linewidth=1)
    
    datetime_tick(ax)
    # ax.set_xlim(datetime.date(1990,1,1),datetime.date(2020,1,1))
    ax.set_xlabel('')
    ax.set_title(name)
plt.xticks(rotation=0,ha='center')
# add colorbar
cbax = ax[-1].inset_axes([0.1, 0.8, 0.3, 0.05])
cb = plt.colorbar(sm,cax=cbax,ticks=np.arange(0,55,5),orientation='horizontal', 
                 boundaries=np.arange(2.5,57.5,5))
cb.ax.tick_params(labelsize=6)
cb.ax.minorticks_off()
cbax.xaxis.tick_top()
print_path(printpath,fig)
finished_plot(fig,'./plots/ECOMON/annual_mean_salinity_timeseries_boxes.jpg')

# per season
fig,ax = plt.subplots(figsize=(5,8),nrows=4,sharex=True,sharey=False)
for name,i in zip(list(woa_box.keys()),range(4)):
    dummy = woa_box_am[name]['sal'].groupby_bins('pres',np.arange(0,55,5)).mean(skipna=True)
    (dummy).plot(ax=ax,hue='pres_bins',
                                        add_legend=False,
                                        linewidth=1)
    
    datetime_tick(ax)
    # ax.set_xlim(datetime.date(1990,1,1),datetime.date(2020,1,1))
    ax.set_xlabel('')
    ax.set_title(name)
plt.xticks(rotation=0,ha='center')
# add colorbar
cbax = ax[-1].inset_axes([0.1, 0.8, 0.3, 0.05])
cb = plt.colorbar(sm,cax=cbax,ticks=np.arange(0,55,5),orientation='horizontal', 
                 boundaries=np.arange(2.5,57.5,5))
cb.ax.tick_params(labelsize=6)
cb.ax.minorticks_off()
cbax.xaxis.tick_top()
print_path(printpath,fig)
finished_plot(fig,'./plots/ECOMON/annual_mean_salinity_timeseries_boxes.jpg')


#%%############3
# Temperature
fig,ax = plt.subplots(figsize=(5,8),nrows=4,sharex=True,sharey=False)
for name,i in zip(list(woa_box.keys()),range(4)):
    dummy = woa_box_am[name]['temp'].groupby_bins('pres',np.arange(0,55,5)).mean()
    (dummy).plot(ax=ax,hue='pres_bins',
                                        add_legend=False,
                                        linewidth=1)
    
    datetime_tick(ax)
    # ax.set_xlim(datetime.date(1990,1,1),datetime.date(2020,1,1))
    ax.set_xlabel('')
    ax.set_title(name)
plt.xticks(rotation=0,ha='center')
# add colorbar
cbax = ax[-1].inset_axes([0.1, 0.8, 0.3, 0.05])
cb = plt.colorbar(sm,cax=cbax,ticks=np.arange(0,55,5),orientation='horizontal', 
                 boundaries=np.arange(2.5,57.5,5))
cb.ax.tick_params(labelsize=6)
cb.ax.minorticks_off()
cbax.xaxis.tick_top()
print_path(printpath,fig)
# finished_plot(fig,'./plots/ECOMON/annual_mean_temperature_timeseries_boxes.jpg')


#%% ########################
# Count per year
fig,ax = plt.subplots(figsize=(5,8),nrows=4,sharex=True,sharey=True)
for name,i in zip(list(woa_box.keys()),range(4)):
    woa_box_am_count[name]['sal_count'][0,:].groupby('time.year').sum().to_series().plot.bar(ax=ax)
    ax.set_title(name)
print_path(printpath,fig)
# finished_plot(fig,'./plots/ECOMON/annual_mean_count_timeseries_boxes.jpg')
























#%%  
# Testing
# =============================================================================

loc = [-72,-70]
presrange = [0,5]
varloc='lon'
fig,ax,data = plot_boxplot(subset_ds(botdepth=[0,200],presrange=presrange,loc=loc,varloc=varloc,
                                     timerange=[1990,2000],
                                     var='temp'),
                           subset_ds([0,200],presrange=presrange,loc=loc,varloc=varloc,
                                     timerange=[2010,2020],
                                     var='temp'),
                           label=['1990-2000','2010-2019'],
                           var='temp')

#%% mean profiles

test = subset_ds(botdepth=[50,100],presrange=[0,100],loc=[-72,-70],varloc='lon',timerange=[1990,2000],var='temp')
test2 = subset_ds(botdepth=[50,100],presrange=[0,100],loc=[-72,-70],varloc='lon',timerange=[2010,2023],var='temp')

plt.figure()
plt.plot(np.nanmean(test,axis=1),np.arange(0,100))
plt.plot(np.nanmean(test2,axis=1),np.arange(0,100))
plt.gca().invert_yaxis()

# check how many values for each level
dummy =test.values
dummy[np.isfinite(dummy)]=1
np.nansum(dummy,axis=1)


# =============================================================================
#%% Mean seasonal cycle Salinity in different time periods - depth levels 5m interval
# =============================================================================

varloc = 'lon'  # 'lat'
loc =  [-72,-70] #[37,39] # 
# initialize empty list to count number of values that go into average
inshore_count = []
offshore_count = []

font_for_print()
fig,ax = plt.subplots(nrows=2,sharex=True,figsize=(10*cm,8*cm),sharey=True)
cols = plt.get_cmap('blues',13)
norm = colors.Normalize(vmin=0,vmax=50)
sm = plt.cm.ScalarMappable(cmap=cols, norm=norm)
sm.set_array([])
j=0

# add maps
cols2 = plt.get_cmap('tab10')
ax2 = fig.add_axes([0.85, 0.4, 0.1, 0.1])
ax2.plot(woa['lon'],woa['lat'],marker='.',linestyle='',color='gray')

ax3 = fig.add_axes([0.85, 0.85, 0.1, 0.1])
ax3.plot(woa['lon'],woa['lat'],marker='.',linestyle='',color='gray')


for i in np.arange(0,55,5):
    sal_inshore = subset_ds(botdepth=[0,200],presrange=[i,i+4],loc=loc,varloc=varloc,timerange=[1990,2000])
    sal_inshore.mean('pres',skipna=True).groupby('time.month').mean(skipna=True).plot(ax=ax[0],label=str(i),
                                                                                      color=cols(j+1),linewidth=0.5,
                                                                                      linestyle='-')
    inshore_count.append(sal_inshore.mean('pres',skipna=True).groupby('time.month').count().interp(month=np.arange(1,13),
                                                                                         kwargs={"fill_value": 0}).values)
    ax3.plot(sal_inshore.mean('pres',skipna=True).dropna(dim='profile')['lon'],
             sal_inshore.mean('pres',skipna=True).dropna(dim='profile')['lat'],
             marker='.',linestyle='',color=cols2(0),markersize=1)
    sal_offshore = subset_ds(botdepth=[0,200],presrange=[i,i+4],loc=loc,varloc=varloc,timerange=[2010,2019])
    sal_offshore.mean('pres',skipna=True).groupby('time.month').mean(skipna=True).plot(ax=ax[1],label=str(i),
                                                                                       color=cols(j+1),linewidth=0.5,
                                                                                       linestyle='-')
    offshore_count.append(sal_offshore.mean('pres',skipna=True).groupby('time.month').count().interp(month=np.arange(1,13),
                                                                                         kwargs={"fill_value": 0}).values)
    j=j+1

ax[0].set_xticks(np.arange(1,13))
labels=[month_converter(i+1) for i in range(12)]
ax[0].set_xticklabels(labels)
# ax[0].set_xlabel('')
# ax[1].set_xlabel('')
# ax[0].set_title('Inshore (0-50m)')
# ax[1].set_title('Offshore (50-150m)')
ax[0].set_title('1990-2005')
ax[1].set_title('2010-2019')
# ax[0].grid(False)
# ax[1].grid(False)

for axh in ax:
    axh.grid(False)
    axh.set_xlabel('')
    axh.set_ylim([31,35])

# add colorbar
cbax = fig.add_axes([0.15, 0.17, 0.3, 0.03])
cb = plt.colorbar(sm,cax=cbax,ticks=np.arange(0,55,5),orientation='horizontal', 
                 boundaries=np.arange(2.5,57.5,5))
cb.ax.tick_params(labelsize=6)
cbax.xaxis.tick_top()

# add data to map
# ax2.plot(sal_inshore.dropna(dim='profile')['lon'],
#          sal_inshore.dropna(dim='profile')['lat'],
#          marker='.',linestyle='',color=cols(0),markersize=1)
ax2.plot(sal_offshore.dropna(dim='profile')['lon'],
          sal_offshore.dropna(dim='profile')['lat'],
          marker='.',linestyle='',color=cols2(1),markersize=1)

for axh in [ax2,ax3]:
    axh.set_xticklabels('')
    axh.set_yticklabels('')
    axh.set_xlim(-76.5,-65.8)
    axh.set_ylim(35,42.1)
    axh.grid(False)
    
    
count1 = np.round(np.mean(inshore_count[0:6],axis=0),0).astype(int)
count2 = np.round(np.mean(offshore_count,axis=0),0).astype(int)

# plt.text(0.2,0.01,'# profiles: ' + str(np.round(np.mean(offshore_count,axis=0),0).astype(int)),
#          transform = fig.transFigure,fontsize=5)
# plt.text(0.2,0.52,'# profiles: ' + str(np.round(np.mean(inshore_count[0:6],axis=0),0).astype(int)),
#          transform = fig.transFigure,fontsize=5)
# ax[1].set_xticklabels([f'{months[i]}' for i in range(12)] )#\n ({count1[i]}) \n ({count2[i]})
cols=plt.get_cmap('tab10')
ticks = np.arange(1,13)
for i, x in enumerate(ax[1].get_xticks()):
    ax[1].text(ticks[i], 30.1, f"({count1[i]})", size=5, ha='center',color=cols(0))
    ax[1].text(ticks[i], 29.8, f"({count2[i]})", size=5, ha='center',color=cols(1))
fig.subplots_adjust(bottom=0.15)     # Add space at bottom
# finished_plot(fig,'./plots/woa_sal_mean_seasonal_cycle_37-39N_depth_levels_p1_vs_p2.jpg')# 
# finished_plot(fig,'./plots/woa_sal_mean_seasonal_cycle_70-72W_depth_levels_p1_vs_p2.jpg')


 
# =============================================================================
#%% Mean seasonal cycle Temperature in different time periods - depth levels 5m interval
# =============================================================================

varloc = 'lon'  # 'lat'
loc =  [-72,-70] #[37,39] # 
# initialize empty list to count number of values that go into average
inshore_count = []
offshore_count = []

font_for_print()
fig,ax = plt.subplots(nrows=2,sharex=True,figsize=(10*cm,8*cm),sharey=True)
cols = plt.get_cmap('Blues',13)
norm = colors.Normalize(vmin=0,vmax=50)
sm = plt.cm.ScalarMappable(cmap=cols, norm=norm)
sm.set_array([])
j=0

# add maps
cols2 = plt.get_cmap('tab10')
ax2 = fig.add_axes([0.85, 0.4, 0.1, 0.1])
ax2.plot(woa['lon'],woa['lat'],marker='.',linestyle='',color='gray')

ax3 = fig.add_axes([0.85, 0.85, 0.1, 0.1])
ax3.plot(woa['lon'],woa['lat'],marker='.',linestyle='',color='gray')


for i in np.arange(0,55,5):
    sal_inshore = subset_ds(botdepth=[0,200],presrange=[i,i+4],loc=loc,varloc=varloc,timerange=[1990,2005],var='temp')
    sal_inshore.mean('pres',skipna=True).groupby('time.month').mean(skipna=True).plot(ax=ax[0],label=str(i),
                                                                                      color=cols(j+1),linewidth=0.5,
                                                                                      linestyle='-')
    inshore_count.append(sal_inshore.mean('pres',skipna=True).groupby('time.month').count().interp(month=np.arange(1,13),
                                                                                     kwargs={"fill_value": 0}).values)
    ax3.plot(sal_inshore.mean('pres',skipna=True).dropna(dim='profile')['lon'],
             sal_inshore.mean('pres',skipna=True).dropna(dim='profile')['lat'],
             marker='.',linestyle='',color=cols2(0),markersize=1)
    sal_offshore = subset_ds(botdepth=[0,200],presrange=[i,i+4],loc=loc,varloc=varloc,timerange=[2010,2019],var='temp')
    sal_offshore.mean('pres',skipna=True).groupby('time.month').mean(skipna=True).plot(ax=ax[1],label=str(i),
                                                                                       color=cols(j+1),linewidth=0.5,
                                                                                       linestyle='-')
    offshore_count.append(sal_offshore.mean('pres',skipna=True).groupby('time.month').count().interp(month=np.arange(1,13),
                                                                                     kwargs={"fill_value": 0}).values)
    j=j+1
    
ax[0].set_xticks(np.arange(1,13))
labels=[month_converter(i+1) for i in range(12)]
ax[0].set_xticklabels(labels)
# ax[0].set_xlabel('')
# ax[1].set_xlabel('')
# ax[0].set_title('Inshore (0-50m)')
# ax[1].set_title('Offshore (50-150m)')
ax[0].set_title('1990-2005')
ax[1].set_title('2010-2019')
# ax[0].grid(False)
# ax[1].grid(False)

for axh in ax:
    axh.grid(False)
    axh.set_xlabel('')
    # axh.set_ylim([31,35])
    axh.axhline(10,color='k',alpha=0.3,zorder=0,linewidth=0.5)

# add colorbar
cbax = fig.add_axes([0.15, 0.35, 0.3, 0.03])
cb = plt.colorbar(sm,cax=cbax,ticks=np.arange(0,55,5),orientation='horizontal', 
                 boundaries=np.arange(2.5,57.5,5))
cb.ax.tick_params(labelsize=6)
cbax.xaxis.tick_top()

# add data to map
# ax2.plot(sal_inshore.dropna(dim='profile')['lon'],
#          sal_inshore.dropna(dim='profile')['lat'],
#          marker='.',linestyle='',color=cols(0),markersize=1)
ax2.plot(sal_offshore.dropna(dim='profile')['lon'],
          sal_offshore.dropna(dim='profile')['lat'],
          marker='.',linestyle='',color=cols2(1),markersize=1)

for axh in [ax2,ax3]:
    axh.set_xticklabels('')
    axh.set_yticklabels('')
    axh.set_xlim(-76.5,-65.8)
    axh.set_ylim(35,42.1)
    axh.grid(False)
    
    
count1 = np.round(np.mean(inshore_count[0:6],axis=0),0).astype(int) # only count upper 6 depth levels
count2 = np.round(np.mean(offshore_count[0:6],axis=0),0).astype(int) # only count upper 6 depth levels

# plt.text(0.2,0.01,'# profiles: ' + str(np.round(np.mean(offshore_count,axis=0),0).astype(int)),
#          transform = fig.transFigure,fontsize=5)
# plt.text(0.2,0.52,'# profiles: ' + str(np.round(np.mean(inshore_count[0:6],axis=0),0).astype(int)),
#          transform = fig.transFigure,fontsize=5)
# ax[1].set_xticklabels([f'{months[i]}' for i in range(12)] )#\n ({count1[i]}) \n ({count2[i]})
cols=plt.get_cmap('tab10')
ticks = np.arange(1,13)
for i, x in enumerate(ax[1].get_xticks()):
    ax[1].text(ticks[i], -0.3, f"({count1[i]})", size=5, ha='center',color=cols(0))
    ax[1].text(ticks[i], -1.5, f"({count2[i]})", size=5, ha='center',color=cols(1))
fig.subplots_adjust(bottom=0.15)     # Add space at bottom
# finished_plot(fig,'./plots/woa_sal_mean_seasonal_cycle_37-39N_depth_levels_p1_vs_p2.jpg')# 
# finished_plot(fig,'./plots/woa_sal_mean_seasonal_cycle_70-72W_depth_levels_p1_vs_p2.jpg')

# =============================================================================
#%% Mean seasonal cycle pot density in different time periods - depth levels 5m interval
# =============================================================================

varloc = 'lon'  # 'lat'
loc =   [-72,-70] #[37,39] #
# initialize empty list to count number of values that go into average
inshore_count = []
offshore_count = []

font_for_print()
fig,ax = plt.subplots(nrows=2,sharex=True,figsize=(10*cm,8*cm),sharey=True)
cols = plt.get_cmap('blues',13)
norm = colors.Normalize(vmin=0,vmax=50)
sm = plt.cm.ScalarMappable(cmap=cols, norm=norm)
sm.set_array([])
j=0

# add maps
cols2 = plt.get_cmap('tab10')
ax2 = fig.add_axes([0.85, 0.4, 0.1, 0.1])
ax2.plot(woa['lon'],woa['lat'],marker='.',linestyle='',color='gray')

ax3 = fig.add_axes([0.85, 0.85, 0.1, 0.1])
ax3.plot(woa['lon'],woa['lat'],marker='.',linestyle='',color='gray')


for i in np.arange(0,55,5):
    sal_inshore = subset_ds(botdepth=[20,200],presrange=[i,i+4],loc=loc,varloc=varloc,timerange=[1990,2005],var='pden')
    sal_inshore.mean('pres',skipna=True).groupby('time.month').mean(skipna=True).plot(ax=ax[0],label=str(i),
                                                                                      color=cols(j+1),linewidth=0.5,
                                                                                      linestyle='-')
    inshore_count.append(sal_inshore.mean('pres',skipna=True).groupby('time.month').count().interp(month=np.arange(1,13),
                                                                                     kwargs={"fill_value": 0}).values)
    ax3.plot(sal_inshore.mean('pres',skipna=True).dropna(dim='profile')['lon'],
             sal_inshore.mean('pres',skipna=True).dropna(dim='profile')['lat'],
             marker='.',linestyle='',color=cols2(0),markersize=1)
    sal_offshore = subset_ds(botdepth=[20,200],presrange=[i,i+4],loc=loc,varloc=varloc,timerange=[2010,2019],var='pden')
    sal_offshore.mean('pres',skipna=True).groupby('time.month').mean(skipna=True).plot(ax=ax[1],label=str(i),
                                                                                       color=cols(j+1),linewidth=0.5,
                                                                                       linestyle='-')
    offshore_count.append(sal_offshore.mean('pres',skipna=True).groupby('time.month').count().interp(month=np.arange(1,13),
                                                                                     kwargs={"fill_value": 0}).values)
    j=j+1
    
ax[0].set_xticks(np.arange(1,13))
labels=[month_converter(i+1) for i in range(12)]
ax[0].set_xticklabels(labels)
# ax[0].set_xlabel('')
# ax[1].set_xlabel('')
# ax[0].set_title('Inshore (0-50m)')
# ax[1].set_title('Offshore (50-150m)')
ax[0].set_title('1990-2005')
ax[1].set_title('2010-2019')
# ax[0].grid(False)
# ax[1].grid(False)

for axh in ax:
    axh.grid(False)
    axh.set_xlabel('')
    axh.set_ylim([20.5,26.7])
    # axh.axhline(10,color='k',alpha=0.3,zorder=0,linewidth=0.5)

# add colorbar
cbax = fig.add_axes([0.15, 0.17, 0.3, 0.03])
cb = plt.colorbar(sm,cax=cbax,ticks=np.arange(0,55,5),orientation='horizontal', 
                 boundaries=np.arange(2.5,57.5,5))
cb.ax.tick_params(labelsize=6)
cbax.xaxis.tick_top()

# add data to map
# ax2.plot(sal_inshore.dropna(dim='profile')['lon'],
#          sal_inshore.dropna(dim='profile')['lat'],
#          marker='.',linestyle='',color=cols(0),markersize=1)
ax2.plot(sal_offshore.dropna(dim='profile')['lon'],
          sal_offshore.dropna(dim='profile')['lat'],
          marker='.',linestyle='',color=cols2(1),markersize=1)

for axh in [ax2,ax3]:
    axh.set_xticklabels('')
    axh.set_yticklabels('')
    axh.set_xlim(-76.5,-65.8)
    axh.set_ylim(35,42.1)
    axh.grid(False)
    
    
count1 = np.round(np.mean(inshore_count[0:6],axis=0),0).astype(int) # only count upper 6 depth levels
count2 = np.round(np.mean(offshore_count[0:6],axis=0),0).astype(int) # only count upper 6 depth levels

# plt.text(0.2,0.01,'# profiles: ' + str(np.round(np.mean(offshore_count,axis=0),0).astype(int)),
#          transform = fig.transFigure,fontsize=5)
# plt.text(0.2,0.52,'# profiles: ' + str(np.round(np.mean(inshore_count[0:6],axis=0),0).astype(int)),
#          transform = fig.transFigure,fontsize=5)
# ax[1].set_xticklabels([f'{months[i]}' for i in range(12)] )#\n ({count1[i]}) \n ({count2[i]})
cols=plt.get_cmap('tab10')
ticks = np.arange(1,13)
for i, x in enumerate(ax[1].get_xticks()):
    ax[1].text(ticks[i], 19.3, f"({count1[i]})", size=5, ha='center',color=cols(0))
    ax[1].text(ticks[i], 19, f"({count2[i]})", size=5, ha='center',color=cols(1))
fig.subplots_adjust(bottom=0.15)     # Add space at bottom
# finished_plot(fig,'./plots/woa_sal_mean_seasonal_cycle_37-39N_depth_levels_p1_vs_p2.jpg')# 
finished_plot(fig,'./plots/woa_pden_mean_seasonal_cycle_70-72W_depth_levels_p1_vs_p2.jpg')

# =============================================================================
#%% Subplot - Mean seaonal profiles for different periods
# =============================================================================

varloc = {}
varloc['north'] = 'lon'
varloc['south'] =  'lat'
loc = {}
loc['north'] = [-72,-70] 
loc['south'] = [37,39] 

botdepth = {}
botdepth['inshore'] =  [0,50] 
botdepth['offshore'] = [50,150] 
presrange = [0,40]
varnames = ['sal','temp','pden']


#################
# function to plot different variables in different axes
def plot_mean_monthly_profiles(ds,ax,botdepth,presrange,varloc,var,loc):
    sal_inshore_p1 = subset_ds(ds,botdepth=botdepth,presrange=presrange,loc=loc,varloc=varloc,timerange=[1990,2005],var=var)
    sal_inshore_p1_month = sal_inshore_p1.groupby('time.month').mean(skipna=True).interp(month=np.arange(1,13))
    count1 = sal_inshore_p1.groupby('time.month').count().interp(month=np.arange(1,13))
    count1 = np.round(np.mean(count1[0:6],axis=0),0).values.astype(int) 
    
    sal_inshore_p2 = subset_ds(ds,botdepth=botdepth,presrange=presrange,loc=loc,varloc=varloc,timerange=[2010,2020],var=var)
    sal_inshore_p2_month = sal_inshore_p2.groupby('time.month').mean(skipna=True).interp(month=np.arange(1,13))
    # print(sal_inshore_p2_month)
    count2 = sal_inshore_p2.groupby('time.month').count().interp(month=np.arange(1,13))
    count2 = np.round(np.mean(count2[0:6],axis=0),0).values.astype(int) 
    
    # monthly mean value for plotting       
    mean_monthly = sal_inshore_p2.groupby('time.month').mean().mean('pres').interp(month=np.arange(1,13))
    
    
    ## plotting
    cols=plt.get_cmap('tab10',12)
    if var=='temp':
        inc=5
    elif var=='sal':
        inc=1.5
    else:
        inc=2
    offset = np.arange(0,inc*12+inc,inc)

    tick=[]
    ax.grid(False)
    for i in sal_inshore_p2_month.month:
        (sal_inshore_p1_month+offset[i-1]-mean_monthly[i-1]).sel(month=i).plot(ax=ax,y='pres',
                                                                               yincrease=False,label='1990-2005',
                                                                               color=cols(i-1),alpha=0.3)
        (sal_inshore_p2_month+offset[i-1]-mean_monthly[i-1]).sel(month=i).plot(ax=ax,y='pres',
                                                                               yincrease=False,label='2010-2019',
                                                                               color=cols(i-1))
        tick.append(offset[i-1])
 
        
    ax.set_xticks(tick)
    ax.set_xticklabels([month_converter(i) for i in np.arange(1,13)])
    ax.set_xlabel('')
    ax.set_title('')
    ax.set_ylim(40,0)
    ax.set_xlim(offset[0]-inc*0.5,offset[-1]+inc*0.1)
    ax.set_title(var)
    
    ticks = np.arange(1,13)
    for i in range(12):
        ax.text(tick[i], ax.get_yticks()[-1]+ax.get_yticks()[-1]*0.22, f"({count1[i]})", size=9, ha='center',color=cols(i),alpha=0.5)
        ax.text(tick[i], ax.get_yticks()[-1]+ax.get_yticks()[-1]*0.3, f"({count2[i]})", size=9, ha='center',color=cols(i))
    

#%%

### plotting 
font_for_pres()                     
for box in ['box1','box2','box3','box4']:
    for posshelf in ['inshore','offshore']:
        fig,axn = plt.subplots(figsize=(10,8),nrows=3)
        plt.subplots_adjust(hspace=0.4)
        for i in range(3):
            plot_mean_monthly_profiles(woa_box[box],ax=axn[i],varloc='lat',
                                        loc=[36,42],
                                        presrange=[0,50],
                                        botdepth=botdepth[posshelf],
                                        var=varnames[i])
        plt.suptitle(f"{box}-{posshelf}")
        # plt.legend(loc='upper right', bbox_to_anchor=(0.7, 1.1),label=['t','3'])
        # finished_plot(fig,f"./plots/ECOMON/woa_profiles_mean_seasonal_cycle_{poslat}_{posshelf}_p1_vs_p2.jpg")




# =============================================================================
#%% Mean seaonal profiles for different periods
# =============================================================================

varloc = 'lat'  # 'lat'
loc =   [35,46] #[37,39] #
botdepth =  [0,50] # [50,150] #
presrange = [0,30]
var='sal'


sal_inshore_p1 = subset_ds(woa_box['box3'],botdepth=botdepth,presrange=presrange,loc=loc,varloc=varloc,timerange=[1990,2005],var=var)
sal_inshore_p1_month = sal_inshore_p1.groupby('time.month').mean(skipna=True)
count1 = sal_inshore_p1.groupby('time.month').count()
count1 = np.round(np.mean(count1[0:6],axis=0),0).values.astype(int) 

sal_inshore_p2 = subset_ds(woa_box['box3'],botdepth=botdepth,presrange=presrange,loc=loc,varloc=varloc,timerange=[2010,2020],var=var)
sal_inshore_p2_month = sal_inshore_p2.groupby('time.month').mean(skipna=True)
count2 = sal_inshore_p2.groupby('time.month').count()
count2 = np.round(np.mean(count2[0:6],axis=0),0).values.astype(int) 

mean_monthly = sal_inshore_p2.groupby('time.month').mean().mean('pres')


### plotting 
font_for_pres()
                        
fig,ax = plt.subplots(figsize=(10,4))
cols=plt.get_cmap('tab10')
j=0
tick=[]
ax.grid(False)
for i in np.arange(1,12):
    (sal_inshore_p1_month+j-mean_monthly[i-1]).sel(month=i).plot(y='pres',yincrease=False,label='1990-2005',color='gray')
    (sal_inshore_p2_month+j-mean_monthly[i-1]).sel(month=i).plot(y='pres',yincrease=False,label='2010-2019',color=cols(0))
    tick.append(j)
    j=j+2
    
ax.set_xticks(tick)
ax.set_xticklabels([month_converter(i) for i in np.arange(1,12)])
ax.set_xlabel('')
ax.set_title('')
ax.set_ylim(30,0)
# plt.legend()

cols=plt.get_cmap('tab10')
ticks = np.arange(1,12)
for i, x in enumerate(ax.get_xticks()):
    ax.text(tick[i], ax.get_yticks()[-1]+ax.get_yticks()[-1]*0.11, f"({count1[i]})", size=9, ha='center',color='gray')
    ax.text(tick[i], ax.get_yticks()[-1]+ax.get_yticks()[-1]*0.15, f"({count2[i]})", size=9, ha='center',color=cols(0))
fig.subplots_adjust(bottom=0.2)     # Add space at bottom

# finished_plot(fig,f"./plots/woa_{var}_profiles_mean_seasonal_cycle_{loc[0]}-{loc[1]}W_offshore_p1_vs_p2.jpg")


# =============================================================================
#%% monthly salinity distribution
# =============================================================================

font_for_pres()
extent = [-77,-69,35.5,43]
fig,ax = plt.subplots(figsize=(8,10),nrows=4,ncols=3,subplot_kw = dict(projection=proj),
                      sharex=True,sharey=True,constrained_layout=True)
axh=ax.flatten()

for i in range(12):
    gl = plot_map_NWA(axh[i],extent,c=[45,100],plotbathy='contour',gulfstream=False)
    dummy = woa['sal'].sel(pres=slice(0,5)).mean('pres',skipna=True).where(woa['time.month']==i+1,drop=True)
    cc=axh[i].scatter(dummy.lon,dummy.lat,c=dummy,transform=ccrs.PlateCarree(),cmap='Spectral_r',
                   vmin=31,vmax=35,s=2,zorder=10)
    axh[i].add_feature(cartopy.feature.RIVERS)
    axh[i].text(-76,41,month_converter(i+1),transform=ccrs.PlateCarree())
    # gl.xlocator = mticker.FixedLocator(np.arange(-80,60))
    # gl.ylocator = mticker.FixedLocator(np.arange(35,44))
# plt.legend()
fig.colorbar(cc, ax=ax[:],shrink=0.8)

finished_plot(fig,'./plots/ECOMON/woa_monthly_map_scatter_salinity.jpg')


#################
# Density
#
fig,ax = plt.subplots(figsize=(8,10),nrows=4,ncols=3,subplot_kw = dict(projection=proj),
                      sharex=True,sharey=True,constrained_layout=True)
axh=ax.flatten()

for i in range(12):
    gl = plot_map_NWA(axh[i],extent,c=[45,100],plotbathy='contour',gulfstream=False)
    dummy = woa['pden'].sel(pres=slice(0,5)).mean('pres',skipna=True).where(woa['time.month']==i+1,drop=True)
    cc=axh[i].scatter(dummy.lon,dummy.lat,c=dummy,transform=ccrs.PlateCarree(),cmap='Spectral_r',s=2,zorder=10)
    axh[i].add_feature(cartopy.feature.RIVERS)
    axh[i].text(-76,41,month_converter(i+1),transform=ccrs.PlateCarree())
    # gl.xlocator = mticker.FixedLocator(np.arange(-80,60))
    # gl.ylocator = mticker.FixedLocator(np.arange(35,44))
# plt.legend()
# fig.colorbar(cc, ax=ax[:],shrink=0.8)

# finished_plot(fig,'./plots/ECOMON/woa_monthly_map_scatter_density.jpg')

#%% Test
# =============================================================================
import seaborn as sns
g = sns.JointGrid(data=test,
                  x='sal',y='pres',hue='pres',palette=cmocean.cm.dense,marginal_ticks=True,ratio=1)
g.plot_joint(sns.scatterplot,legend=False)
g.plot_marginals(sns.histplot,bins=30)
# sns.histplot(x=CT_contr.isel(profile=slice(0,100)).to_dataset(name='sal').to_dataframe()['sal'],
#              ax=g.ax_marg_x,bins=30)
g.fig.axes[0].invert_yaxis()
g.fig.axes[0].grid(False)
g.fig.axes[1].grid(False)
g.fig.axes[2].grid(False)
# g.fig.axes[0].axvline(0)

# =============================================================================
#%%   Check out density/salinity in 2012
# =============================================================================

#%%% Boxplot
plt.close('all')
presrange = [0,30]
loc = [35,45] # first number has to be smaller than second
varloc = 'lat'
var='sal'
# fig,ax = plt.subplots(figsize=(8,8),sharex=True,nrows=4,sharey=True)
for i,name in zip(range(4),['box1','box2','box3','box4']):
    dummy = woa_box[name].where((woa_box[name].time>= np.datetime64('2019-01-01')) &
                                  (woa_box[name].time<np.datetime64('2019-12-31'))
                                  ,drop=True).sel(pres=slice(*presrange)).mean('pres',skipna=True)

    
    fig,ax = plt.subplots(figsize=(8,4))
    plot_boxplot(ax,
                 subset_ds(ds=woa_box[name],botdepth=[0,100],
                                presrange=presrange,loc=loc,varloc=varloc,
                                timerange=[1990,2020],var=var),
                 label=['1990-2020','2012'],
                 var=var)
    for i in np.arange(1,13):
        dummy_plot = dummy[var].where(dummy.time.dt.month==i,drop=True).values.flatten()
        ax.plot(np.ones(len(dummy_plot))*i,dummy_plot,
                marker='*',color='r',linestyle='None',zorder=100)

    ax.set_title(f"{name}:{presrange[0]}-{presrange[1]}m")
    # finished_plot(fig,f'./plots/ECOMON/box_plot_sal_{name}_{presrange[0]}-{presrange[1]}m_all_vs_2011.jpg')
    
#%%% Minipap
plt.close('all')
def mini_map():
    plt.rcParams.update({'font.size': 11})
    fig,ax = plt.subplots(figsize=(5,3),sharex=True,sharey=True,
                              subplot_kw = dict(projection=proj_rot),constrained_layout=False)
    # plt.subplots_adjust(left=None, right=None, wspace=0.01, hspace=0.01) 
    extent = [-77,-69.5,37,41.5]
    gl = plot_map_NWA(ax,extent,
                      c=[50,75,100], #np.arange(0,200,10),#
                      plotbathy='contour',gulfstream=False) 
    return fig,ax
    
# mini_map()

# load mursst data to plot on minimaps
mursst,mursst_anomaly = load_mursst()

## sea surface salinity data preprocesses for three products
sss_monthly = np.load('./data/sss_monthly.npy',allow_pickle=True)
#%%% plot individual profiles month by month

for year in ['2011','2012'][:1]:
    for month in np.arange(1,13)[:1]:
        plt.close('all')
        font_for_pres()
        # plt.close('all')
        fig,ax = plt.subplots(figsize=(15,8),ncols=5,nrows=4,sharex=False,sharey=True,
                              constrained_layout=True)
        axh=ax.flatten()
        cols=plt.get_cmap('tab10')
        
        # plot minimap with profile positions
        # font_for_print()
        fig2,ax2=mini_map()
        (mursst_anomaly.sel(time=f'{int(year)}-{month}').squeeze()).plot.contourf(levels=np.arange(-4,4.2,0.2),ax=ax2,
                                                                                  transform=ccrs.PlateCarree(),
                                                                                  vmin=-4,vmax=4,cmap='RdBu_r')
        # fig3,ax3=mini_map()
        # (sss_monthly['smos'].sel(time=f'{int(year)}-{month}').squeeze()).plot.contourf(levels=30,ax=ax3,transform=ccrs.PlateCarree(),cmap='Spectral_r',vmin=31.5,vmax=36)
        # fig4,ax4=mini_map()
        # (sss_monthly['oisss'].sel(time=f'{int(year)}-{month}').squeeze()).plot.contourf(levels=30,ax=ax4,transform=ccrs.PlateCarree(),cmap='Spectral_r',vmin=31.5,vmax=36)
        
        # ax.grid(False)
        for name,j in zip(['box1','box2','box3','box4'],range(4)):
            # name='box4'
            maxdepth = 70
            dummy = woa_box[name].where((woa_box[name].time>= np.datetime64(f'{year}-01-01')) &
                                          (woa_box[name].time<np.datetime64(f'{year}-12-31')) &
                                          (woa_box[name].botdepth<maxdepth)
                                          ,drop=True).isel(pres=slice(0,80))
            dummy_all = woa_box[name].where((woa_box[name].time>= np.datetime64('1990-01-01')) &
                                          (woa_box[name].time<np.datetime64('2021-12-31')) &
                                          (woa_box[name].botdepth<maxdepth)
                                          ,drop=True).isel(pres=slice(0,80))
            
    
            # chose month
            # month = 1
            dummy_mm = dummy.where(dummy.time.dt.month==month,drop=True)
            # add profiles to map
            ax2.plot(dummy_mm.lon,dummy_mm.lat,transform=ccrs.PlateCarree(),
                     marker='.',linestyle='None',markersize=5,color=cols(j))
            # # add profiles to map
            # ax3.plot(dummy_mm.lon,dummy_mm.lat,transform=ccrs.PlateCarree(),
            #          marker='.',linestyle='None',markersize=5,color=cols(j))
            # ax4.plot(dummy_mm.lon,dummy_mm.lat,transform=ccrs.PlateCarree(),
            #          marker='.',linestyle='None',markersize=5,color=cols(j))
            
            # check if profiles available for given month
            if np.shape(dummy_mm['temp'].values)[1]>0: # only plot of profiles present, to avoid errors
                dummy_mm_all = dummy_all.where(dummy_all.time.dt.month==month,drop=True)
                
                #loop through variables
                for var,i in zip(['sal','temp','pden','SA_contr','CT_contr'],range(5)):
                # dummy_mm_all[var].T.plot(hue='profile',color='gray')
                    # ax[i,j].set_color_cycle(plt.get_cmap('Spectral',12))
                    dummy_mm_all[var].T.plot(ax=ax[j,i],hue='profile',y='pres',color='gray',add_legend=False,alpha=0.5,linewidth=0.5)
                    dummy_mm[var].T.plot(ax=ax[j,i],hue='profile',y='pres',add_legend=False,alpha=0.5)
                    dummy_mm_all[var].median('profile',skipna=True).T.plot(ax=ax[j,i],y='pres',color='k')
                    dummy_mm[var].mean('profile',skipna=True).T.plot(ax=ax[j,i],y='pres',color='r')
                    std = dummy_mm_all[var].std('profile',skipna=True)
                    (dummy_mm_all[var].median('profile',skipna=True)+std).T.plot(ax=ax[j,i],y='pres',color='k',
                                                                                alpha=0.5,linestyle='dashed',
                                                                                linewidth=0.1)
                    (dummy_mm_all[var].median('profile',skipna=True)-std).T.plot(ax=ax[j,i],y='pres',color='k',
                                                                                alpha=0.5,linestyle='dashed',
                                                                                linewidth=0.1)
                    ax[j,i].spines[['right', 'top']].set_visible(False)
                    
        # random things
        axh[-1].invert_yaxis()
        [axh.grid(False) for axh in ax.flatten()]
        fig.suptitle(year+'-'+month_converter(month))
        
        # # save profiles
        # finished_plot(fig,f'./plots/ECOMON/profiles/profiles_'+month_converter(month)+'_'+year+'_boxes.jpg')
        # save maps
        # finished_plot(fig2,f'./plots/ECOMON/profiles/map_MURSSTa_profiles_'+month_converter(month)+'_'+year+'.jpg')
        # finished_plot(fig3,f'./plots/ECOMON/profiles/map_SMOS_profiles_'+month_converter(month)+'_'+year+'.jpg')
        # finished_plot(fig4,f'./plots/ECOMON/profiles/map_OISSS_profiles_'+month_converter(month)+'_'+year+'.jpg')

#%%%% targeted plots for OSM2024
year = 2012
month = 5
name = 'box4'

# select data
maxdepth = 70
dummy = woa_box[name].where((woa_box[name].time>= np.datetime64(f'{year}-0{month}-01')) &
                              (woa_box[name].time<=np.datetime64(f'{year}-0{month}-31')) &
                              (woa_box[name].botdepth<maxdepth)
                              ,drop=True).isel(pres=slice(0,80))
dummy_all = woa_box[name].where((woa_box[name].time>= np.datetime64('1990-01-01')) &
                              (woa_box[name].time<np.datetime64('2021-12-31')) &
                              (woa_box[name].botdepth<maxdepth)
                              ,drop=True).isel(pres=slice(0,80))

font_for_pres()
##### temp, sal and dens ###################
fig,ax = plt.subplots(figsize=(8,3),ncols=3,constrained_layout=True,sharey=False)       
# check if profiles available for given month
if np.shape(dummy['temp'].values)[1]>0: # only plot of profiles present, to avoid errors
    dummy_mm_all = dummy_all.where(dummy_all.time.dt.month==month,drop=True)
    
#loop through variables
for var,i in zip(['sal','temp','pden'],range(3)): #,'SA_contr','CT_contr'
# dummy_mm_all[var].T.plot(hue='profile',color='gray')
    # ax[i,j].set_color_cycle(plt.get_cmap('Spectral',12))
    dummy_mm_all[var].T.plot(ax=ax[i],hue='profile',y='pres',color='gray',add_legend=False,alpha=0.1,linewidth=0.5)
    dummy[var].T.plot(ax=ax[i],hue='profile',y='pres',add_legend=False,alpha=0.5)
    dummy_mm_all[var].median('profile',skipna=True).T.plot(ax=ax[i],y='pres',color='k')
    dummy[var].mean('profile',skipna=True).T.plot(ax=ax[i],y='pres',color='r')
    std = dummy_mm_all[var].std('profile',skipna=True)
    (dummy_mm_all[var].median('profile',skipna=True)+std).T.plot(ax=ax[i],y='pres',color='k',
                                                                alpha=0.5,linestyle='dashed',
                                                                linewidth=0.1)
    (dummy_mm_all[var].median('profile',skipna=True)-std).T.plot(ax=ax[i],y='pres',color='k',
                                                                alpha=0.5,linestyle='dashed',
                                                                linewidth=0.1)
    ax[i].spines[['right', 'top']].set_visible(False)
                
    # random things
    # ax[-1].invert_yaxis()
    [axh.grid(False) for axh in ax.flatten()]
    [axh.invert_yaxis() for axh in ax.flatten()]

    ax[0].set_xlabel('salinity')
    ax[1].set_xlabel('temperature [$^\circ$C]')
    ax[2].set_xlabel('pot. density [kg/m$^3$]')
    ax[0].set_ylabel('pressure [db]')
    ax[0].set_xticks([31,32,33,34])
    for axh in ax[1::]:
        axh.spines[['left']].set_visible(False)
        axh.set_ylabel('')
        axh.set_yticks([])
        
finished_plot(fig,f'./plots/ECOMON/profiles/ecomon_profiles_05_2012_T_S_rho.jpg')

        
        
# ##### T & S controbutions  ###################
# fig,ax = plt.subplots(figsize=(6,3),ncols=2,constrained_layout=True,sharey=False)       
# # check if profiles available for given month
# if np.shape(dummy['temp'].values)[1]>0: # only plot of profiles present, to avoid errors
#     dummy_mm_all = dummy_all.where(dummy_all.time.dt.month==month,drop=True)
    
# #loop through variables
# for var,i in zip(['SA_contr','CT_contr'],range(3)): #,
# # dummy_mm_all[var].T.plot(hue='profile',color='gray')
#     # ax[i,j].set_color_cycle(plt.get_cmap('Spectral',12))
#     dummy_mm_all[var].T.plot(ax=ax[i],hue='profile',y='pres',color='gray',add_legend=False,alpha=0.1,linewidth=0.5)
#     dummy[var].T.plot(ax=ax[i],hue='profile',y='pres',add_legend=False,alpha=0.5)
#     dummy_mm_all[var].median('profile',skipna=True).T.plot(ax=ax[i],y='pres',color='k')
#     dummy[var].mean('profile',skipna=True).T.plot(ax=ax[i],y='pres',color='r')
#     std = dummy_mm_all[var].std('profile',skipna=True)
#     (dummy_mm_all[var].median('profile',skipna=True)+std).T.plot(ax=ax[i],y='pres',color='k',
#                                                                 alpha=0.5,linestyle='dashed',
#                                                                 linewidth=0.1)
#     (dummy_mm_all[var].median('profile',skipna=True)-std).T.plot(ax=ax[i],y='pres',color='k',
#                                                                 alpha=0.5,linestyle='dashed',
#                                                                 linewidth=0.1)
#     ax[i].spines[['right', 'top']].set_visible(False)
                
    # # random things
    # # ax[-1].invert_yaxis()
    # [axh.grid(False) for axh in ax.flatten()]
    # [axh.invert_yaxis() for axh in ax.flatten()]
    
    # ax[0].set_xlabel('salinity')
    # ax[1].set_xlabel('temperature [$^\circ$C]')
    # ax[0].set_ylabel('pressure [db]')
    # for axh in ax[1::]:
    #     axh.spines[['left']].set_visible(False)
    #     axh.set_ylabel('')
    #     axh.set_yticks([])
    # fig.suptitle(year+'-'+month_converter(month))
#%% box plot for profiles
def plot_boxplot_profile(ax,ds,label,var='sal'):
    import matplotlib.patches as mpatches
    cols = plt.get_cmap('tab10')
    font_for_pres()
    # dummy = ds.where(ds.time.dt.month==i).mean(dim='pres',skipna=True).stack(z=['profile'])
    for i in range(10):
        dummy = ds.mean(dim='pres',skipna=True).stack(z=['profile'])

        if i==1:
            pdbox = pd.DataFrame({f'{i}': dummy})
            # pdbox2 = pd.DataFrame({f'{i}': dummy})
        else:
            pdbox[f'{i}'] = dummy
        
        
    if var=='sal':
        ylim=[29,35.8]
        ylabel='salinity'
    elif var=='temp':
        ylim=[3,28] 
        ylabel='temperature [\N{DEGREE SIGN}C]'
        
    # plotting
    bp_dict = pdbox.plot(kind='box',ylim=ylim,ax=ax,xticks=np.arange(1,13),
                         patch_artist=True,return_type='both',sym='.')
    months = [month_converter(i+1) for i in range(12)]
    count1 = pdbox.count().values
    ticks = np.arange(1,13)

    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{months[i]}' for i in range(12)])#\n ({count1[i]}) \n ({count2[i]})
    for i in np.arange(0,12):
        ax.text(ticks[i], 28.2, f"({count1[i]})", size=7, ha='center',color=cols(0))
    fig.subplots_adjust(bottom=0.15)     # Add space at bottom
    
    # for row_key, (ax,row) in bp_dict.iteritems():
    for i in np.arange(0,len(bp_dict[1]['boxes'])):
        bp_dict[1]['boxes'][i].set_facecolor(cols(1))
        # bp_dict[1]['fliers'][i].set_color('blue')
        bp_dict[1]['medians'][i].set_color(cols(0))
    # for i in np.arange(4,len(bp_dict[1]['whiskers']),6):
    #     bp_dict[1]['fliers'][i].set_color('lightgray')
    # bp_dict[1]['fliers'][0].set_color('gray')
    ax.set_ylabel(ylabel)
    
    # #add legend
    # patch1 = mpatches.Patch(color=cols(0),label=label[0])
    # # patch2 = mpatches.Patch(color=cols(1),label=label[1])
    # plt.legend(handles=patch1,loc='lower left',fontsize=6)

    #add map
    ax2 = ax.inset_axes([0.8, 0.1, 0.15, 0.15])
    ax2.plot(woa['lon'],woa['lat'],marker='.',linestyle='',color='gray')
    ax2.plot(ds.dropna(dim='profile')['lon'],ds.dropna(dim='profile')['lat'],
              marker='.',linestyle='',color=cols(3),markersize=1)
    ax2.set_xticks([-75,-70])
    ax2.set_xticklabels(['75W','70W'],fontsize=7)
    ax2.set_yticklabels(['35N','40N'],fontsize=7)
    ax2.minorticks_off()
    ax2.set_xlim(-76.5,-65.8)
    ax2.set_ylim(35,42.1)
    return fig,ax,pdbox

#%%

name='box3'
dummy_2011 = woa_box[name].where((woa_box[name].time>= np.datetime64('2011-01-01')) &
                              (woa_box[name].time<np.datetime64('2011-12-31')) &
                              (woa_box[name].botdepth<100)
                              ,drop=True).isel(pres=slice(0,50))
dummy_all = woa_box[name].where((woa_box[name].time>= np.datetime64('1990-01-01')) &
                              (woa_box[name].time<np.datetime64('2022-12-31')) &
                              (woa_box[name].botdepth<100)
                              ,drop=True).isel(pres=slice(0,50))
dummy_mm = dummy_2011.where(dummy_2011.time.dt.month<6,drop=True)
dummy_mm_all = dummy_all.where(dummy_all.time.dt.month==10,drop=True)


cols = plt.get_cmap('tab10')
font_for_pres()
# dummy = ds.where(ds.time.dt.month==i).mean(dim='pres',skipna=True).stack(z=['profile'])
fig,ax = plt.subplots()
ds = dummy_mm_all
for i in np.arange(0,100,10):
    dummy = ds['sal'].sel(pres=slice(i,i+10)).mean(dim='pres',skipna=True).stack(z=['profile'])
    if i==0:
        pdbox = pd.DataFrame({f'{i}': dummy})
        # pdbox2 = pd.DataFrame({f'{i}': dummy})
    else:
        pdbox[f'{i}'] = dummy
# pdbox = pdbox.T
# pdbox['index'] = ['0','10','20','30','40','50','60','70','80','90']
# pdbox = pdbox.set_index('index')
  # plotting
bp_dict = pdbox.plot(kind='box',ax=ax,
                     patch_artist=True,return_type='both',sym='.',vert=False)
ax.invert_yaxis()

j=1
for i in np.arange(0,100,10):
    dd = dummy_mm['sal'].sel(pres=slice(i,i+10)).mean(dim='pres',skipna=True)
    ax.plot(dd,np.ones(len(dd))*j,
            marker='*',color='r',linestyle='None',zorder=100)
    j=j+1