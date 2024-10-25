'''
 Plotting utilities for the Northwest Atlantic
'''
'''
    1) plot_map          Plot bathymetry map of NWA
    2) finished_plot     Either show plot or save if name is provided
    3) calculate_ticks   Align ticks for twinx or twiny
    4) find_bounds       Find bounds for vmin,vmax based in max,min values in dataset
    5) anomaly           Plot anomlay timeseries with anomalies shaded
    6) find_common_max   For comparison find common values for vmin and vmax so that colormap is equal
    7) plot_ipo_bar      Plot bar indicating IPO phases on top of panel
    8) latlon_labels     Turn ticklabels to geocoordinates
    9) ts_appen          Append TS diagrams with density contours and axis labels
    10) datetime_ticks   Set ticks on datetime axis
    11) print_path       Print file path on figure
'''

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import matplotlib.dates as dates
import matplotlib.dates as mdates
import xarray as xr
import seawater as sw
import scipy.io as sc
from utils.xarray_utilities import cut2fullyears

# load gulf stream data
da = xr.open_dataset('/home/sryan/python/utils/gulfstream_position_25cm_xarray.nc')['__xarray_dataarray_variable__']



#
#------------------------------------------------------------------------
# 1) plot bathymetry map

def plot_map_NWA(ax,extent,c=np.arange(0,6000,500),plotbathy=None,gulfstream=False,gscolor='k'):
    """
    Produces map with coastlines using cartopy.
    
    INPUT:
    ax       : axis handle 
    extent   : [lon1,lon2,lat1,lat2]
    
    OUTPUT:
    gl       : grid handle
    
    """
    
    ## load bathymetry
    bathy = xr.open_dataset('/vast/clidex/data/bathymetry/ETOPO2v2/ETOPO2v2c_f4.nc').sel(x=slice(-85,-40),y=slice(25,55)).load()

    
    # projection
    proj = ccrs.PlateCarree()
    
    ## plot bathymetry
    if plotbathy=='contourf':
        cc=ax.contourf(bathy.x,
                    bathy.y,
                    (bathy.z*(-1)),
                    levels=c,cmap=plt.get_cmap('cmo.deep',len(c)),
                    transform=ccrs.PlateCarree())
        plt.colorbar(cc,shrink=0.8)
    elif plotbathy=='contour':
        ax.contour(bathy.x,bathy.y,(bathy.z*(-1)),
                   levels=c,colors='k',transform=ccrs.PlateCarree(),
                   linewidths=0.5,
                   alpha=0.5)
        
    ## add mean Gulf Stream path
    if gulfstream:
        gs = sc.loadmat('/mnt/data/GS_monthly_yearly_two_sat_contours_1993_to_2018.mat')
        gs_time = gs['time_monthly']
        gs_lon = gs['lon_cell'][0,0][0,:]
        gs_lat = gs['lat_cell'][0,0][0,:]
        plt.plot(gs_lon,gs_lat,
                  transform=ccrs.PlateCarree(),color=gscolor,linewidth=1)
        # gs_lat_seas = cut2fullyears(da.sel(time=slice('2010-01','2021-12'))).mean('time')
        # gs_lat_seas.plot(color='k',ax=ax,transform=ccrs.PlateCarree(),linewidth=1);

    ## formatting
    ax.set_extent(extent)
#     ax.coastlines(resolution='50m',color='gray')
    gl = ax.gridlines(crs=proj,draw_labels=True,linestyle='-',alpha=0.1,
                      xlocs=range(-120,40,5),ylocs=range(0,80,5),x_inline=False)
    ax.add_feature(cartopy.feature.GSHHSFeature(scale='intermediate'),
                   facecolor='lightgray',edgecolor='None')
#     gl.xlocator = mticker.FixedLocator(np.arange(-80,-40,5))
#     gl.ylocator = mticker.FixedLocator([10, 20, 30, 40, 50])
    gl.rotate_labels = False
    
    gl.ylabels_right = False
    gl.xlabels_top = False
    
    return gl
    

#
#------------------------------------------------------------------------
# 2) If a figure name is defined, save the figure to that file. Otherwise, display the figure on screen.
def finished_plot (fig, fig_name=None, dpi=300, printpath=None):
    
    if printpath:
        print_path(printpath,fig)
    
    if fig_name is not None:
        print('Saving ' + fig_name)
        fig.savefig(fig_name, dpi=dpi,bbox_inches='tight')
    else:
        fig.show()
        
        
#
#
#-------------------------------------------------------------------------
# 3) align ticks for twin axes
def calculate_ticks(ax, ticks,axis = 'x', round_to=0.1, center=False):
    if axis=='x':
        upperbound = np.ceil(ax.get_xbound()[1]/round_to)
        lowerbound = np.floor(ax.get_xbound()[0]/round_to)
    elif axis=='y':
        upperbound = np.ceil(ax.get_ybound()[1]/round_to)
        lowerbound = np.floor(ax.get_ybound()[0]/round_to)
    dy = upperbound - lowerbound
    fit = np.floor(dy/(ticks - 1)) + 1
    dy_new = (ticks - 1)*fit
    if center:
        offset = np.floor((dy_new - dy)/2)
        lowerbound = lowerbound - offset
    values = np.linspace(lowerbound, lowerbound + dy_new, ticks)
    return values*round_to


#
#
#-------------------------------------------------------------------------
# 4) Find bounds for vmin,vmax based in max,min values in dataset
def find_bounds(ds,decimals,anomaly=False):
    vmin = np.round(ds.min(skipna=True).values- 0.5 * 10**(-decimals),
                    decimals=int(decimals)) # work around, as floor&ceil don't take decimals
    vmax = np.round(ds.max(skipna=True).values+ 0.5 * 10**(-decimals),
                    decimals=int(decimals))
    if anomaly is True:
        vmax = np.max(np.abs((vmin,vmax)))
        vmin = vmax*(-1)
    return vmin,vmax


#--------------------------------------------------------------------------
# 5) anomaly plot
def anomaly(ax,x,y,offset):
    ax.fill_between(x, 0+offset[0], y, where = y>0+offset[0],
                    facecolor='indianred', interpolate=True)
    ax.fill_between(x, 0-offset[1], y, where = y<0-offset[1],
                    facecolor='royalblue', interpolate=True)

#
#
#--------------------------------------------------------------------------
# 6) for comparison find common values for vmin and vmax so that colormap is equal
def find_common_cmax(data1,data2):    
    vmin1 = np.round(np.min((data1.fillna(0).values,data1.fillna(0).values)),decimals=2)
    vmax1 = np.round(np.max((data1.fillna(-100000).values,data1.fillna(-100000).values)),decimals=2)
    vmin2 = np.round(np.min((data2.fillna(0).values,data2.fillna(0).values)),decimals=2)
    vmax2 = np.round(np.max((data2.fillna(-100000).values,data2.fillna(-100000).values)),decimals=2)    
    vval = np.max(np.abs((vmin1,vmax1,vmin2,vmax2)))
    return vval

#
#
#
#------------------------------------------------------------------------
# 7) add ipo phase as bar
# Create rectangle x coordinates
def add_ipo_bar(ax):
    import matplotlib.dates as mdates
    indices = xr.open_dataset('/vortex/clidex/data/obs/climate_indices/indices_noaa_psl_May_13_2020.nc')
    cols = ['dodgerblue','indianred','dodgerblue','indianred']
    text = ['-IPO','+IPO','-IPO','']
    years = [1956,1977,1999,2013,2020]
    for i in range(len(years)-1):
        startTime = indices.sel(Month=slice(str(years[i]) + '-01-01',str(years[i]) + '-01-31'))['Month'].values
        endTime =  indices.sel(Month=slice(str(years[i+1]-1) + '-12-01',str(years[i+1]-1) + '-12-31'))['Month'].values#startTime + timedelta(seconds = 1)
        # convert to matplotlib date representation
        start = mdates.date2num(startTime)
        end = mdates.date2num(endTime)
        width = end - start
        middle = (width/2)+start
        
        ulim = ax.get_ybound()[1]
        llim = ax.get_ybound()[1] - ax.get_ybound()[1]/6 
        rect = Rectangle((start[0], llim), width, ulim, color=cols[i],alpha=0.5)
        ax.text(middle,(ulim-llim)/2+llim,text[i],fontsize=8,fontweight='bold',verticalalignment='center',
               horizontalalignment='center')
        ax.add_patch(rect)

        

#------------------------------------------------------------------------
# 8) Convert tick labels to geocoordinates 
# 
def latlon_label(ax,coord='longitude',axis='x',ticks=None,fmt='%1.2f'):
    """
    Function gets ticks of current axis and creates string with geocoordinates and sets ticklabels
    
    INPUT
    ax:      axis handle
    coord:   string specifying whether latitude or longitude 
    axis:    specifiy if x or y axis
    ticks:   array of tick values if they should be change (!not working so far)
    
    """
    labels=[]
#     ax.set_xticks = ticks
#     if axis=='x': ax.set_xticks = ticks
#     if axis=='y': ax.set_yticks = ticks
    
    # get labels if not specified
    if ticks is None:
        if axis=='x':
            ticks = ax.get_xticks()
        elif axis =='y':
            ticks = ax.get_yticks() 
   
    # create geolabels
    for tick in ticks:
        if (tick>0 and coord=='longitude'):
            labels.append(f"{fmt}\N{DEGREE SIGN}E" % tick)
        elif (tick<0 and coord=='longitude'):
            labels.append(f"{fmt}\N{DEGREE SIGN}W" % abs(tick))
        elif (tick>0 and coord=='latitude'):
            labels.append(f"{fmt}\N{DEGREE SIGN}N" % tick)
        elif (tick<0 and coord=='latitude'):
            labels.append(f"{fmt}\N{DEGREE SIGN}S" % abs(tick))
            
    #set axis labels
    if axis=='x': ax.set_xticklabels(labels)
    elif axis=='y': ax.set_yticklabels(labels)
        

        
#------------------------------------------------------------------------
# 9) Append TS plot with density contours 
# 
        
def ts_append(levels,axh=None):
    """
    ts_append(levels,axh=None)
    
    function to append TS plots with density contours and axis labels
    
    INPUT:
    levels:    array of density levels that should be plotted
    axh:       axis handle if desired
    
    """
    
    # get axis limits of current plot
    if axh:
        tlim = np.ceil(axh.get_ylim())
        slim = np.ceil(axh.get_xlim())
    else: 
        tlim = np.ceil(plt.gca().get_ylim())
        slim = np.ceil(plt.gca().get_xlim())
        
    # create meshgrid and derive  potential density
    sm,tm = np.meshgrid(np.arange(slim[0]-7,slim[1]+7),np.arange(tlim[0]-7,tlim[1]+7))
    densm = sw.pden(sm,tm,0,0)
    # plot contours
    cc=axh.contour(sm,tm,densm-1000,colors='k',alpha=0.3,zorder=0,levels=levels,
                   linestyle='dashed',linewidths=0.5)
    
    # set axis limits to original
    axh.set_ylim(tlim)
    axh.set_xlim(slim)
    
    plt.clabel(cc,fmt='%2.1f',fontsize=8)

    axh.set_xlabel('salinity')
    axh.set_ylabel('temperature [°C]')
    axh.grid(False)
    
    
#------------------------------------------------------------------------
# 10) Set ticks on datetime axis 
# 
def datetime_tick(ax,axis='xaxis',interval='year'):
    eval(f"ax.{axis}.set_major_locator(mdates.YearLocator(1))")
    eval(f"ax.{axis}.set_major_formatter(mdates.DateFormatter('%Y'))")

#------------------------------------------------------------------------
# 11) print file path on figure
#
def print_path(path,fig):
    plt.text(0.01, 0.01, path, horizontalalignment='left',
             verticalalignment='center', transform=fig.transFigure,
             fontsize=5)