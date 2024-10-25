import os
#os.environ["DISPLAY"] = "localhost:10.0"
# 14, 16,  work
#standart libraries:
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
import pandas as pd
import string
import sys
import imp
import xarray as xr
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import scipy.io as sc
import datetime as dt
import warnings 
warnings.filterwarnings("ignore") 
import cmocean as cm
import scipy.signal as signal
import seawater as sw
import gsw
import proplot as pplt
import matplotlib.dates as mdates
from scipy.interpolate import griddata
import mat73
# import scipy as sc


#xr.set_options(display_width=80, display_style='text')
# xr.set_options(display_style='text')

## conversion for figures sizes from inches to cm
cm=1/2.54

## projection for Northwest Atlantic
# proj = ccrs.LambertConformal(central_latitude = 40, 
#                               central_longitude = -71, 
#                               standard_parallels = (25, 25))

## projection for Northwest Atlantic
proj = ccrs.LambertConformal(central_latitude = 40, 
                              central_longitude = -65, 
                              standard_parallels = (25, 25))


## load bathymetry
bathy = xr.open_dataset('/vast/clidex/data/bathymetry/ETOPO2v2/ETOPO2v2c_f4.nc').sel(x=slice(-85,-40),y=slice(25,55)).load()

## pioneer locations
pioneer_loc = {}
pioneer_loc['inshore'] = [40.367,-70.8817]
pioneer_loc['ups_inshore'] = [40.367,-70.7842]
pioneer_loc['centr_inshore'] = [40.2250,-70.89]
pioneer_loc['central'] = [40.13325,-70.77832]
pioneer_loc['centr_offshore'] = [40.0980,-70.8804]
pioneer_loc['offshore'] = [39.9350,-70.8858]
pioneer_loc['ups_offshore'] = [39.9397,-70.7705]

# =============================================================================
# Some functions
# =============================================================================

## projection for Northwest Atlantic rotated
proj_rot = ccrs.LambertConformal(central_latitude = 40, 
                              central_longitude = 0, 
                              standard_parallels = (25, 25))
#
#
## map preparation 
def fig_map(ax,extent,figsize=(5,5)):
    """
    Produces map with coastlines using cartopy.
    
    INPUT:
    ax       : axis handle 
    extent   : [lon1,lon2,lat1,lat2]
    figsize  : specify size (default (5,5))
    
    OUTPUT:
    fig      : Figure handle
    ax       : Axis handle
    
    """
    # modules
    import cartopy.crs as ccrs
    import cartopy
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.ticker as mticker
    import cartopy.feature

    rivers = cartopy.feature.NaturalEarthFeature(
        category='physical', name='rivers_lake_centerlines',
        scale='10m', facecolor='none', edgecolor='cornflowerblue')

    
    # projection
    proj = ccrs.PlateCarree()

    # formatting
    ax.set_extent(extent)
#     ax.coastlines(resolution='50m',color='gray')
    gl = ax.gridlines(crs=proj,draw_labels=True,linestyle='-',alpha=0.1,linewidth=0.5,color='gray',
                      xlocs=range(-120,40,5),ylocs=range(0,80,5),x_inline=False)
    ax.add_feature(cartopy.feature.GSHHSFeature(scale='intermediate'), facecolor='lightgray',edgecolor='gray')
    # ax.add_feature(cartopy.feature.RIVERS)
    ax.add_feature(rivers, linewidth=1)
    gl.xlocator = mticker.FixedLocator(np.arange(-80,-40,5))
    gl.ylocator = mticker.FixedLocator([35, 37.5, 40, 42.5, 45])
    gl.rotate_labels = False
    
    gl.ylabels_right = False
    gl.xlabels_top = False
    return gl
#
#
## plot manual clabels
def manual_clabel(x,y,ax,cc,proj):
    for ind in range(len(x)):
        xy_trans = proj.transform_point(x[ind],y[ind],ccrs.PlateCarree())
        manual_location = [(xy_trans[0],xy_trans[1])]
        ax.clabel(cc,fontsize=10,fmt='%1dm',inline=1, manual=manual_location)
#
#        
## cut time series to only include full years (important when deriving mean seasonal cycle)
def cut2fullyears(ds):
    years = ds['time.year']
    year1 = np.min(years)
    yearend = np.max(years)
    if len(ds['time'].where(years==year1,drop=True))<365:
        year1 = year1+1
    if len(ds['time'].where(years==yearend,drop=True))<365:
        yearend=yearend-1
#         print(f'yearend: {yearend}')
    return ds.sel(time=slice(f'{year1.values}-01',f'{yearend.values}-01'))
#
#
## convert numeric month to string
# def month_converter(month):
#     months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
#     return months[month-1]
def month_converter(i):
    import calendar 
    # get month name
    # return calendar.month_name[i]
    return calendar.month_abbr[i]

# function for monthly ticks
def ticks_month(ax):
    ax.set_xticks(np.arange(1,13))
    ax.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']);
    ax.set_xlabel('')

#
#
#
# =============================================================================
# Gulf Stream data
# =============================================================================
gs = sc.loadmat('/mnt/data/GS_monthly_yearly_two_sat_contours_1993_to_2018.mat')
gs_time = gs['time_monthly']
gs_lon = gs['lon_cell'][0,0][0,:]
gs_lat = gs['lat_cell'][0,0][0,:]
#
#
##  derive mean seasonal cycle to plot
# 1) convert
def matlab2datetime(matlab_datenum):
    day = dt.datetime.fromordinal(int(matlab_datenum))
    dayfrac = dt.timedelta(days=matlab_datenum%1) - dt.timedelta(days = 366)
    return day + dayfrac

# apply function to every element
convert = np.vectorize(matlab2datetime)
gs_datetime = convert(gs['time_monthly'][0])
#
#
# 2) have to interpolate onto regular lon grid, in order to derive a mean positions. 
# This is probelmatic when large meanders are present, however in the western portion 
# the GS seems to be stable enough so that it doesn't make a difference
lon_int = np.arange(-78,-50,0.1)+360
lat_int=[]

for i in range(len(gs['lat_cell_monthly'])):
    lon_monthly = gs['lon_cell_monthly'][i][0][0]
    lat_monthly = gs['lat_cell_monthly'][i][0][0]
    lat_int.append(np.interp(lon_int,lon_monthly,lat_monthly))
#
#
# 3) create xarray
da = xr.DataArray(
    data=lat_int,
    dims=["time", "lon"],
    coords=dict(
        lon=(["lon"], lon_int-360),
        time=gs_datetime))
#
#
# 4) derive monthly mean
gs_lat_seas = cut2fullyears(da.sel(time=slice('2015-01','2021-12'))).groupby('time.month').mean('time')

# # my own libraries:
# #import m_general as M

# #import AA_plot_base as AA
# def json_load(name, path, verbose=False):
#     import json
#     full_name= (os.path.join(path,name+ '.json'))

#     with open(full_name, 'r') as ifile:
#         data=json.load(ifile)
#     if verbose:
#         print('loaded from: ',full_name)
#     return data

# # mconfig=json_load('config','/home/mhell/2021_ICESat2_tracks/config/')

# # # add project depenent libraries
# # sys.path.append(mconfig['paths']['local_script'])
# # sys.path.append(mconfig['paths']['local_script'] +'/ICEsat2_SI_tools/')

# # import m_colormanager_ph3 as M_color
# # import m_tools_ph3 as MT
# # import m_general_ph3 as M

# #load colorscheme
# # col=M_color.color(path=mconfig['paths']['analysis']+'../config/', name='color_def')


# lstrings =iter([i+') ' for i in list(string.ascii_lowercase)])
# # define journal fig sizes
# # fig_sizes = mconfig['fig_sizes']['AMS']


SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
# #csfont = {'fontname':'Comic Sans MS'}
# # legend_properties = {'weight':'bold'}
# #font.family: sans-serif
# #font.sans-serif: Helvetica Neue

# #import matplotlib.font_manager as font_manager
# #font_dirs = ['/home/mhell/HelveticaNeue/', ]
# #font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
# #font_list = font_manager.createFontList(font_files)
# #font_manager.fontManager.ttflist.extend(font_list)

plt.rc('font', size=SMALL_SIZE, serif='Helvetica Neue', weight='normal')          # controls default text sizes
# #plt.rc('font', size=SMALL_SIZE, serif='DejaVu Sans', weight='light')
# plt.rc('text', usetex='false')
# plt.rc('axes', titlesize=MEDIUM_SIZE, labelweight='normal')     # fontsize of the axes title
# plt.rc('axes', labelsize=SMALL_SIZE, labelweight='normal') #, family='bold')    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE, frameon=False)    # legend fontsize
# plt.rc('figure', titlesize=MEDIUM_SIZE, titleweight='bold', autolayout=True) #, family='bold')  # fontsize of the figure title
plt.rc('figure', facecolor='white') # white background on figure

# #figure.autolayout : False
# #matplotlib.rcParams['pdf.fonttype'] = 42
# #matplotlib.rcParams['ps.fonttype'] = 42


# plt.rc('path', simplify=True)

# plt.rcParams['figure.figsize'] = (8, 4)#(20.0, 10.0) #inline


# ### TICKS
# # see http://matplotlib.org/api/axis_api.html#matplotlib.axis.Tick
# #xtick.top            : False   # draw ticks on the top side
# #xtick.bottom         : True   # draw ticks on the bottom side
# #xtick.major.size     : 3.5      # major tick size in points
# #xtick.minor.size     : 2      # minor tick size in points
# #xtick.major.width    : .8    # major tick width in points
# #xtick.minor.width    : 0.6    # minor tick width in points
# #xtick.major.pad      : 3.5      # distance to major tick label in points
# #xtick.minor.pad      : 3.4      # distance to the minor tick label in points
# #xtick.color          : k      # color of the tick labels
# #xtick.labelsize      : medium # fontsize of the tick labels
# #xtick.direction      : out    # direction: in, out, or inout
# #xtick.minor.visible  : False  # visibility of minor ticks on x-axis
# #xtick.major.top      : True   # draw x axis top major ticks
# #xtick.major.bottom   : True   # draw x axis bottom major ticks
# #xtick.minor.top      : True   # draw x axis top minor ticks
# #xtick.minor.bottom   : True   # draw x axis bottom minor ticks

# #ytick.left           : True   # draw ticks on the left side
# #ytick.right          : False  # draw ticks on the right side
# #ytick.major.size     : 3.5      # major tick size in points
# #ytick.minor.size     : 2      # minor tick size in points
# #ytick.major.width    : 0.8    # major tick width in points
# #ytick.minor.width    : 0.6    # minor tick width in points
# #ytick.major.pad      : 3.5      # distance to major tick label in points
# #ytick.minor.pad      : 3.4      # distance to the minor tick label in points
# #ytick.color          : k      # color of the tick labels
# #ytick.labelsize      : medium # fontsize of the tick labels
# #ytick.direction      : out    # direction: in, out, or inout
# #ytick.minor.visible  : False  # visibility of minor ticks on y-axis
# #ytick.major.left     : True   # draw y axis left major ticks
# #ytick.major.right    : True   # draw y axis right major ticks
# #ytick.minor.left     : True   # draw y axis left minor ticks
# #ytick.minor.right    : True   # draw y axis right minor ticks
# plt.tick_params(axis='x', which='both', direction='out', top=True, bottom=True,length=3)
# 

plt.rc('xtick.minor',top=False,bottom=False)
plt.rc('ytick.minor',left=False,right=False)

# plt.rc('xtick.major', size= 4, width=1 )
# plt.rc('ytick.major', size= 3.8, width=1 )

# #axes.facecolor      : white   # axes background color
# #axes.edgecolor      : black   # axes edge color
# #axes.linewidth      : 0.8     # edge linewidth
# #axes.grid           : False   # display grid or not
# #axes.titlesize      : large   # fontsize of the axes title
# #axes.titlepad       : 6.0     # pad between axes and title in points
# #axes.labelsize      : medium  # fontsize of the x any y labels
# #axes.labelpad       : 4.0     # space between label and axis
# #axes.labelweight    : normal  # weight of the x and y labels
# #axes.labelcolor     : black

# plt.rc('axes', labelsize= MEDIUM_SIZE, labelweight='normal')




# # axes.spines.left   : True   # display axis spines
# # axes.spines.bottom : True
# # axes.spines.top    : True
# # axes.spines.right  : True
# plt.rc('axes.spines', top= False, right=False )



def font_for_print():

    SMALL_SIZE = 6
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 10
    #csfont = {'fontname':'Comic Sans MS'}
    legend_properties = {'weight':'bold'}
    #font.family: sans-serif
    #font.sans-serif: Helvetica Neue

    #import matplotlib.font_manager as font_manager
    #font_dirs = ['/home/mhell/HelveticaNeue/', ]
    #font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
    #font_list = font_manager.createFontList(font_files)
    #font_manager.fontManager.ttflist.extend(font_list)

    plt.rc('font', size=SMALL_SIZE, serif='Helvetica Neue', weight='normal')          # controls default text sizes
    #plt.rc('font', size=SMALL_SIZE, serif='DejaVu Sans', weight='light')
    plt.rc('text', usetex='false')
    plt.rc('axes', titlesize=MEDIUM_SIZE, labelweight='normal')     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE, labelweight='normal') #, family='bold')    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE, frameon=False)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE, titleweight='bold', autolayout=False) #, family='bold')  # fontsize of the figure title

    #figure.autolayout : False
    #matplotlib.rcParams['pdf.fonttype'] = 42
    #matplotlib.rcParams['ps.fonttype'] = 42


    #plt.rc('xtick.major', size= 4, width=1 )
    #plt.rc('ytick.major', size= 3.8, width=1 )


    plt.rc('axes', labelsize= SMALL_SIZE, labelweight='normal')
    
def font_medium():

    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12
    #csfont = {'fontname':'Comic Sans MS'}
    legend_properties = {'weight':'bold'}
    #font.family: sans-serif
    #font.sans-serif: Helvetica Neue

    #import matplotlib.font_manager as font_manager
    #font_dirs = ['/home/mhell/HelveticaNeue/', ]
    #font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
    #font_list = font_manager.createFontList(font_files)
    #font_manager.fontManager.ttflist.extend(font_list)

    plt.rc('font', size=SMALL_SIZE, serif='Helvetica Neue', weight='normal')          # controls default text sizes
    #plt.rc('font', size=SMALL_SIZE, serif='DejaVu Sans', weight='light')
    plt.rc('text', usetex='false')
    plt.rc('axes', titlesize=MEDIUM_SIZE, labelweight='normal')     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE, labelweight='normal') #, family='bold')    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE, frameon=False)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE, titleweight='bold', autolayout=False) #, family='bold')  # fontsize of the figure title

    #figure.autolayout : False
    #matplotlib.rcParams['pdf.fonttype'] = 42
    #matplotlib.rcParams['ps.fonttype'] = 42


    #plt.rc('xtick.major', size= 4, width=1 )
    #plt.rc('ytick.major', size= 3.8, width=1 )


    plt.rc('axes', labelsize= SMALL_SIZE, labelweight='normal')
    
    

def font_for_pres():

    SMALL_SIZE = 10
    MEDIUM_SIZE = 12
    BIGGER_SIZE = 14
    #csfont = {'fontname':'Comic Sans MS'}
    legend_properties = {'weight':'bold'}
    #font.family: sans-serif
    #font.sans-serif: Helvetica Neue

    #import matplotlib.font_manager as font_manager
    #font_dirs = ['/home/mhell/HelveticaNeue/', ]
    #font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
    #font_list = font_manager.createFontList(font_files)
    #font_manager.fontManager.ttflist.extend(font_list)

    plt.rc('font', size=SMALL_SIZE, serif='Helvetica Neue', weight='normal')          # controls default text sizes
    #plt.rc('font', size=SMALL_SIZE, serif='DejaVu Sans', weight='light')
    plt.rc('text', usetex='false')
    plt.rc('axes', titlesize=SMALL_SIZE, labelweight='normal')     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE, labelweight='normal') #, family='bold')    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE, frameon=False)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE, titleweight='bold', autolayout=False) #, family='bold')  # fontsize of the figure title

    #figure.autolayout : False
    #matplotlib.rcParams['pdf.fonttype'] = 42
    #matplotlib.rcParams['ps.fonttype'] = 42


    #plt.rc('xtick.major', size= 4, width=1 )
    #plt.rc('ytick.major', size= 3.8, width=1 )


    plt.rc('axes', labelsize= SMALL_SIZE, labelweight='normal')
    
    
def font_medium():

    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 12
    #csfont = {'fontname':'Comic Sans MS'}
    legend_properties = {'weight':'bold'}
    #font.family: sans-serif
    #font.sans-serif: Helvetica Neue

    #import matplotlib.font_manager as font_manager
    #font_dirs = ['/home/mhell/HelveticaNeue/', ]
    #font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
    #font_list = font_manager.createFontList(font_files)
    #font_manager.fontManager.ttflist.extend(font_list)

    plt.rc('font', size=SMALL_SIZE, serif='Helvetica Neue', weight='normal')          # controls default text sizes
    #plt.rc('font', size=SMALL_SIZE, serif='DejaVu Sans', weight='light')
    plt.rc('text', usetex='false')
    plt.rc('axes', titlesize=MEDIUM_SIZE, labelweight='normal')     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE, labelweight='normal') #, family='bold')    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE, frameon=False)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE, titleweight='bold', autolayout=False) #,



# add project depenent libraries
#sys.path.append(config['paths']['local_script'])
plt.tick_params(axis='x', which='both', direction='out', top=True, bottom=True,length=3)
