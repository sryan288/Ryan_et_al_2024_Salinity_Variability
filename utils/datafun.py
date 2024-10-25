'''
 Data manipulation utilities
'''
'''
    1) butterworth_lowpass_filter          Low-pass filter along given axis
    2) matlab2datetime                     Convert Matlab time (datenum) to datetime
    3) datetime2matlab                     Convert datetime to Matlab time (datenum)
    4) datetim2numeric                     convert mdatetime64 to numeric value to perform arithmetrics
    5) columndate2datetime                 convert array with columns of yy,mm,dd,HH,MM,SS to datetime object

'''

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from scipy import signal
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
import xarray as xr
import pandas as pd
import datetime as dt



#
#------------------------------------------------------------------------
# 1) butterworth_lowpass_filter
       
def butterworth_lowpass_filter(data, order=2, cutoff_freq=1.0/10.0, axis=0):
    """Filter input data.
    
    For unfiltered data, use `cutoff_freq=1`.
    
    Currently, this returns a numpy array.
    """
    B, A = signal.butter(order, cutoff_freq, output="ba")
    return signal.filtfilt(B, A, data, axis=axis)

#
#------------------------------------------------------------------------
# 2) convert matlab time to datetime
def matlab2datetime(matlab_datenum):
    day = dt.datetime.fromordinal(int(matlab_datenum))
    dayfrac = dt.timedelta(days=matlab_datenum%1) - dt.timedelta(days = 366)
    return day + dayfrac

#
#------------------------------------------------------------------------
# 3) convert datetime to Matlab's datenum
def datetime2matlab(time):
    """
    INPUT
    time: timevectors as datetime64 fromat
    
    OUTPUT
    timenum: time in Matlab's datenum format
    """
    timenum = mdates.date2num(time)+datetime.date(1971, 1, 2).toordinal()
    return timenum


#
#------------------------------------------------------------------------
# 4) convert mdatetime64 to numeric value to perform arithmetrics
def datetime2numeric(time):
    # convert xarray datetime to pandas datime (still a string)
    time_pd = pd.to_datetime(time.values)
    # convert to int64 (numeric value - not sure yet what exactly the numbers mean, can to simple arithmetric like datenum in matlba)
    time_np = time_pd.astype(np.int64)
    average_time_np = np.average(time_np)
    # convert bacl to datetime64
    average_time_pd = np.datetime64(pd.to_datetime(average_time_np))
    
#
#
# 5) convert column dates to datetime object
def columndate2datetime(date_array):
    """
    date_array = [[2011,05,11,12,03,30],
                  [2012,05,11,12,03,30],
                  [2011,06,11,12,03,30]]
    
    returns array with datetime objects
    """
    return [datetime(year, month, day,HH,MM,SS) for year, month, day,HH,MM,SS in date_array]