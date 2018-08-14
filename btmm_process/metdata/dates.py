import numpy as np
import pandas as pd
from datetime import timedelta
import datetime


def matlabdn2datetime(matlabDatenum):
    # Convert MATLAB datenum to numpy datetime 64
    npDatetime = [np.datetime64(datetime.fromordinal(int(dt))
                                + timedelta(days=float(np.remainder(dt, 1)))
                                - timedelta(days=366))
                  for dt in matlabDatenum]
    return(npDatetime)


def dparse(y, mon, d, h, minute, s, ms):
    # Function to parse the date column within read_csv
    x = y + '/' + mon + '/' + d + ' ' + h + ':' + minute + ':' + s + '.' + ms
    d = pd.datetime.strptime(x, '%Y/%m/%d %H:%M:%S.%f')
    return d
