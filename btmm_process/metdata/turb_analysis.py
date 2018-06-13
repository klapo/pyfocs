import xarray as xr

def quad_analy(quant1, quant2):
    # Upper right (+ and +)
    quad1 = ((quant1 > 0) & (quant2 > 0)).mean() * 100
    # Upper left (- and +)
    quad2 = ((quant1 < 0) & (quant2 > 0)).mean() * 100
    # Lower right (+ and -)
    quad3 = ((quant1 > 0) & (quant2 < 0)).mean() * 100
    # Lower left (- and -)
    quad4 = ((quant1 < 0) & (quant2 < 0)).mean() * 100

    return([quad1, quad2, quad3, quad4])

def reynolds_decomp(xr_ds, field_name):
    halfhour_index = (xr_ds['time.hour'] / 0.5)
    halfhour_index.name = 'halfHourly'

    bar = xr_ds[field_name].groupby(halfhour_index).mean()
    out_data_array = xr_ds[field_name].groupby(halfhour_index) - bar

    return(out_data_array)

def standardize(x):
    return (x - x.mean()) / x.std()
