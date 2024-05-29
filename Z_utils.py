import numpy as np
from scipy.optimize import curve_fit

def exponential_func(x, a, b):
    return a * np.exp(-b * x)

def exponential_cdf_func(x, a):
     return 1-np.exp(-a*x)

def get_binding_data(res, dt):
    dt_sum = np.cumsum(dt)
    a=np.where(res==0)[0]
    Dt = dt_sum[a]
    Diff = Dt[1:]-Dt[:-1]
    return Diff

def get_stepping_data(res, dt):
    dt_sum = np.cumsum(dt)
    a=np.where(res==2)[0]
    Dt = dt_sum[a]
    Diff = Dt[1:]-Dt[:-1]
    return Diff

def normalised_cdf_histogramm(data, bins =100, density=True):
    counts, bin_edges = np.histogram(data, bins=bins, density=True)
    cdf = np.cumsum(counts)
    M = cdf[-1]
    return bin_edges[1:], cdf/M

def normalised_pdf_histogramm(data, bins =100, density=True):
    counts, bin_edges = np.histogram(data, bins=bins, density=True)
    S = np.sum(counts)    
    return (bin_edges[:-1] + bin_edges[1:]) / 2, counts/S



def fit (res, dt, func):   
    # hist, bins = np.histogram(func(res, dt), bins=1000, density=True) 
    # bin_centers = (bins[:-1] + bins[1:]) / 2

    x,y = normalised_pdf_histogramm(func(res, dt))
    params, covariance = curve_fit(exponential_func, x,y)
    return params[1], covariance[1][0]

def fit_cdf(res, dt, func):
    x,y = normalised_cdf_histogramm(func(res, dt), bins=1000)
    params, covariance = curve_fit(exponential_cdf_func, x,y)
    return params[0], covariance[0][0]