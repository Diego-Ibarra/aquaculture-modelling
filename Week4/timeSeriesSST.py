# Imports for the entire module
import pandas as pd

def get_data(lat,lon,minday,maxday):
    isub = 1

    base_url='http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdAGsstamday_LonPM180.csv?'
    query='sst[('+minday+'):'+str(isub)+':('+maxday+')][(0.0):'+str(isub)+':(0.0)][('+str(lat)+'):'+str(isub)+':('+str(lat)+')][('+str(lon)+'):'+str(isub)+':('+str(lon)+')]'
    url = base_url+query
    print('Downloading Satellite Data from: POES, AVHRR and GAC (Time-series)...')

    data = pd.read_csv(url)
    sst = data['sst'][1:].astype(float)
    return sst

def make_forcing(days,dt,sst):
    # Resize array and Eliminate Gaps by Interpolating in between
    import scipy.interpolate as intrp
    import numpy as np

    months = len(sst)

    x = np.arange(0, months)

    f = intrp.interp1d(x, sst, kind='linear', fill_value='extrapolate' )

    NoSTEPS = int(days / dt)
    newx = np.linspace(0,days/30,NoSTEPS) # Makes and vector array of equally spaced numbers from zero to "days"

    new_sst = f(newx)

    # Get rid of nans
    nans, x= np.isnan(new_sst), lambda z: z.nonzero()[0]
    new_sst[nans]= np.interp(x(nans), x(~nans), new_sst[~nans])
    return new_sst

def plot_data(sst,title='SST'):
    import matplotlib.pyplot as plt
    fig, (ax) = plt.subplots(1,1)
    ax.plot(sst,'b-')
    #ax = sst.plot()
    ax.set_xlabel('Number of records')
    ax.set_ylabel('Sea Surface Temperature (oC)')
    ax.set_title(title)
    return

          
          

if __name__ == "__main__":

    lat = 43
    lon = -62
    minday='2010-01-01T12:00:00Z'
    maxday='2012-12-31T12:00:00Z'
    days = 365 * 3
    dt = 0.01
    
    sst = get_data(lat,lon,minday,maxday)
    new_sst = make_forcing(days,dt,sst)

    plot_data(sst, title='Orginal sst')

    plot_data(new_sst, title='New INTERPOLATED sst')
    