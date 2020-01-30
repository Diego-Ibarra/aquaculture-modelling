import urllib
import netCDF4 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW LINE !!!!! <<<<<<

year = 2016
month = 5
minlat = 38
maxlat = 48
minlon = -67
maxlon = -50
isub = 1 # Stride

minday = str(year)+'-'+str(month).zfill(2)+'-16T12:00:00Z'
maxday = str(year)+'-'+str(month).zfill(2)+'-16T12:00:00Z'

base_url='https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdAGsstamday_LonPM180.nc?'
query='sst[('+minday+'):'+str(isub)+':('+maxday+')][(0.0):'+str(isub)+':(0.0)][('+str(minlat)+'):'+str(isub)+':('+str(maxlat)+')][('+str(minlon)+'):'+str(isub)+':('+str(maxlon)+')]'
URL = base_url+query

URL2 = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdAGsstamday_LonPM180.nc?sst[(2016-05-16T12:00:00Z):1:(2016-05-16T12:00:00Z)][(0.0):1:(0.0)][(38):1:(48)][(-67):1:(-50)]'

# Download data and store it in NetCDF file
fileName='sst2.nc' 
urllib.urlretrieve (URL, fileName)

#%% open NetCDF data from file # <<<<<<<<<<<<<<<<<<< NEW LINE !!!!! <<<<<<
nc = netCDF4.Dataset(fileName)   # <<<<<<<<<<<<<<<<< NEW LINE !!!!! <<<<<<
ncv = nc.variables # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW LINE !!!!! <<<<<<