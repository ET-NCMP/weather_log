import numpy as np
from netCDF4 import Dataset,num2date
from datetime import datetime
import matplotlib.pyplot as plt

def mslp_climate(years):
    """
    Calculates monthly MSLP climatology over specified years.
    Data are NH on 0.25 degree grid but grid can be specified.
    """
    ## Dimensions of grid
    nmonths=12; nlats= 361; nlons= 1440
    ## Array to store output from files
    mslp_store = np.zeros((len(years),nmonths,nlats,nlons))
    ## Loop through years and get MSLP
    for y, year in enumerate(years):
        print(year)
        f = Dataset("/storage/basic/arise/fh004579/weather_log/climate_files/era5_monthly_mslp_"+str(year)+".nc",'r')
        mslp_store[y]= f.variables['msl'][:]/100 ## convert to hPa
        if y==0: ## open the grid
            lats = f.variables['latitude'][:]
            lons = f.variables['longitude'][:]

        f.close()

    ## Calculate climate
    climate_mslp = np.mean(mslp_store,axis=0)

    return climate_mslp, lats, lons

years=np.arange(1981,2011,1)
climate_mslp, lats,lons = mslp_climate(years)

## time axis which is the month
time=np.arange(1,13,1)

## Now output to netCDF
# ### output to a netCDF
f = Dataset('mslp_climate_era5_'+str(years[0])+'_'+str(years[-1])+'.nc','w', format='NETCDF4')
f.createDimension('latitude', len(lats))
f.createDimension('longitude', len(lons))
f.createDimension('time', len(time))

latid=f.createVariable('latitude','f4','latitude')
lonid=f.createVariable('longitude','f4','longitude')
timeid=f.createVariable('time','f4','time')
climid=f.createVariable('mslp_climate','f8',('time','latitude','longitude'))

latid[:]=lats
lonid[:]=lons
timeid[:]=time
timeid.units='month of year'

climid[:]=climate_mslp
climid.units='hPa'
climid.long_name="Monthly mean MSLP"+str(years[0])+"-"+str(years[-1])+" ERA-5"

f.close()
