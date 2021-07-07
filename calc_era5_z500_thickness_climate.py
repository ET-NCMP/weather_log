import numpy as np
from netCDF4 import Dataset,num2date
from datetime import datetime
import matplotlib.pyplot as plt

def get_climate(years):
    """
    Calculates monthly Z500 and thickness climatology over specified years.
    Data are NH on 0.25 degree grid but grid can be specified.
    """
    ## Dimensions of grid
    nmonths=12; nlats= 361; nlons= 1440
    ## Array to store output from files
    z500_store = np.zeros((len(years),nmonths,nlats,nlons))
    thickness_store = np.zeros((len(years),nmonths,nlats,nlons))
    ## Loop through years and get MSLP
    for y, year in enumerate(years):
        print(year)

        f = Dataset("/storage/basic/arise/fh004579/weather_log/climate_files/era5_monthly_z500_z1000_"+str(year)+".nc",'r')

        if y==0: ## open the grid
            lats = f.variables['latitude'][:]
            lons = f.variables['longitude'][:]
            levs = f.variables['level'][:]

        z500 = f.variables['z'][:,np.where(levs==500)[0][0],:,:]/9.81
        z1000 = f.variables['z'][:,np.where(levs==1000)[0][0],:,:]/9.81

        z500_store[y]= z500
        thickness_store[y]=z500-z1000

        f.close()

    ## Calculate climate
    climate_z500 = np.mean(z500_store,axis=0)
    climate_thickness = np.mean(thickness_store,axis=0)

    return climate_z500, climate_thickness, lats, lons

years=np.arange(1981,2011,1)
climate_z500,climate_thickness, lats,lons = get_climate(years)

## time axis which is the month
time=np.arange(1,13,1)

## Now output to netCDF
# ### output to a netCDF
f = Dataset('z500_thickness_climate_era5_'+str(years[0])+'_'+str(years[-1])+'.nc','w', format='NETCDF4')
f.createDimension('latitude', len(lats))
f.createDimension('longitude', len(lons))
f.createDimension('time', len(time))

latid=f.createVariable('latitude','f4','latitude')
lonid=f.createVariable('longitude','f4','longitude')
timeid=f.createVariable('time','f4','time')
climidz=f.createVariable('z500_climate','f8',('time','latitude','longitude'))
climidthick=f.createVariable('thickness_climate','f8',('time','latitude','longitude'))

latid[:]=lats
lonid[:]=lons
timeid[:]=time
timeid.units='month of year'

climidz[:]=climate_z500
climidz.units='gpm'
climidz.long_name="Monthly mean z500 height "+str(years[0])+"-"+str(years[-1])+" ERA-5"

climidthick[:]=climate_thickness
climidthick.units='gpm'
climidthick.long_name="Monthly mean 500-1000 hPa thickness "+str(years[0])+"-"+str(years[-1])+" ERA-5"

f.close()
