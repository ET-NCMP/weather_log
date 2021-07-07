#!/usr/bin/env python3
import cdsapi
import numpy as np
import os

## User specify year & month to download
year='2021'
month='05'

##=============================================================================
os.chdir("/storage/basic/arise/fh004579/weather_log/current_month")

c = cdsapi.Client()

# MSLP (single level data), Z500 (pressure level), Z1000 (pressure level)

## Pressure level data
c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'product_type': 'monthly_averaged_reanalysis',
        'variable': ['geopotential','u_component_of_wind','v_component_of_wind'],
        'pressure_level': [
            '250','500', '1000',
        ],
        'year': str(year),
        'month': str(month),
        'time': '00:00',
        'area': [
            90, -180, 0,
            180,
        ],
        'format': 'netcdf',
    },
    'era5_monthly_weatherlog_plev_'+str(year)+'_'+str(month)+'.nc')

## Single level data
c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': 'mean_sea_level_pressure',
        'year': str(year),
        'month': str(month),
        'time': '00:00',
        'area': [
            90, -180, 0,
            180,
        ],
    },
    'era5_monthly_weatherlog_mslp_'+str(year)+'_'+str(month)+'.nc')
