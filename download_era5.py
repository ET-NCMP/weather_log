#!/usr/bin/env python3
import argparse
import cdsapi
import os


DATADIR = os.path.join(os.environ['DATADIR'], 'WeatherLog')


def download_era(year, month):

    c = cdsapi.Client()

    # MSLP (single level data), Z500 (pressure level), Z1000 (pressure level)
    plev_file = os.path.join(DATADIR, f'era5_monthly_weatherlog_plev_{year}_{month:02d}.nc')
    mslp_file = os.path.join(DATADIR, f'era5_monthly_weatherlog_mslp_{year}_{month:02d}.nc')

    # Pressure level data
    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'product_type': 'monthly_averaged_reanalysis',
            'variable': ['geopotential', 'u_component_of_wind', 'v_component_of_wind'],
            'pressure_level': [
                '250', '500', '1000',
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
        plev_file)

    # Single level data
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
        mslp_file)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Weather Log Data Getter')
    parser.add_argument('-y', type=int, default=2021, help='Year to plot')
    parser.add_argument('-m', type=int, default=6, help='Month to plot')

    args = parser.parse_args()

    download_era(args.y, args.m)
