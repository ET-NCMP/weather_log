import argparse
import os
import numpy as np
from netCDF4 import Dataset, num2date
from datetime import datetime
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib import rcParams, colors, cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


DATADIR = os.path.join(os.environ['DATADIR'], 'WeatherLog')


class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def get_climate_files(month):
    print("Getting climate files for month", month)
    moidx = month - 1  # pythonic

    # MSLP
    f1 = Dataset(os.path.join(DATADIR, "mslp_climate_era5_1991_2020.nc"), 'r')
    mslp_climate = f1.variables['mslp_climate'][moidx]
    f1.close()

    # Z500 and thickness
    f2 = Dataset(os.path.join(DATADIR, "z500_thickness_climate_era5_1991_2020.nc"), 'r')
    z500_climate = f2.variables['z500_climate'][moidx]
    thickness_climate = f2.variables['thickness_climate'][moidx]
    f2.close()

    return mslp_climate, z500_climate, thickness_climate


def get_current_data(year, month):
    print("Getting data for month", year, month)
#    month = str(month).zfill(2)

    # MSLP
    f1 = Dataset(os.path.join(DATADIR, f"era5_monthly_weatherlog_mslp_{year}_{month:02d}.nc"), 'r')
    lats = f1.variables['latitude'][:]
    lons = f1.variables['longitude'][:]
    mslp = f1.variables['msl'][0] / 100
    f1.close()

    # Z500 and thickness
    f2 = Dataset(os.path.join(DATADIR, f"era5_monthly_weatherlog_plev_{year}_{month:02d}.nc"), 'r')
    levs = f2.variables['level'][:]
    z250 = f2.variables['z'][0, np.where(levs == 250)[0][0]] / 9.81
    z500 = f2.variables['z'][0, np.where(levs == 500)[0][0]] / 9.81
    z1000 = f2.variables['z'][0, np.where(levs == 1000)[0][0]] / 9.81
    u250 = f2.variables['u'][0, np.where(levs == 250)[0][0]]
    v250 = f2.variables['v'][0, np.where(levs == 250)[0][0]]
    f2.close()
    thickness = z500 - z1000
    wind_250 = np.sqrt(u250 ** 2 + v250 ** 2)

    return mslp, z250, z500, thickness, wind_250, lons, lats


def plot_map(year_plot, month_plot):

    dpi = 600

    mslp_climate, z500_climate, thickness_climate = get_climate_files(month_plot)
    mslp_month, z250_month, z500_month, thickness_month, wind_month, lons, lats = get_current_data(year_plot,
                                                                                                   month_plot)

    mslp_anom = mslp_month - mslp_climate
    thickness_anom = thickness_month - thickness_climate
    z500_anom = z500_month - z500_climate

    # PLot!
    proj = ccrs.NorthPolarStereo()
    trans = ccrs.PlateCarree()

    # contouring levels
    # MSLP actual
    msl_actual_levs = np.arange(880, 1084, 4)
    # MSLP anom
    # msl_anom_levs = np.concatenate((np.arange(-20,0,2),np.arange(2,22,2)))
    msl_anom_levs = np.arange(-20, 22, 2)
    msl_cb_ticks = np.arange(-20, 24, 4)
    # Thickness actual
    thickness_actual_levs = np.arange(492, 600, 6)
    thick_cb_ticks = np.arange(-14, 16, 4)
    # Wind
    wind_levs = np.arange(20, 82, 2)
    wind_cb_ticks = np.arange(20, 90, 10)
    z250_clevs = np.arange(900, 1110, 12)
    # Thickness anom
    thickness_clevs = np.concatenate((np.arange(-14, 0, 2), np.arange(2, 16, 2)))
    # thickness_clevs=np.arange(-14,16,2)
    # boundary for map
    extent = [-180, 180, 30, 90]

    # Lines for dateline/Greenwich Meridian
    dateline_lat, dateline_lon = 30, -180
    pole_lat, pole_lon = 30, 0
    geodetic = ccrs.Geodetic()
    dateline_lon_t, dateline_lat_t = proj.transform_point(dateline_lon, dateline_lat, geodetic)
    pole_lon_t, pole_lat_t = proj.transform_point(pole_lon, pole_lat, geodetic)
    # Line for 90W-90E
    west_lat, west_lon = 30, -90
    east_lat, east_lon = 30, 90
    geodetic = ccrs.Geodetic()
    west_lon_t, west_lat_t = proj.transform_point(west_lon, west_lat, geodetic)
    east_lon_t, east_lat_t = proj.transform_point(east_lon, east_lat, geodetic)

    # linestyle for these
    ls = '-'
    lw = 0.4
    c = 'gray'
    fontsize_title = 6.5
    # colourmaps for the three subplots
    cmap_u = 'RdYlBu_r'
    nlevs = len(msl_anom_levs)
    cmap_thick = plt.get_cmap('RdBu_r')
    # create colourmap for MSL
    bottom = cm.get_cmap('YlOrRd', 22)  # 128
    top = cm.get_cmap('Blues_r', 22)

    newcolors = np.vstack((top(np.linspace(0, 1, 22)), bottom(np.linspace(0, 1, 22))))
    white = np.array([1.0, 1.0, 1.0, 1.0])
    newcolors[20:24, :] = white
    cmap_msl = ListedColormap(newcolors, name='BuYlOrRd')

    print("Begin figure")
    fig, axs = plt.subplots(3, 1, tight_layout=True, figsize=(3, 9), subplot_kw=dict(projection=proj))

    # Z250
    cf3 = axs[0].contourf(lons, lats, wind_month, cmap=cmap_u, transform=trans, levels=wind_levs, extend='max')
    th_cr = axs[0].contour(lons, lats, z250_month / 10, colors='k', linewidths=0.25, levels=z250_clevs, transform=trans)
    axs[0].clabel(th_cr, fmt='%1.0f', fontsize=4, inline_spacing=0.5)
    axs[0].set_extent(extent, ccrs.PlateCarree())
    axs[0].coastlines(color='gray')

    axs[0].plot([dateline_lon_t, pole_lon_t], [dateline_lat_t, pole_lat_t], color=c, linewidth=lw, linestyle=ls)
    axs[0].plot([west_lon_t, east_lon_t], [west_lat_t, east_lat_t], color=c, linewidth=lw, linestyle=ls)

    axs[0].set_title("250hPa height (contours, dam) and wind (filled)", loc='left', fontsize=fontsize_title)

    cax = fig.add_axes([0.88, 0.65, 0.03, 0.23])
    cb3 = plt.colorbar(cf3, cax=cax, spacing='proportional', orientation='vertical', ticks=wind_cb_ticks)
    cb3.set_label("m s$^{-1}$", fontsize=8)
    cb3.ax.tick_params(labelsize=7)

    # Thickness
    cf2 = axs[1].contourf(lons, lats, thickness_anom / 10, cmap=cmap_thick, transform=trans, levels=thickness_clevs,
                          extend='both')
    z_cr = axs[1].contour(lons, lats, thickness_month / 10, colors='k', linewidths=0.25, levels=thickness_actual_levs,
                          transform=trans)
    axs[1].clabel(z_cr, fmt='%1.0f', fontsize=4, inline_spacing=0.5)
    axs[1].set_extent(extent, ccrs.PlateCarree())
    axs[1].coastlines(color='gray')

    axs[1].set_title("1000-500hPa thickness (contours) and anomalies (filled)", loc='left', fontsize=fontsize_title)
    axs[1].plot([dateline_lon_t, pole_lon_t], [dateline_lat_t, pole_lat_t], color=c, linewidth=lw, linestyle=ls)
    axs[1].plot([west_lon_t, east_lon_t], [west_lat_t, east_lat_t], color=c, linewidth=lw, linestyle=ls)

    cax = fig.add_axes([0.88, 0.38, 0.03, 0.23])
    cb2 = plt.colorbar(cf2, cax=cax, spacing='proportional', orientation='vertical', ticks=thick_cb_ticks)
    cb2.set_label("dam", fontsize=8)
    cb2.ax.tick_params(labelsize=7)

    # MSLP
    cf1 = axs[2].contourf(lons, lats, mslp_anom, cmap=cmap_msl, transform=trans, levels=msl_anom_levs, extend='both')
    # cmap_msl._lut[int((nlevs/2)-1):int((nlevs/2)+1)] = [1.,1.,1.,1.]
    msl_cr = axs[2].contour(lons, lats, mslp_month, colors='k', linewidths=0.25, levels=msl_actual_levs,
                            transform=trans)
    axs[2].clabel(msl_cr, fmt='%1.0f', fontsize=4, inline_spacing=0.5)
    axs[2].set_extent(extent, ccrs.PlateCarree())
    axs[2].coastlines(color='gray')

    axs[2].set_title("Mean MSLP (contours) and anomalies (filled)", loc='left', fontsize=fontsize_title)
    axs[2].plot([dateline_lon_t, pole_lon_t], [dateline_lat_t, pole_lat_t], color=c, linewidth=lw, linestyle=ls)
    axs[2].plot([west_lon_t, east_lon_t], [west_lat_t, east_lat_t], color=c, linewidth=lw, linestyle=ls)

    cax = fig.add_axes([0.88, 0.11, 0.03, 0.23])
    cb = plt.colorbar(cf1, cax=cax, spacing='proportional', orientation='vertical', ticks=msl_cb_ticks)
    cb.set_label("hPa", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    figfile = os.path.join(DATADIR, f"weather_log_front_page_{year_plot}_{month_plot:02d}.png")
    plt.savefig(figfile, dpi=dpi, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Weather Log Map Plotter')
    parser.add_argument('-y', type=int, default=2021, help='Year to plot')
    parser.add_argument('-m', type=int, default=6, help='Month to plot')

    args = parser.parse_args()

    plot_map(args.y, args.m)
