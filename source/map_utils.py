"""
FileName: map_utils.py
Description: Shared mapping library for research-grade satellite visualization.
Author: Mingyu Shin (mingyu911@pukyong.ac.kr)
"""
import math
import os

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

# Map / colorbar axes position in figure coordinates: [left, bottom, width, height]
AX_RECT = [0.12, 0.08, 0.85, 0.80]
CAX_RECT = [0.12, 0.91, 0.85, 0.03]

GADM_ROOT = "/home/mingyu911/S2_CalAC/data/GADM"
GADM_COUNTRIES = ['KOR', 'JPN', 'PRK', 'CHN', 'RUS', 'TWN']

PLATE = ccrs.PlateCarree()


def apply_research_style():
    """Sets a professional style for academic plots."""
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'DejaVu Serif'],
        'mathtext.fontset': 'cm',
        'axes.titlesize': 25,
        'axes.labelsize': 15,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'axes.linewidth': 1.2,
        'lines.linewidth': 1.2,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'savefig.dpi': 300,
        'savefig.format': 'png',
    })


def get_grid_step(gap):
    """Picks a 1 / 2.5 / 5 x 10^n step giving about 2-3 grid lines over `gap`."""
    target = gap / 2.5
    base = 10 ** math.floor(math.log10(target))
    rel = target / base

    if rel < 1.5:
        return 1 * base
    if rel < 3.5:
        return 2.5 * base
    return 5 * base


def get_tick_format(step):
    """Returns a number format with just enough decimals for the given tick step."""
    dec = 0
    while abs(round(step, dec) - step) > 1e-10 and dec < 6:
        dec += 1
    return f'.{dec}f'


def get_cbar_step(gap):
    """Picks a 1 / 2 / 5 / 10 x 10^n step giving about 5 colorbar intervals over `gap`."""
    target = gap / 5.0
    if target == 0:
        return 1

    base = 10 ** math.floor(math.log10(target))
    rel = target / base

    if rel < 1.5:
        return 1 * base
    if rel < 3.5:
        return 2 * base
    if rel < 7.5:
        return 5 * base
    return 10 * base


def add_boundaries(ax, edgecolor='black', linewidth=1.0, level=0):
    """Draws GADM country boundaries (only the shapefiles that exist)."""
    for code in GADM_COUNTRIES:
        shp_path = os.path.join(GADM_ROOT, f"gadm41_{code}_shp", f"gadm41_{code}_{level}.shp")
        if os.path.exists(shp_path):
            ax.add_geometries(shpreader.Reader(shp_path).geometries(), PLATE,
                              facecolor='none', edgecolor=edgecolor, linewidth=linewidth)


def draw_map(cmap, lon, lat, lon_range, lat_range, map_para, para_range, map_name, file_name, save_path):
    """Standardized mapping function called by main scripts."""
    apply_research_style()

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_axes(AX_RECT, projection=PLATE)
    ax.set_extent([*lon_range, *lat_range], crs=PLATE)
    ax.set_aspect('auto')

    # Parameter map
    norm = mcolors.Normalize(vmin=para_range[0], vmax=para_range[1])
    im = ax.pcolormesh(lon, lat, map_para, cmap=cmap, norm=norm, transform=PLATE, shading='auto')

    # Top colorbar
    cax = fig.add_axes(CAX_RECT)
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    cb.set_label(map_name, labelpad=15, fontsize=20)
    cax.xaxis.set_label_position('top')

    add_boundaries(ax)

    # Grid labels only (no grid lines)
    x_step = get_grid_step(lon_range[1] - lon_range[0])
    y_step = get_grid_step(lat_range[1] - lat_range[0])

    gl = ax.gridlines(draw_labels=True, linewidth=0)
    gl.top_labels = gl.right_labels = False
    gl.xlocator = mticker.MultipleLocator(x_step)
    gl.ylocator = mticker.MultipleLocator(y_step)
    gl.xformatter = LongitudeFormatter(number_format=get_tick_format(x_step))
    gl.yformatter = LatitudeFormatter(number_format=get_tick_format(y_step))

    # Axis labels, centered on the map axes
    left, bottom, width, height = AX_RECT
    x_center = left + width / 2
    y_center = bottom + height / 2
    fig.text(x_center, 0.01, "Longitude (°)", ha='center', va='bottom', fontsize=17)
    fig.text(0.02, y_center, "Latitude (°)", ha='center', va='center', rotation='vertical', fontsize=17)

    # Colorbar ticks
    cbar_step = get_cbar_step(para_range[1] - para_range[0])
    cax.xaxis.set_major_locator(mticker.MultipleLocator(cbar_step))
    cax.xaxis.set_minor_locator(mticker.MultipleLocator(cbar_step / 2))
    cax.tick_params(axis='x', which='major', direction='in', length=10, labelsize=15)
    cax.tick_params(axis='x', which='minor', direction='in', length=5)
    cax.xaxis.set_ticks_position('bottom')

    os.makedirs(save_path, exist_ok=True)
    fig.savefig(os.path.join(save_path, f"{file_name}.png"), dpi=300)
    plt.close(fig)