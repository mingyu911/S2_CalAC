"""
FileName: map_rgb.py
Description: Generates RGB quick images.
Author: Mingyu Shin (mingyu911@pukyong.ac.kr)
"""
import math
import os

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

# Map axes position in figure coordinates: [left, bottom, width, height]
AX_RECT = [0.12, 0.08, 0.85, 0.80]

GADM_ROOT = "/home/mingyu911/S2_CalAC/data/GADM"
GADM_COUNTRIES = ['KOR', 'JPN', 'PRK', 'CHN', 'RUS', 'TWN']

GAMMA = 0.5               # RGB gamma stretch
SHOW_BOUNDARIES = False   # True: draw GADM country boundaries

PLATE = ccrs.PlateCarree()


def apply_research_style():
    """Sets a professional style for academic plots."""
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'DejaVu Serif'],
        'mathtext.fontset': 'cm',
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'axes.linewidth': 1.2,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'savefig.dpi': 300,
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


def add_boundaries(ax, edgecolor='black', linewidth=0.5, level=0):
    """Draws GADM country boundaries (only the shapefiles that exist)."""
    for code in GADM_COUNTRIES:
        shp_path = os.path.join(GADM_ROOT, f"gadm41_{code}_shp", f"gadm41_{code}_{level}.shp")
        if os.path.exists(shp_path):
            ax.add_geometries(shpreader.Reader(shp_path).geometries(), PLATE,
                              facecolor='none', edgecolor=edgecolor, linewidth=linewidth)


def draw_rgb_map(lon, lat, lon_range, lat_range, refl_r, refl_g, refl_b, map_name, file_name, save_path):
    apply_research_style()

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_axes(AX_RECT, projection=PLATE)
    ax.set_extent([*lon_range, *lat_range], crs=PLATE)
    ax.set_aspect('auto')

    # RGB image: clip to [0, 1] and gamma-stretch
    rgb = np.clip(np.dstack([refl_r, refl_g, refl_b]), 0, 1) ** GAMMA
    rgb = rgb[::-1, :, :]  # same row flip as the original index_map[::-1, :]

    ax.pcolormesh(lon, lat, rgb, transform=PLATE, shading='nearest')

    if SHOW_BOUNDARIES:
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

    # Axis labels and title, centered on the map axes
    left, bottom, width, height = AX_RECT
    x_center = left + width / 2
    y_center = bottom + height / 2
    fig.text(x_center, 0.01, "Longitude (°)", ha='center', va='bottom', fontsize=17)
    fig.text(0.02, y_center, "Latitude (°)", ha='center', va='center', rotation='vertical', fontsize=17)
    fig.text(x_center, 0.92, map_name, ha='center', va='bottom', fontsize=20)

    os.makedirs(save_path, exist_ok=True)
    fig.savefig(os.path.join(save_path, f"{file_name}.png"))
    plt.close(fig)