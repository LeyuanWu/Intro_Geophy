##################################################################################
# # ! IGRF-13 计算得到中国主要城市主磁场近百年变化
# # ! Plot variations of Magnetic (Z, F, D, I) from 1900 to 2025
# # ! for 6 major Chinese cities using the IGRF-13 model (2×2 layout)
##################################################################################
# %% # ! Setup
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from LW_MAG import *

# Plot style settings
from matplotlib import font_manager as fm
font_path = 'NotoSansCJKsc-Regular.otf'
fm.fontManager.addfont(font_path)
plt.rcParams['font.sans-serif'] = ['Noto Sans CJK SC', 'DejaVu Sans']
# %%
# # ! 城市
CITIES = {
    '北京'    : (39.9042, 116.4074),
    '杭州'    : (30.2741, 120.1551),
    '西安'    : (34.3416, 108.9398),
    '广州'    : (23.1291, 113.2644),
    '成都'    : (30.5728, 104.0668),
    '乌鲁木齐': (43.8256, 87.6168)
}

print(f"{'City':<12} {'Latitude':>12} {'Longitude':>12}", flush=True)
lats = []; lons = []
for city, (lat, lon) in CITIES.items():
    lats.append(lat); lons.append(lon)
    print(f'{city:12} {lat:12.4f} {lon:12.4f}', flush=True)
Lat = np.array(lats); Lon = np.array(lons)

# %%
# # ! 时间
DATES = [datetime(year, 1, 1) for year in range(1900, 2026, 5)]
print(f'\n{"Epoch":^12}', flush=True)
for date in DATES:
    str_date = f'{date:%Y-%m-%d}'
    print(f'{str_date:^12}', flush=True)

# %%
# # ! 地磁场模型计算 <偏角、倾角、总场强度>
# TODO ******************** Parameters ******************** #
h = 0.0
shc_model = 'IGRF13'
# TODO ********************** end ************************* #
n_row = len(DATES); n_col = len(CITIES)
BN = np.zeros((n_row, n_col))
BE = np.zeros((n_row, n_col))
BD = np.zeros((n_row, n_col))
for i, date in enumerate(DATES): 
    Bn, Be, Bd = SHS_gd(Lon, Lat, h, date, '%s.shc' % shc_model)
    BN[i, :] = Bn; BE[i, :] = Be; BD[i, :] = Bd
BX, BY, BZ, BD, BI, BH, BF = ComputeMagB7(BN, BE, BD)

# %%
# # ! Pandas table of cities with coordinates and magnetic elements (latest epoch)
latest_idx = -1   # 2025

# Latest values
df_latest = pd.DataFrame({
    '城市'       : list(CITIES.keys()),
    '纬度'       : np.round(lats, 4),
    '经度'       : np.round(lons, 4),
    'Z (nT)'    : np.round(BZ[latest_idx, :], 1),
    'F (nT)'    : np.round(BF[latest_idx, :], 1),
    'D (°)'     : np.round(BD[latest_idx, :], 2),
    'I (°)'     : np.round(BI[latest_idx, :], 2)
})

# Secular change 1900 -> 2025
df_change = pd.DataFrame({
    '城市'       : list(CITIES.keys()),
    'ΔZ (nT)'   : np.round(BZ[latest_idx, :] - BZ[0, :], 1),
    'ΔF (nT)'   : np.round(BF[latest_idx, :] - BF[0, :], 1),
    'ΔD (°)'    : np.round(BD[latest_idx, :] - BD[0, :], 2),
    'ΔI (°)'    : np.round(BI[latest_idx, :] - BI[0, :], 2)
})

print('\n' + '=' * 80)
print('Magnetic Elements (Z, F, D, I) for Chinese Major Cities (Epoch: 2025-01-01)')
print('=' * 80)
print(df_latest.to_string(index=False))
print('=' * 80)
print('\nSecular Change from 1900 to 2025')
print(df_change.to_string(index=False))
print('=' * 80)

# %%
# # ! 绘图
# Distinct colors and markers for each city
colors = plt.cm.tab10(np.linspace(0, 1, len(CITIES)))
markers = ['o', 's', '^', 'D', 'v', 'p']

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 6),
                        dpi=300, layout='constrained')
# fig.suptitle('Secular Variation of Geomagnetic Elements (Z, F, D, I) (1900–2025)\n'
#              'IGRF-13 Model for Major Chinese Cities', fontsize=15, fontweight='bold')

# Helper function to plot each panel with consistent styling
def plot_panel(ax, data, ylabel, unit):
    for i, city in enumerate(CITIES.keys()):
        ax.plot(DATES, data[:, i], color=colors[i], linestyle='-', linewidth=2.0,
                marker=markers[i], markersize=5, markevery=2, zorder=3, label=city)
    ax.xaxis.set_major_locator(mdates.YearLocator(10))
    ax.xaxis.set_minor_locator(mdates.YearLocator(5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.tick_params(axis='x', rotation=45, labelsize=10)
    ax.minorticks_on()
    ax.grid(True, which='major', linestyle='--', alpha=0.6)
    ax.grid(True, which='minor', linestyle=':', alpha=0.3)
    ax.set_xlabel('Epoch')
    ax.set_ylabel(f'{ylabel} ({unit})')
    ax.legend(loc='lower left', ncol=2, frameon=True, 
              facecolor="#F0ED8E", edgecolor='black', framealpha=0.5)

######## *垂直分量
plot_panel(axs[0, 0], BZ, '$B_Z$', 'nT')

######## *总场强度
plot_panel(axs[0, 1], BF, '$F$', 'nT')

######## *地磁偏角
plot_panel(axs[1, 0], BD, '$D$', '°')

######## *地磁倾角
plot_panel(axs[1, 1], BI, '$I$', '°')

plt.savefig('IGRF13_City.png', dpi=300, bbox_inches='tight')
plt.show()