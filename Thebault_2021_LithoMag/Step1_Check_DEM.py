##################################################################################
# # ! 1050阶岩石圈磁场模型 -> 区域曲面及平面磁场三分量
# Step 1: 分析地形情况
##################################################################################
# %% # ! Setup
import numpy as np
import pygmt
import xarray as xr
# %%
# # ! 下载中国西部区域 15s 分辨率 DEM 数据并保存为 NetCDF 文件
region = [70, 110, 15, 55]
# dem = pygmt.datasets.load_earth_relief(
#     resolution="15s",
#     region=region,
#     registration="pixel"
# )
# dem.to_netcdf("west_china_dem_15s_p.nc")
# print("=" * 40)
# print(f"DEM saved to: west_china_dem_15s_p.nc")
# print(f"Region: {region}")
# print(f"Grid shape: {dem.shape}")
# print(f"Resolution: 15 arc-seconds (~450 m)")
# print("=" * 40)
# %%
# # ! 读取保存的 DEM 文件并使用 pygmt 绘图
xr_dem = xr.open_dataarray("west_china_dem_15s_p.nc")
arr_dem = xr_dem.data
h_max = arr_dem.max()
h_min = arr_dem.min()
h_cmin = np.round(h_min, decimals=-2)
h_cmax = np.round(h_max, decimals=-2)

fig_dem = pygmt.Figure()
fig_dem.grdimage(
    grid=xr_dem,
    region=region,
    projection="Q12c",
    frame="a5f5g1",
    cmap="geo"
)
fig_dem.coast(shorelines="1/0.5p,black", borders="1/0.5p,black")
fig_dem.colorbar(
    frame=["x+lElevation", "y+lm"],
    position="JBC+o0c/0.8c+w12c/0.5c+h"
)
# %%
# # ! 计算每个 1°×1° 地形块的最大-最小高度差值，并将数值绘制在地图上
lons = xr_dem.lon.data  # 1D 经度坐标
lats = xr_dem.lat.data  # 1D 纬度坐标
lon_edges = np.arange(int(np.floor(lons.min())), int(np.ceil(lons.max())))
lat_edges = np.arange(int(np.floor(lats.min())), int(np.ceil(lats.max())))
n_lon = len(lon_edges)
n_lat = len(lat_edges)

# 存储每个 1°×1° patch 的最大-最小高差 (km)
dh = np.full((n_lat, n_lon), np.nan)
for i, lon0 in enumerate(lon_edges):
    ix = np.where((lons >= lon0) & (lons <= lon0 + 1))[0]
    for j, lat0 in enumerate(lat_edges):
        jx = np.where((lats >= lat0) & (lats <= lat0 + 1))[0]
        patch = arr_dem[np.ix_(jx, ix)]
        dh[j, i] = (patch.max() - patch.min()) / 1000.0

lon_c = lon_edges + 0.5
lat_c = lat_edges + 0.5

for j in range(n_lat):
    for i in range(n_lon):
        if np.isnan(dh[j, i]):
            continue
        if dh[j, i] > 1.0:
            font = "5p,Helvetica-Bold,red"
        else:
            font = "4p,Helvetica,black"
        fig_dem.text(
            x=lon_c[i],
            y=lat_c[j],
            text=f"{dh[j, i]:.1f}",
            font=font,
            justify="CM"
        )
fig_dem.savefig("west_china_dem_15s_p.png", dpi=300)
fig_dem.show(dpi=300, width=800)
