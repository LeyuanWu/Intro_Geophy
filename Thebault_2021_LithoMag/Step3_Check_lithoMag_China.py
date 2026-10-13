##################################################################################
# # ! 1050阶岩石圈磁场模型 -> 区域曲面及平面磁场三分量
# Step 3: 分析岩石圈磁场标量差 (dt) 分布情况 - 中国西部
##################################################################################
# %% # ! Setup
import numpy as np
import pygmt
import xarray as xr
# %%
# # ! 读取岩石圈磁场数据文件并使用 pygmt 绘制 dt 分布图
region = [70, 110, 15, 55]
xr_mag = xr.open_dataset("west_china_mag_01m_g.nc")
xr_dt = xr_mag["dt"]
arr_dt = xr_dt.data

fig_dt = pygmt.Figure()
pygmt.makecpt(cmap="polar", series=[-300, 300])
fig_dt.grdimage(
    grid=xr_dt,
    region=region,
    projection="Q12c",
    frame="a5f5g1"
)
fig_dt.coast(shorelines="1/0.5p,black", borders="1/0.5p,black")
fig_dt.colorbar(
    frame=["x+lLithospheric dt", "y+lnT"],
    position="JBC+o0c/0.8c+w12c/0.5c+h"
)
# %%
# # ! 计算每个 1°×1° 块的最大-最小 dt 差值，并将数值绘制在地图上
lons = xr_dt.lon.data  # 1D 经度坐标
lats = xr_dt.lat.data  # 1D 纬度坐标
lon_edges = np.arange(int(np.floor(lons.min())), int(np.ceil(lons.max())))
lat_edges = np.arange(int(np.floor(lats.min())), int(np.ceil(lats.max())))
n_lon = len(lon_edges)
n_lat = len(lat_edges)

# 存储每个 1°×1° patch 的最大-最小 dt 差值 (nT)
ddif = np.full((n_lat, n_lon), np.nan)
for i, lon0 in enumerate(lon_edges):
    ix = np.where((lons >= lon0) & (lons <= lon0 + 1))[0]
    for j, lat0 in enumerate(lat_edges):
        jx = np.where((lats >= lat0) & (lats <= lat0 + 1))[0]
        patch = arr_dt[np.ix_(jx, ix)]
        ddif[j, i] = patch.max() - patch.min()

lon_c = lon_edges + 0.5
lat_c = lat_edges + 0.5

for j in range(n_lat):
    for i in range(n_lon):
        if np.isnan(ddif[j, i]):
            continue
        if ddif[j, i] > 300.0:
            font = "3p,Helvetica-Bold,black"
        else:
            font = "3p,Helvetica,black"
        fig_dt.text(
            x=lon_c[i],
            y=lat_c[j],
            text=f"{ddif[j, i]:.0f}",
            font=font,
            justify="CM"
        )
fig_dt.savefig("west_china_lithoMag_dt.png", dpi=300)
fig_dt.show(dpi=300, width=800)

