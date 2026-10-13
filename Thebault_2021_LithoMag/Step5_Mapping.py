##################################################################################
# Step 5: 读取 curve / plane 磁场数据并绘制 2×2 子图
#   a) DEM
#   b) dt on curve surface
#   c) dt on plane surface
#   d) dt curve - dt plane
##################################################################################
# %% # ! Setup
import numpy as np
import pygmt
import xarray as xr

# %% # ! ROI
ROI = [106, 107, 33, 34]
lon_min, lon_max, lat_min, lat_max = ROI

# %% # ! Load DEM
xr_dem = xr.open_dataarray("west_china_dem_15s_p.nc")
mask_lon = (xr_dem.lon >= lon_min) & (xr_dem.lon <= lon_max)
mask_lat = (xr_dem.lat >= lat_min) & (xr_dem.lat <= lat_max)
xr_dem_roi = xr_dem.where(mask_lon & mask_lat, drop=True)

# %% # ! Helper: load txt -> xarray DataArray
def load_dt_to_xarray(filepath):
    """
    读取 Step4 输出的 txt 文件，将 lon/lat/dt 转换为 xarray.DataArray。
    坐标直接按 Step4 的 meshgrid 参数生成，避免文本回读带来的浮点精度问题。
    """
    data = np.loadtxt(filepath, comments="#")
    dt = data[:, 12]   # dt 为第 13 列 (0-based index 12)

    res = 0.25 / 60
    n_lon = int((lon_max - lon_min) / res) + 1
    n_lat = int((lat_max - lat_min) / res) + 1
    lons = np.linspace(lon_min, lon_max, n_lon)
    lats = np.linspace(lat_min, lat_max, n_lat)

    # Step4 中按 LON_gd.ravel(), LAT_gd.ravel() 写入，shape 为 (n_lat, n_lon)
    dt_grid = dt.reshape(n_lat, n_lon)

    da = xr.DataArray(
        dt_grid,
        coords=[("lat", lats), ("lon", lons)],
        name="dt",
    )
    return da

# %% # ! Load curve & plane dt
xr_dt_curve = load_dt_to_xarray("lithoMag_data_curve.txt")
xr_dt_plane = load_dt_to_xarray("lithoMag_data_plane.txt")

# 计算差值
dt_diff = xr_dt_curve - xr_dt_plane
dt_diff.name = "dt_diff"

print(f"dt_curve  range: {float(xr_dt_curve.min()):.2f} ~ {float(xr_dt_curve.max()):.2f} nT")
print(f"dt_plane  range: {float(xr_dt_plane.min()):.2f} ~ {float(xr_dt_plane.max()):.2f} nT")
print(f"dt_diff   range: {float(dt_diff.min()):.2f} ~ {float(dt_diff.max()):.2f} nT")


# %% # ! Helper: round up to nearest multiple
def round_up_to_multiple(value, multiple):
    return int(np.ceil(value / multiple)) * multiple


# %% # ! Make symmetric CPTs for polar colormap
# b/c 共享同一色标，保证两者可直接对比
dt_min = min(float(xr_dt_curve.min()), float(xr_dt_plane.min()))
dt_max = max(float(xr_dt_curve.max()), float(xr_dt_plane.max()))
abs_max_dt = max(abs(dt_min), abs(dt_max))
rounded_max_dt = round_up_to_multiple(abs_max_dt, 20)
pygmt.makecpt(cmap="polar", series=[-rounded_max_dt, rounded_max_dt], output="dt_shared.cpt")

# d) 差值单独配色，仍令 0 位于色标中心
abs_max_diff = max(abs(float(dt_diff.min())), abs(float(dt_diff.max())))
rounded_max_diff = round_up_to_multiple(abs_max_diff, 5)
pygmt.makecpt(cmap="polar", series=[-rounded_max_diff, rounded_max_diff], output="dt_diff.cpt")

# %% # ! Plot 2×2 subplots
fig = pygmt.Figure()

with fig.subplot(
    nrows=2,
    ncols=2,
    figsize=("16c", "18c"),
    frame="a1f0.25g0.25",
    autolabel="a)+o0.1c/0.3c",
    margins=["1.2c","0.4c"],
):
    # ---------- a) DEM ----------
    with fig.set_panel(panel=0):
        fig.grdimage(
            grid=xr_dem_roi,
            region=ROI,
            projection="M?",
            cmap="geo",
        )
        fig.colorbar(
            frame=["xa500+lElevation", "y+lm"],
            position="JMR+o0.4c/0c+w7c/0.3c",
        )

    # ---------- b) dt on curve ----------
    with fig.set_panel(panel=1):
        fig.grdimage(
            grid=xr_dt_curve,
            region=ROI,
            projection="M?",
            cmap="dt_shared.cpt",
        )
        fig.colorbar(
            frame=["xa20+lLithospheric dt (curve)", "y+lnT"],
            position="JMR+o0.4c/0c+w7c/0.3c",
        )

    # ---------- c) dt on plane ----------
    with fig.set_panel(panel=2):
        fig.grdimage(
            grid=xr_dt_plane,
            region=ROI,
            projection="M?",
            cmap="dt_shared.cpt",
        )
        fig.colorbar(
            frame=["xa20+lLithospheric dt (plane)", "y+lnT"],
            position="JMR+o0.4c/0c+w7c/0.3c",
        )

    # ---------- d) dt curve - dt plane ----------
    with fig.set_panel(panel=3):
        fig.grdimage(
            grid=dt_diff,
            region=ROI,
            projection="M?",
            cmap="dt_diff.cpt",
        )
        fig.colorbar(
            frame=["xa5+lCurve - Plane", "y+lnT"],
            position="JMR+o0.4c/0c+w7c/0.3c",
        )

fig.savefig(f"lithoMag_mapping.png", dpi=300)
fig.show(dpi=300, width=800)
print("Figure saved to lithoMag_mapping.png")
