##################################################################################
# # ! 1050阶岩石圈磁场模型 -> 区域曲面及平面磁场三分量
# Step 4: 生成地形上数据和平面上数据
##################################################################################
# %% # ! Setup
import numpy as np
import time
import pygmt
import pyshtools as pysh
import xarray as xr
from LW_CordTrans import *
# %%
# # ! ROI & DEM & dt
# ROI
ROI = [106, 107, 33, 34]
lon_min, lon_max, lat_min, lat_max = ROI
# DEM
xr_dem = xr.open_dataarray("west_china_dem_15s_p.nc")
mask_lon = (xr_dem.lon >= lon_min) & (xr_dem.lon <= lon_max)
mask_lat = (xr_dem.lat >= lat_min) & (xr_dem.lat <= lat_max)
xr_dem_roi = xr_dem.where(mask_lon & mask_lat, drop=True)
arr_dem_roi = xr_dem_roi.data
print(f"DEM height diff: {np.nanmax(arr_dem_roi)-np.nanmin(arr_dem_roi)} m")
# dt
xr_mag = xr.open_dataset("west_china_lithomag_01m_g.nc")
xr_dt = xr_mag["dt"]
mask_lon = (xr_dt.lon >= lon_min) & (xr_dt.lon <= lon_max)
mask_lat = (xr_dt.lat >= lat_min) & (xr_dt.lat <= lat_max)
xr_dt_roi = xr_dt.where(mask_lon & mask_lat, drop=True) 
arr_dt_roi = xr_dt_roi.data
print(f"dt diff: {np.nanmax(arr_dt_roi)-np.nanmin(arr_dt_roi)} nT")
# %%
# # ! DEM -> ECEF -> Geocentric
str_calc = "plane" # "curve" or "plane"
res = 0.25/60 # degree
dh_curve = 500     # height above DEM (m)
dh_plane = 200     # height above max_dem (m)
n_lon = int((lon_max-lon_min)/res)+1
n_lat = int((lat_max-lat_min)/res)+1
Lon_gd = np.linspace(lon_min, lon_max, n_lon);
Lat_gd = np.linspace(lat_min, lat_max, n_lat);
LON_gd, LAT_gd = np.meshgrid(Lon_gd, Lat_gd);
xr_dem_interp = xr_dem.interp(lon=Lon_gd, lat=Lat_gd, method="linear")
if str_calc == "curve":
    H_m = xr_dem_interp.data + dh_curve
else:   
    max_dem = np.nanmax(xr_dem_interp.data)
    H_m = np.full(shape=LON_gd.shape, fill_value=max_dem) + dh_plane
H_km = H_m / 1.e3
a, b, f, ecc = GetEllipsoid('wgs84')
ecef_x, ecef_y, ecef_z = ell2ecef(LON_gd, LAT_gd, H_km, 
                                  a=a, ecc=ecc, degrees=True)
LON_gc, LAT_gc, Rs_gc = ecef2sph(ecef_x, ecef_y, ecef_z, degrees=True)
# %%
# # ! Earth’s <main> magnetic field model
calc_lon   = LON_gc.flatten()
calc_lat   = LAT_gc.flatten()
calc_rs_m  = 1.e3 * Rs_gc.flatten()
######## * SHMagCoeffs
clm_mag_main = pysh.SHMagCoeffs.from_file('igrf13coeffs.txt', 
                                           format='igrf', 
                                           r0=6371.2e3,
                                           lmax=13, 
                                           year=2015, 
                                           file_units='nT',
                                           name='IGRF_13 (Earth)', 
                                           units='nT',
                                           encoding='utf-8')
######## * SHMagGrid
t_start = time.perf_counter()
arr_mag_main = clm_mag_main.expand(r=calc_rs_m, 
                                   lat=calc_lat, 
                                   lon=calc_lon,
                                   lmax_calc=clm_mag_main.lmax)
print(f"Time cost for main field expand: {time.perf_counter()-t_start:.3f} s")
Bx_main_gc = -1. * arr_mag_main[:, 1].reshape(LON_gc.shape)
By_main_gc =       arr_mag_main[:, 2].reshape(LON_gc.shape)
Bz_main_gc = -1. * arr_mag_main[:, 0].reshape(LON_gc.shape)
ecef_Bx, ecef_By, ecef_Bz = ned2ecef_v(LON_gc, LAT_gc, 
                                       Bx_main_gc, By_main_gc, Bz_main_gc, 
                                       degrees=True);
Bx_main_gd, By_main_gd, Bz_main_gd = ecef2ned_v(ecef_x, ecef_y, ecef_z, 
                                                ecef_Bx, ecef_By, ecef_Bz, 
                                                a, ecc);
# %%
# # ! Earth’s <lithospheric> magnetic field model
######## * SHMagCoeffs
clm_mag_litho = pysh.SHMagCoeffs.from_file('Thebault_1050.txt', 
                                           r0=6371.2e3,
                                           lmax=1050, 
                                           skip=13,
                                           header=False,  
                                           file_units='nT',
                                           name='Thebault_1050 (Earth)', 
                                           units='nT',
                                           encoding='utf-8')
######## * SHMagGrid
t_start = time.perf_counter()
arr_mag_litho = clm_mag_litho.expand(r=calc_rs_m, 
                                     lat=calc_lat, 
                                     lon=calc_lon,
                                     lmax_calc=clm_mag_litho.lmax)
print(f"Time cost for lithospheric field expand: {time.perf_counter()-t_start:.3f} s")
Bx_litho_gc = -1. * arr_mag_litho[:, 1].reshape(LON_gc.shape)
By_litho_gc =       arr_mag_litho[:, 2].reshape(LON_gc.shape)
Bz_litho_gc = -1. * arr_mag_litho[:, 0].reshape(LON_gc.shape)
ecef_Bx, ecef_By, ecef_Bz = ned2ecef_v(LON_gc, LAT_gc, 
                                       Bx_litho_gc, By_litho_gc, Bz_litho_gc, 
                                       degrees=True);
Bx_litho_gd, By_litho_gd, Bz_litho_gd = ecef2ned_v(ecef_x, ecef_y, ecef_z, 
                                                   ecef_Bx, ecef_By, ecef_Bz, 
                                                   a, ecc);
# %%
# # ! Compute dt & Dt (scalar difference between total and main magnetic field)
Bx_tot = Bx_main_gd + Bx_litho_gd
By_tot = By_main_gd + By_litho_gd
Bz_tot = Bz_main_gd + Bz_litho_gd
B_tot   = np.sqrt(Bx_tot**2 + By_tot**2 + Bz_tot**2)
B0_main = np.sqrt(Bx_main_gd**2 + By_main_gd**2 + Bz_main_gd**2)
# Projection to B0
dt = (Bx_litho_gd*Bx_main_gd + By_litho_gd*By_main_gd + Bz_litho_gd*Bz_main_gd)/B0_main 
# Scalar difference 
Dt = B_tot - B0_main 
# Diff between dt & Dt
Dt_m_dt = Dt - dt
# %% 
# # ! Write data to .txt file
# Planar approximation: center lon/lat as origin, h=0 as origin
R_Earth = 6371.2e3  # Earth's mean radius in meters (consistent with magnetic model reference radius)
mean_lon = (lon_min + lon_max) / 2.0
mean_lat = (lat_min + lat_max) / 2.0
x = R_Earth * np.deg2rad(LAT_gd - mean_lat)
y = R_Earth * np.cos(np.deg2rad(mean_lat)) * np.deg2rad(LON_gd - mean_lon)
z = - H_m 
dx = R_Earth * np.deg2rad(res)
dy = R_Earth * np.cos(np.deg2rad(mean_lat)) * np.deg2rad(res)

with open(f"lithomag_data_{str_calc}.txt", "w") as f:
    # Header comments
    f.write("# ==============================================================================\n")
    f.write("# Lithospheric and Main Magnetic Field Data\n")
    f.write("# ==============================================================================\n")
    f.write(f"# Calculation type: {str_calc}\n")
    f.write(f"# Region: lon [{lon_min}, {lon_max}] deg, lat [{lat_min}, {lat_max}] deg\n")
    f.write(f"# Resolution: {res*60:.2f} arc-minute ({res:.6f} deg)\n")
    if str_calc == "curve":
        f.write(f"# Height: DEM + {dh_curve} m (curve surface)\n")
    else:
        f.write(f"# Height: max(DEM) + {dh_plane} m (plane surface)\n")
    f.write("# Main field model: IGRF-13 (lmax=13, year=2015)\n")
    f.write("# Lithospheric field model: Thebault_1050 (lmax=1050)\n")
    f.write("# Coordinate system: geodetic (WGS84) for lon/lat/H; NED for Bx/By/Bz\n")
    f.write("# Planar approximation (x,y,z): origin at (mean_lon, mean_lat, h=0)\n")
    f.write(f"#   mean_lon = {mean_lon:.6f} deg, mean_lat = {mean_lat:.6f} deg\n")
    f.write(f"#   dx (grid spacing, lat direction) = {dx:.3f} m\n")
    f.write(f"#   dy (grid spacing, lon direction) = {dy:.3f} m\n")
    f.write("# The following planar approx. equations lon/lat in radians\n")
    f.write("#   x = R_Earth * (lat - mean_lat)                 [northward, m]\n")
    f.write("#   y = R_Earth * cos(mean_lat) * (lon - mean_lon) [eastward, m]\n")
    f.write("#   z = -H                                         [downward, m]\n")
    f.write("# Units: lon/lat in degree, H/x/y/z in meter, B in nT\n")
    f.write("# ==============================================================================\n")
    f.write("# Column descriptions:\n")
    f.write("#   1: lon    - geodetic longitude (deg)\n")
    f.write("#   2: lat    - geodetic latitude (deg)\n")
    f.write("#   3: H      - height above ellipsoid (m)\n")
    f.write("#   4: x      - planar x-coordinate (northward) (m)\n")
    f.write("#   5: y      - planar y-coordinate (eastward) (m)\n")
    f.write("#   6: z      - planar z-coordinate (downward) (m)\n")
    f.write("#   7: Bx_main - main field X-component (northward) (nT)\n")
    f.write("#   8: By_main - main field Y-component (eastward) (nT)\n")
    f.write("#   9: Bz_main - main field Z-component (downward) (nT)\n")
    f.write("#  10: Bx_litho - lithospheric field X-component (northward) (nT)\n")
    f.write("#  11: By_litho - lithospheric field Y-component (eastward) (nT)\n")
    f.write("#  12: Bz_litho - lithospheric field Z-component (downward) (nT)\n")
    f.write("#  13: dt     - projection of lithospheric field onto main field (nT)\n")
    f.write("#  14: Dt     - scalar difference |B_total| - |B_main| (nT)\n")
    f.write("#  15: Dt-dt  - difference between scalar and projection approaches (nT)\n")
    f.write("# ==============================================================================\n")
    f.write("# lon(deg)    lat(deg)      H(m)"
            "      x(m)        y(m)       z(m)"
            "   Bx_main   By_main   Bz_main"
            "  Bx_litho  By_litho  Bz_litho"
            "    dt(nT)    Dt(nT)     Dt-dt\n")
    data = np.column_stack([
        LON_gd.ravel(), LAT_gd.ravel(), H_m.ravel(),
        x.ravel(), y.ravel(), z.ravel(),
        Bx_main_gd.ravel(), By_main_gd.ravel(), Bz_main_gd.ravel(),
        Bx_litho_gd.ravel(), By_litho_gd.ravel(), Bz_litho_gd.ravel(),
        dt.ravel(), Dt.ravel(), Dt_m_dt.ravel()
    ])
    np.savetxt(f, data,
               fmt="%10.6f %10.6f %10.3f %10.3f %10.3f %10.3f"
                   " %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f")
    print(f"Data written to lithomag_data_{str_calc}.txt ({data.shape[0]} rows x {data.shape[1]} cols)")