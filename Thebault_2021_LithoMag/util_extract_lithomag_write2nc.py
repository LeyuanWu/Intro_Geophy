##################################################################################
# # ! 1050阶岩石圈磁场模型 -> 区域曲面及平面磁场三分量
# Step 2: 岩石圈磁场情况 - 全球
# Reference: 
# Thébault, E., Hulot, G., Langlais, B., & Vigneron, P. (2021). 
# A Spherical Harmonic Model of Earth’s Lithospheric Magnetic Field up to Degree 1050. 
# Geophysical Research Letters, 48(21), e2021GL095147. https://doi.org/10.1029/2021GL095147
##################################################################################
# %% # ! Setup
import numpy as np
import time
import pygmt
import pyshtools as pysh
import xarray as xr
# %%
# # ! Earth’s <main> magnetic field model
res = 1./60 # degree
print(f"Grid resolution: {int(60*res)} min")
######## * SHMagCoeffs
clm_mag_main = pysh.SHMagCoeffs.from_file('igrf13coeffs.txt', 
                                           format='igrf', 
                                           r0=6371.2e3,
                                           lmax=13, 
                                           year=2025, 
                                           file_units='nT',
                                           name='IGRF_13 (Earth)', 
                                           units='nT',
                                           encoding='utf-8')
######## * SHMagGrid
t_start = time.perf_counter()
grd_mag_main = clm_mag_main.expand(a=clm_mag_main.r0, 
                                   f=0.0,
                                   lmax=int(90/res)-1,
                                   lmax_calc=clm_mag_main.lmax,
                                   extend=False)
print(f"Time cost for main field expand: {time.perf_counter()-t_start:.3f} s")
Bx_main = -1. * grd_mag_main.theta.data
By_main =       grd_mag_main.phi.data
Bz_main = -1. * grd_mag_main.rad.data
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
grd_mag_litho = clm_mag_litho.expand(a=clm_mag_litho.r0, 
                                     f=0.0,
                                     lmax=int(90/res)-1,
                                     lmax_calc=clm_mag_litho.lmax,
                                     extend=False)
print(f"Time cost for lithospheric field expand: {time.perf_counter()-t_start:.3f} s")
Bx_litho = -1. * grd_mag_litho.theta.data
By_litho =       grd_mag_litho.phi.data
Bz_litho = -1. * grd_mag_litho.rad.data
xr_mag_litho = grd_mag_litho.to_xarray()
# %%
# # ! Compute dt & Dt (scalar difference between total and main magnetic field)
Bx_tot = Bx_main + Bx_litho
By_tot = By_main + By_litho
Bz_tot = Bz_main + Bz_litho
B_tot   = np.sqrt(Bx_tot**2 + By_tot**2 + Bz_tot**2)
B0_main = np.sqrt(Bx_main**2 + By_main**2 + Bz_main**2)
dt = (Bx_litho*Bx_main + By_litho*By_main + Bz_litho*Bz_main)/B0_main # Projection to B0
Dt = B_tot - B0_main # Scalar difference 
# %% # ! Save Bx_litho, By_litho, Bz_litho, dt, Dt to the existing xarray Dataset
# Use xr_mag_litho directly so coordinates are guaranteed aligned
xr_litho_anomaly = xr_mag_litho.copy()

# Assign the corrected components (with proper sign conventions)
xr_litho_anomaly['Bx_litho'] = (('lat', 'lon'), Bx_litho)
xr_litho_anomaly['By_litho'] = (('lat', 'lon'), By_litho)
xr_litho_anomaly['Bz_litho'] = (('lat', 'lon'), Bz_litho)
xr_litho_anomaly['dt'] = (('lat', 'lon'), dt)
xr_litho_anomaly['Dt'] = (('lat', 'lon'), Dt)

# Drop the temporary original variables if desired
xr_litho_anomaly = xr_litho_anomaly.drop_vars(
    ['theta', 'phi', 'radial', 'total', 'potential'], errors='ignore'
)

xr_litho_anomaly.attrs.update({
    'description': 'Lithospheric Magnetic Field Components and Scalar Differences',
    'reference': 'Thébault et al. (2021), GRL',
    'r_calc_m': clm_mag_main.r0,
    'units': 'nT'
})

# %% # ! Extract data belonging to the ROI and save to .nc 
str_reg = "NEPacific"  # or "NorthAtlantic"
reg_NEPacific = [220., 235., 38., 52.]       # Juan de Fuca Ridge
reg_NorthAtlantic = [325., 340., 57., 63.]   # Reykjanes Ridge

if str_reg == "NEPacific":
    lon_min, lon_max, lat_min, lat_max = reg_NEPacific
elif str_reg == "NorthAtlantic":
    lon_min, lon_max, lat_min, lat_max = reg_NorthAtlantic

# Use .sel() — more reliable than .where() for box extraction
# Handle both ascending and descending latitude coordinates automatically
if xr_litho_anomaly.lat[0] > xr_litho_anomaly.lat[-1]:
    # lat is descending (90 -> -90)
    lat_slice = slice(lat_max, lat_min)
else:
    # lat is ascending (-90 -> 90)
    lat_slice = slice(lat_min, lat_max)

xr_roi = xr_litho_anomaly.sel(
    lat=lat_slice,
    lon=slice(lon_min, lon_max)
)

xr_roi.attrs['region'] = f'lon:[{lon_min},{lon_max}], lat:[{lat_min},{lat_max}]'
fname_nc = f"{str_reg}_lithomag_{int(res*60):02d}m_g.nc"
xr_roi.to_netcdf(fname_nc)
print(f"ROI data saved to: {fname_nc}")
print(f"ROI shape: {xr_roi.dims}")
print(f"ROI Bx_litho min/max: {xr_roi.Bx_litho.min().values:.3f} / {xr_roi.Bx_litho.max().values:.3f}")
