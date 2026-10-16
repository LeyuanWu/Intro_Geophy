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
# %%
# # ! Save Bx_litho, By_litho, Bz_litho, dt, Dt to a single xarray.DataArray
lat = xr_mag_litho.lat
lon = xr_mag_litho.lon
data_vars = {
    'dt': (('lat', 'lon'), dt)
}   
xr_litho_anomaly = xr.Dataset(
    data_vars=data_vars,
    coords={'lat': lat,'lon': lon},
    attrs={'description': 'Lithospheric Magnetic dt',
           'reference': 'Thébault et al. (2021), GRL',
           'r_calc_m': clm_mag_main.r0,
           'units': 'nT'})

fname_nc = f"world_lithomag_{int(res*60):02d}m_g.nc"
xr_litho_anomaly.to_netcdf(fname_nc)
print(f"Data saved to: {fname_nc}")