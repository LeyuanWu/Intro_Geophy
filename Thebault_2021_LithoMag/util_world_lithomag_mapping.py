##################################################################################
# # ! 1050阶岩石圈磁场模型 -> 全球成图
##################################################################################
# %% # ! Setup
import numpy as np
import pygmt
import xarray as xr
# %% # ! Load data
res = 1./60 # degree
xr_lithomag = xr.open_dataset(f"world_lithomag_{int(res*60):02d}m_g.nc")
xr_dt = xr_lithomag.dt
# %% # ! Mapping
str_rgn= "northpole" # "northpole" or "southpole"
cut_mag = 150.
if str_rgn == "northpole":
    xr_dt = xr_dt.sel(lat=slice(90, 0))
    rgn = [0, 360, 60, 90]
    proj = 'S0/90/30/10c';
elif str_rgn == "southpole":
    xr_dt = xr_dt.sel(lat=slice(0, -90))
    rgn = [0, 360, -90, -60]
    proj = 'S0/-90/30/10c';

fig = pygmt.Figure()
pygmt.makecpt(cmap="polar", series=[-cut_mag, cut_mag])
fig.grdimage(
    grid=xr_dt,
    region=rgn,
    projection=proj,
    frame=["xa30f30g30", "ya15f15g15"])
fig.coast(shorelines="1/0.5p,black")
fig.colorbar(
    frame=[f"x+lMagnetic anomaly", "y+lnT"],
    position="JBC+o0c/1.0c+w10c/0.3c+h")
fig.savefig(f"lithomag_{str_rgn}_res{int(res*60):02d}min.png", dpi=300)
fig.show(dpi=300, width=800)

