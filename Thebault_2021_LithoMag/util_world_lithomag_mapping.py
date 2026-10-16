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
cut_mag = 150.
fig = pygmt.Figure()
pygmt.makecpt(cmap="polar", series=[-cut_mag, cut_mag])
fig.grdimage(
    grid=xr_dt,
    region=[-139, -123, 46, 55],
    projection="L-132/50/46/54/12c",
    frame=["xa3f3g3", "ya2f2g2"])
fig.coast(shorelines="1/0.5p,black")
fig.colorbar(
    frame=[f"x+lMagnetic anomaly", "y+lnT"],
    position="JBC+o0c/1.0c+w10c/0.3c+h")
fig.show(dpi=300, width=800)