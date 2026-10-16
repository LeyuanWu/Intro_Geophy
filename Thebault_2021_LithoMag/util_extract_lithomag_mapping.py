##################################################################################
# # ! 1050阶岩石圈磁场模型 -> 典型区域磁条带成图
##################################################################################
# %% # ! Setup
import numpy as np
import pygmt
import xarray as xr
# %% # ! Load data
res = 1./60 # degree
str_reg = "EastPacific" # "EastPacific" or "MidAtlantic"
reg_EastPacific = [215, 275, -15., 45.]   
reg_MidAtlantic = [295., 355., -10, 50.] 

if str_reg == "EastPacific":
    region = reg_EastPacific
    lon_min, lon_max, lat_min, lat_max = region
    xr_lithomag = xr.open_dataset(f"{str_reg}_lithomag_{int(res*60):02d}m_g.nc")
elif str_reg == "MidAtlantic":
    region = reg_MidAtlantic
    lon_min, lon_max, lat_min, lat_max = region
    xr_lithomag = xr.open_dataset(f"{str_reg}_lithomag_{int(res*60):02d}m_g.nc")

xr_dt = xr_lithomag.dt

# %% # ! Mapping
cut_mag = 150.
fig = pygmt.Figure()
pygmt.makecpt(cmap="polar", series=[-cut_mag, cut_mag])
fig.grdimage(
    grid=xr_dt,
    region=region,
    projection="M12c",
    frame=["xa10f5g5", "ya10f5g5"])
fig.coast(shorelines="1/0.5p,black")
fig.colorbar(
    frame=[f"x+lMagnetic anomaly", "y+lnT"],
    position="JBC+o0c/1.0c+w10c/0.3c+h")
fig.savefig(f"lithomag_{str_reg}_res{int(res*60):02d}min.png", dpi=300)
fig.show(dpi=300, width=800)