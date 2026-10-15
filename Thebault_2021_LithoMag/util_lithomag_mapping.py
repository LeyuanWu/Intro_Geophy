##################################################################################
# # ! 1050阶岩石圈磁场模型 -> 典型区域磁条带成图
##################################################################################
# %% # ! Setup
import numpy as np
import pygmt
import xarray as xr
# %% # ! Mapping
res = 1./60 # degree
str_reg = "NEPacific" # "NorthAtlantic"
reg_NEPacific = [220., 235., 38., 52.]     # Juan de Fuca Ridge
reg_NorthAtlantic = [325., 340., 57., 63.] # Reykjanes Ridge

if str_reg == "NEPacific":
    region = reg_NEPacific
    lon_min, lon_max, lat_min, lat_max = region
    xr_lithomag = xr.open_dataset(f"{str_reg}_lithomag_{int(res*60):02d}m_g.nc")
elif str_reg == "NorthAtlantic":
    region = reg_NorthAtlantic
    lon_min, lon_max, lat_min, lat_max = region
    xr_lithomag = xr.open_dataset(f"{str_reg}_lithomag_{int(res*60):02d}m_g.nc")

xr_dt = xr_lithomag.dt
fig = pygmt.Figure()
pygmt.makecpt(cmap="polar", series=[-150,150])
fig.grdimage(
    grid=xr_dt,
    region=region,
    projection="M12c",
    frame=["xa5f1g1", "ya5f1g1"])
fig.coast(shorelines="1/0.5p,black")
fig.colorbar(
    frame=[f"x+lMagnetic anomaly", "y+lnT"],
    position="JBC+o0c/0.5c+w10c/0.3c+h")
fig.savefig(f"lithomag_{str_reg}_res{int(res*60):02d}min.png", dpi=300)
fig.show(dpi=300, width=800)