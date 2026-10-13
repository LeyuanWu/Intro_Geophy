# # Plot DEM
# fig_dem = pygmt.Figure()
# fig_dem.grdimage(
#     grid=xr_dem_sub,
#     region=region,
#     projection="Q12c",
#     frame="a0.5f0.5g0.5",
#     cmap='geo'
# )
# fig_dem.colorbar(
#     frame=["xa200+lElevation", "y+lm"],
#     position="JBC+o0c/0.8c+w12c/0.5c+h"
# )
# fig_dem.savefig("roi_dem.png", dpi=300)
# fig_dem.show(dpi=300, width=800)
# # Plot dt
# fig_dt = pygmt.Figure()
# fig_dt.grdimage(
#     grid=xr_dt_sub,
#     region=region,
#     projection="Q12c",
#     frame="a0.5f0.5g0.5",
#     cmap='polar'
# )
# fig_dt.colorbar(
#     frame=["x+lLithospheric dt", "y+lnT"],
#     position="JBC+o0c/0.8c+w12c/0.5c+h"
# )
# fig_dt.savefig("roi_dt.png", dpi=300)
# fig_dt.show(dpi=300, width=800)