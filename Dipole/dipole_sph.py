##################################################################################
# # ! Dipole field on the surface of a sphere
##################################################################################
# %%
# # ! Setup
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from magnetic_forward_numba import *
# %%
# # ! Axial Dipole field using IGRF-1990 (Blakely, p. 165)
g1_0 = -29775.0 # IGRF-1990 axial dipole coefficient in nT
a = 6371.2e3 # Earth's mean radius in meters
mu0 = 4 * np.pi * 1e-7 # Permeability of free space in H/m
mx, my = 0., 0. # Dipole moment components in A·m²
mz = 4 * np.pi / mu0 * a**3 * 1.e-9 * g1_0 # Dipole moment in A·m²
dipM = np.array([mx, my, mz])

lon_min, lon_max, n_lon = -180., 180, 37
lat_min, lat_max, n_lat = -90., 90., 37
Lon = np.linspace(lon_min, lon_max, n_lon)
Lat = np.linspace(lat_min, lat_max, n_lat)
LON, LAT = np.meshgrid(Lon, Lat)
obs_lon, obs_lat = LON.ravel(), LAT.ravel()
obs_x, obs_y, obs_z = sph2ecef(obs_lon, obs_lat, a, degrees=True)
obs_xyz = np.column_stack([obs_x, obs_y, obs_z])
obs_bx, obs_by, obs_bz = dipole_field(obs_xyz, dipM, np.array([0, 0, 0]))

# %%
# # ! Plot 3D vector field
obs_mag = np.sqrt(obs_bx**2 + obs_by**2 + obs_bz**2)
norm = plt.Normalize(vmin=obs_mag.min(), vmax=obs_mag.max())
cmap = cm.plasma  # Pick any colormap (e.g., viridis, plasma, coolwarm)
arrow_colors = cmap(norm(obs_mag))
q_colors = np.concatenate([arrow_colors, arrow_colors, arrow_colors])


fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

quiver = ax.quiver(obs_x, obs_y, obs_z, obs_bx, obs_by, obs_bz, 
                   length=0.15 * a, colors=q_colors, normalize=True)

mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
fig.colorbar(mappable, ax=ax, label='nT', shrink=0.6)

ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
ax.set_box_aspect([1, 1, 1])
ax.view_init(elev=0, azim=45)

plt.tight_layout()
plt.show()