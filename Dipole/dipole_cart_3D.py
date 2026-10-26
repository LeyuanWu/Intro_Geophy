##################################################################################
# # ! Magnetic field on a plane above a dipole
##################################################################################
# %%
# # ! Setup
import numpy as np
import pyvista as pv
from magnetic_forward_numba import *
# %%
# # ! Vertical/Horizontal Dipole/Component
mu0 = 4 * np.pi * 1e-7   # Permeability of free space in H/m
M = 100.                 # Dipole moment magnitude in A·m²
dip_orient = 'D'         # 'D' = Downward dipole,    'N' = Northward dipole
comp = 'N'               # 'D' = Downward component, 'N' = Northward component
if dip_orient.upper() == 'D':
    mx, my, mz = 0., 0., M
elif dip_orient.upper() == 'N':
    mx, my, mz = M, 0., 0.
dip_M = np.array([mx, my, mz])
pos_M = np.array([0, 0, 10.])
comp_info = {
    'D': ('bd', r'$B_z$ (nT)'),
    'N': ('bn', r'$B_x$ (nT)')
}
comp_name, comp_label = comp_info[comp.upper()]
# %%
# # ! Corse grid for arrows
nor_min, nor_max, n_nor = -20., 20., 9
est_min, est_max, n_est = -20., 20., 9
dwd_min, dwd_max, n_dwd = -10., 0., 3
Nor = np.linspace(nor_min, nor_max, n_nor)
Est = np.linspace(est_min, est_max, n_est)
Dwd = np.linspace(dwd_min, dwd_max, n_dwd)
NOR, EST, DWD = np.meshgrid(Nor, Est, Dwd)
obs_nor, obs_est, obs_dwd = NOR.ravel(), EST.ravel(), DWD.ravel()
obs_ned = np.column_stack((obs_nor, obs_est, obs_dwd))
obs_bn, obs_be, obs_bd = dipole_field(obs_ned, dip_M, pos_M)
mag_B = np.sqrt(obs_bn**2 + obs_be**2 + obs_bd**2)
obs_B = np.column_stack((obs_bn, obs_be, obs_bd))
# %%
# # ! Fine grid for contour
nor_1, nor_2, Nnor = -20., 20., 81
est_1, est_2, Nest = -20., 20., 81
NorF = np.linspace(nor_1, nor_2, Nnor)
EstF = np.linspace(est_1, est_2, Nest)
NORF, ESTF = np.meshgrid(NorF, EstF)
DWDF = 0.0 * np.ones(NORF.shape)
obs_nor_f, obs_est_f, obs_dwd_f = NORF.ravel('F'), ESTF.ravel('F'), DWDF.ravel('F')
obs_ned_f = np.column_stack((obs_nor_f, obs_est_f, obs_dwd_f))
obs_bn_f, obs_be_f, obs_bd_f = dipole_field(obs_ned_f, dip_M, pos_M)
# Select the scalar array for contour / surface plot based on comp
obs_bc_f = obs_bn_f if comp_name == 'bn' else obs_bd_f
GridZ0 = pv.StructuredGrid(NORF, ESTF, DWDF)
GridZ0.point_data['bn'] = obs_bn_f
GridZ0.point_data['be'] = obs_be_f
GridZ0.point_data['bd'] = obs_bd_f
# %%
# # ! Dipole marker
dir_M = dip_M / np.linalg.norm(dip_M)
arrow_scale = 5.0
dipole_arrow = pv.Arrow(
    start=pos_M - 0.5 * arrow_scale * dir_M, direction=dir_M,
    tip_length=0.25, tip_radius=0.12, shaft_radius=0.04,
    tip_resolution=24, shaft_resolution=24, scale=arrow_scale,
)
# %%
# # ! <Pyvista> 3D plot
mesh_vec = pv.PolyData(obs_ned)
mesh_vec.point_data['vectors'] = obs_B
mesh_vec.point_data['magnitude'] = np.log(2.0 + mag_B)

mag_arrow = pv.Arrow(
    tip_length=0.30,      # Length of the cone tip (fraction of total)
    tip_radius=0.10,      # Thickness of the cone tip
    shaft_radius=0.03,    # Thickness of the cylinder shaft
    tip_resolution=24,    # Smoothness of the tip
    shaft_resolution=24   # Smoothness of the shaft
)

mesh_magvec = mesh_vec.glyph(
    orient='vectors', scale='magnitude', factor=1.0, 
    geom=mag_arrow, progress_bar=False)

pl = pv.Plotter(notebook=True, window_size=[800, 600], image_scale=3)
pl.set_background('white') 

# * Magnetic field arrows
pl.add_mesh(mesh_magvec, color='black', show_scalar_bar=False)

# * Contour
bc_max = np.abs(obs_bc_f).max()
contours = GridZ0.contour(isosurfaces=20, scalars=comp_name, rng=[-bc_max, bc_max])
pl.add_mesh(contours, color='black', line_width=4, render_lines_as_tubes=True)

pl.add_mesh(
    GridZ0, scalars=comp_name, cmap='seismic', clim=[-bc_max, bc_max],
    opacity=0.8, show_scalar_bar=True,
    scalar_bar_args=dict(title=comp_label,
                         title_font_size=16, 
                         label_font_size=14,
                         position_x=0.3, 
                         position_y=0.1, 
                         n_labels=5, fmt='%.1f'))

# Dipole marker
pl.add_mesh(dipole_arrow, color='cyan', show_scalar_bar=False)

# Case label
pl.add_text(f"{dip_orient.upper()}-dipole / {comp.upper()}-component",
            position='upper_left', font_size=28, color='black',
            font='courier', shadow=True)


# View  
if dip_orient == 'D'and comp == 'D':
    az, el = 160, -18
elif dip_orient == 'D'and comp == 'N':
    az, el = 90, -18
elif dip_orient == 'N'and comp == 'D':
    az, el = 90, -18
elif dip_orient == 'N'and comp == 'N':
    az, el = 90, -18
v_az, v_el = np.deg2rad(az), np.deg2rad(el)
pl.view_vector([np.cos(v_el)*np.cos(v_az),
                np.cos(v_el)*np.sin(v_az),
                np.sin(v_el)])
pl.camera.up = (0, 0, -1)   
pl.camera.zoom(1.6)
pl.camera.focal_point = (0., -3., 0.)

######## * Choose for interactive 
# pl.add_axes(color='black')
# pl.show(jupyter_backend='client')

######## * Choose for static png
pl.add_axes(color='black', line_width=6, label_size=(0.04, 0.04),
            viewport=(0.0, 0.0, 0.3, 0.3))
pl.show(jupyter_backend='static')
