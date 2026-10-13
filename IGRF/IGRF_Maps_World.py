##################################################################################
# # ! IGRF-13 地球磁场模型计算得到的地磁要素分布图
# # ! Plot components of Geo-Magnetic Model from 'SHC' file
##################################################################################
# %% # ! Setup
import pygmt;
import xarray as xr;
import numpy as np;
from datetime import datetime;
from LW_MAG import *;
from LW_Stat import *;
# %%
# # ! Estimations of MagB7 at <h> using the code 
# # ! https://github.com/IAGA-VMOD/ppigrf.git
# TODO ******************** Parameters ******************** #
h = 0.0;
lon_min = -180;  lon_max = 180;  n_lon = 361;
lat_min = -89.5; lat_max = 89.5; n_lat = 180;
shc_model = 'IGRF13';
date = datetime(2025, 1, 1);
# TODO ********************** end ************************* #
Lon = np.linspace(lon_min, lon_max, n_lon);
Lat = np.linspace(lat_min, lat_max, n_lat);
LON, LAT = np.meshgrid(Lon, Lat);
HGT = h*np.ones(LON.shape);
Bn, Be, Bd = SHS_gd(LON, LAT, HGT, date, '%s.shc'%(shc_model));
BX, BY, BZ, BD, BI, BH, BF = ComputeMagB7(Bn, Be, Bd);
######## Save to ascii file
file_MagB7 = shc_model + f'_{date:%Y-%m-%d}.txt';
SaveMagB7File(file_MagB7, LON, LAT, HGT, BX, BY, BZ, BD, BI, BH, BF);
######## The 7 components into a dict
B7s = ['BX','BY','BZ','BD','BI','BH','BF'];
MagB7 = {}
for (name, bcomp) in zip(['Lon','Lat','H']+B7s,\
        [LON, LAT, HGT, BX, BY, BZ, BD, BI, BH, BF]):
    MagB7[name] = bcomp;
# %%  
# # ! 地磁要素分布图
# TODO ******************** Parameters ******************** #
mTitle = f'IGRF-13 Epoch: {date.strftime("%Y-%m-%d")}'; 
SUBTTs = ['Main Field Declination (D)', 
          'Main Field Inclination (I)', 
          'Main Field Total Intensity (F)'];
ELEs = ['BD','BI','BF'];
REGNs = ['d'];
PROJs = ['Q18c'];
FRAMEs = [['WSNE', 'xa30f30g30', 'ya15f15g15']];
MapNames = ['Cylinder'];
sf = 4; # Smooth factor
# TODO ********************** end ************************* #
for subtt, ele in zip(SUBTTs, ELEs):
    data = MagB7[ele].reshape(LON.shape);
    match ele:
        case 'BD':
            lvl = 5; anolvl = 10; lbl_plm = 'd5c'; anno = 'f6p+gwhite+o';
        case 'BI':
            lvl = 5; anolvl = 20; lbl_plm = 'd5c'; anno = 'f6p+gwhite+o';
        case 'BX' | 'BY' | 'BZ' | 'BH' | 'BF':
            lvl = 2500; anolvl = 5000; lbl_plm = 'd5c'; anno = 'f6p+gwhite+o';
    # Pens for drawing <Reference & Compare & Difference>
    PENs = {'pos':['athin,red','c0.05p,red'],
            'neg':['athin,blue','c0.05p,blue'],
            'zero':'athin,green'};
    for rgn, proj, frm, mapname in zip(REGNs, PROJs, FRAMEs, MapNames):
        fig = pygmt.Figure();
        with pygmt.config(FONT_ANNOT_PRIMARY='9p,Helvetica',
                          FONT_LABEL='9p,Helvetica-Bold',
                          MAP_FRAME_TYPE='fancy', MAP_FRAME_PEN='thin,black'):
            fig.basemap(region=rgn, projection=proj, frame=frm);
            fig.coast(water='lightblue', land='white', transparency=25);
            grdLat = Lat; grdLon = Lon; grdB = data;
            grid = xr.DataArray(grdB, coords=[('latitude', grdLat), ('longitude', grdLon)]);
            grid.gmt.registration = 0; grid.gmt.gtype = 1;
            if grid.max()>0:
                fig.grdcontour(grid=grid, levels=lvl, annotation=f'{anolvl}+{anno}', 
                                limit='p', label_placement=lbl_plm, resample=sf,
                                pen=PENs['pos']);
            if grid.min()<0:
                fig.grdcontour(grid=grid, levels=lvl, annotation=f'{anolvl}+{anno}', 
                                limit='n', label_placement=lbl_plm, resample=sf,
                                pen=PENs['neg']);
            if grid.min()<0 and grid.max()>0:
                fig.grdcontour(grid=grid, levels=lvl, annotation=f'{anolvl}+{anno}', 
                                limit=[-1e-6,1e-6], label_placement=lbl_plm, resample=sf,
                                pen=PENs['zero']);  
            fig.show(dpi=400, width=800);
            fig.savefig(f'IGRF13_{date.strftime("%Y-%m-%d")}_{ele}_{mapname}.png', dpi=300);