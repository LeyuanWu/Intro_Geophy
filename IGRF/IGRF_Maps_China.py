##################################################################################
# # ! IGRF-13 地球磁场模型计算得到的中国地区地磁要素分布图
# # ! Plot components of Geo-Magnetic Model from 'SHC' file for China Region
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
lon_min = 70;   lon_max = 140;  n_lon = 141;
lat_min = 0;   lat_max = 55;   n_lat = 111;
shc_model = 'IGRF13';
date = datetime(1925, 1, 1);
# TODO ********************** end ************************* #
Lon = np.linspace(lon_min, lon_max, n_lon);
Lat = np.linspace(lat_min, lat_max, n_lat);
LON, LAT = np.meshgrid(Lon, Lat);
HGT = h*np.ones(LON.shape);
Bn, Be, Bd = SHS_gd(LON, LAT, HGT, date, '%s.shc'%(shc_model));
BX, BY, BZ, BD, BI, BH, BF = ComputeMagB7(Bn, Be, Bd);
######## Save to ascii file
file_MagB7 = shc_model + f'_{date:%Y-%m-%d}_China.txt';
SaveMagB7File(file_MagB7, LON, LAT, HGT, BX, BY, BZ, BD, BI, BH, BF);
######## The 7 components into a dict
B7s = ['BX','BY','BZ','BD','BI','BH','BF'];
MagB7 = {}
for (name, bcomp) in zip(['Lon','Lat','H']+B7s,\
        [LON, LAT, HGT, BX, BY, BZ, BD, BI, BH, BF]):
    MagB7[name] = bcomp;
# %%  
# # ! 地磁要素分布图 (中国地区)
# TODO ******************** Parameters ******************** #
mTitle = f'IGRF-13 Epoch: {date.strftime("%Y-%m-%d")} (China Region)'; 
SUBTTs = ['Main Field Declination (D)', 
          'Main Field Inclination (I)', 
          'Main Field Total Intensity (F)'];
ELEs = ['BD','BI','BF'];
REGNs = [[lon_min, lon_max, lat_min, lat_max]];
PROJs = ['Q12c'];
FRAMEs = [['WSNE', 'xa10f5g10', 'ya10f5g10']];
MapNames = ['China'];
sf = 4; # Smooth factor
# TODO ********************** end ************************* #
for subtt, ele in zip(SUBTTs, ELEs):
    data = MagB7[ele].reshape(LON.shape);
    match ele:
        case 'BD':
            lvl = 1; anolvl = 2; lbl_plm = 'd3c'; anno = 'f8p+gwhite+o';
        case 'BI':
            lvl = 2; anolvl = 10; lbl_plm = 'd3c'; anno = 'f8p+gwhite+o';
        case 'BX' | 'BY' | 'BZ' | 'BH' | 'BF':
            lvl = 1000; anolvl = 5000; lbl_plm = 'd3c'; anno = 'f8p+gwhite+o';
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
            fig.coast(water='lightblue', land='white', transparency=25,
                      rivers='1/0.5p,lightblue',
                      shorelines=['1/0.5p,black', '2/0.3p,black']);
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
