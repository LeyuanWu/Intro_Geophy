##################################################################################
# # ! IGRF-13 地球磁场模型计算得到的中国地区主磁场近百年变化图
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
lat_min = 0;    lat_max = 55;   n_lat = 111;
shc_model = 'IGRF13';
date1 = datetime(1925, 1, 1);
date2 = datetime(2025, 1, 1);
# TODO ********************** end ************************* #
Lon = np.linspace(lon_min, lon_max, n_lon);
Lat = np.linspace(lat_min, lat_max, n_lat);
LON, LAT = np.meshgrid(Lon, Lat);
HGT = h*np.ones(LON.shape);
Bn1, Be1, Bd1 = SHS_gd(LON, LAT, HGT, date1, '%s.shc'%(shc_model));
Bn2, Be2, Bd2 = SHS_gd(LON, LAT, HGT, date2, '%s.shc'%(shc_model));
BX1, BY1, BZ1, BD1, BI1, BH1, BF1 = ComputeMagB7(Bn1, Be1, Bd1);
BX2, BY2, BZ2, BD2, BI2, BH2, BF2 = ComputeMagB7(Bn2, Be2, Bd2);
SV_BX = (BX2-BX1)/1.0; SV_BY = (BY2-BY1)/1.0; SV_BZ = (BZ2-BZ1)/1.0;
SV_BD = (BD2-BD1)/1.0; SV_BI = (BI2-BI1)/1.0; SV_BH = (BH2-BH1)/1.0; SV_BF = (BF2-BF1)/1.0;
######## Save to ascii file
file_MagB7 = shc_model + f'_SV_{date2:%Y-%m-%d}_{date1:%Y-%m-%d}_China.txt';
SaveMagB7File(file_MagB7, LON, LAT, HGT, SV_BX, SV_BY, SV_BZ, SV_BD, SV_BI, SV_BH, SV_BF);
######## The 7 components into a dict
B7s = ['BX','BY','BZ','BD','BI','BH','BF'];
SV_B7 = {}
for (name, bcomp) in zip(['Lon','Lat','H']+B7s,\
        [LON, LAT, HGT, SV_BX, SV_BY, SV_BZ, SV_BD, SV_BI, SV_BH, SV_BF]):
    SV_B7[name] = bcomp;
# %%  
# # ! 地磁要素分布图 (中国地区)
# TODO ******************** Parameters ******************** #
mTitle = f'IGRF-13: Epoch {date2.strftime("%Y")} - Epoch {date1.strftime("%Y")} (China Region)'; 
SUBTTs = ['Change of Declination (D)', 
          'Change of Inclination (I)', 
          'Change of Total Intensity (F)'];
ELEs = ['BD','BI','BF'];
rgn = [lon_min, lon_max, lat_min, lat_max]; 
proj = 'Q12c';
frm = ['WSNE', 'xa10f5g10', 'ya10f5g10'];
mapname = 'China';
sf = 4; # Smooth factor
# TODO ********************** end ************************* #
for subtt, ele in zip(SUBTTs, ELEs):
    data = SV_B7[ele];
    match ele:
        case 'BD':
            lvl = 0.5; anolvl = 1; lbl_plm = 'd3c'; anno = 'f8p+gwhite+o';
        case 'BI':
            lvl = 0.5; anolvl = 1; lbl_plm = 'd3c'; anno = 'f8p+gwhite+o';
        case 'BF':
            lvl = 200; anolvl = 1000; lbl_plm = 'd3c'; anno = 'f8p+gwhite+o';
    # Pens for drawing <Reference & Compare & Difference>
    PENs = {'pos':['athick,red','c0.2p,red'],
            'neg':['athick,blue','c0.2p,blue'],
            'zero':'athick,green'};
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
        fig.savefig(f'IGRF13_SV_{date2:%Y-%m-%d}'
                    f'_{date1:%Y-%m-%d}_{ele}_{mapname}.png', dpi=300);