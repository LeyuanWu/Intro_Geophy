##################################################################
# # ! Leyuan Wu
# Magnetic Field Modelling functions
##################################################################

import numpy as np
import pandas as pd
from LW_CordTrans import *

def is_leapyear(year):
    """ Check for leapyear (handles arrays and preserves shape)

    """

    # if array:
    if type(year) is np.ndarray:
        out = np.full_like(year, False, dtype = bool)

        out[ year % 4   == 0] = True
        out[ year % 100 == 0] = False
        out[ year % 400 == 0] = True

        return out

    # if scalar:
    if year % 400 == 0:
        return True

    if year % 100 == 0:
        return False

    if year % 4 == 0:
        return True

    else:
        return False

def yearfrac_to_datetime(fracyear):
    """ 
    Convert fraction of year to datetime 

    Parameters
    ----------
    fracyear : iterable
        Date(s) in decimal year. E.g., 2021-03-28 is 2021.2377
        Must be an array, list or similar.

    Returns
    -------
    datetimes : array
        Array of datetimes
    """

    year = np.uint16(fracyear) # truncate fracyear to get year
    # use pandas to_timedelta to represent time since beginning of year: 
    delta_year = pd.to_timedelta((fracyear - year)*(365 + is_leapyear(year)), unit = 'D')
    # and DatetimeIndex to represent beginning of years:
    start_year = pd.DatetimeIndex(list(map(str, year)))
 
    # adding them produces the datetime:
    return (start_year + delta_year).to_pydatetime()

def get_legendre(theta, keys):
    """ 
    Calculate Schmidt semi-normalized associated Legendre functions

    Calculations based on recursive algorithm found in "Spacecraft Attitude Determination and Control" by James Richard Wertz
    
    Parameters
    ----------
    theta : array
        Array of colatitudes in degrees
    keys: iterable
        list of spherical harmnoic degree and order, tuple (n, m) for each 
        term in the expansion

    Returns
    -------
    P : array 
        Array of Legendre functions, with shape (theta.size, len(keys)). 
    dP : array
        Array of dP/dtheta, with shape (theta.size, len(keys))
    """

    # get maximum N and maximum M:
    n, m = np.array([k for k in keys]).T
    nmax, mmax = np.max(n), np.max(m)

    theta = theta.flatten()[:, np.newaxis]
    theta_rad = np.radians(theta)

    P = {}
    dP = {}
    sinth = np.sin(theta_rad)
    costh = np.cos(theta_rad)

    # Initialize Schmidt normalization
    S = {}
    S[0, 0] = 1.

    # initialize the functions:
    for n in range(nmax +1):
        for m in range(nmax + 1):
            P[n, m] = np.zeros_like(theta, dtype = np.float64)
            dP[n, m] = np.zeros_like(theta, dtype = np.float64)

    P[0, 0] = np.ones_like(theta, dtype = np.float64)
    for n in range(1, nmax +1):
        for m in range(0, min([n + 1, mmax + 1])):
            # do the legendre polynomials and derivatives
            if n == m:
                P[n, n]  = sinth * P[n - 1, m - 1]
                dP[n, n] = sinth * dP[n - 1, m - 1] + costh * P[n - 1, n - 1]
            else:

                if n == 1:
                    Knm = 0.
                    P[n, m]  = costh * P[n -1, m]
                    dP[n, m] = costh * dP[n - 1, m] - sinth * P[n - 1, m]

                elif n > 1:
                    Knm = ((n - 1)**2 - m**2) / ((2*n - 1)*(2*n - 3))
                    P[n, m]  = costh * P[n -1, m] - Knm*P[n - 2, m]
                    dP[n, m] = costh * dP[n - 1, m] - sinth * P[n - 1, m] - Knm * dP[n - 2, m]

            # compute Schmidt normalization
            if m == 0:
                S[n, 0] = S[n - 1, 0] * (2.*n - 1)/n
            else:
                S[n, m] = S[n, m - 1] * np.sqrt((n - m + 1)*(int(m == 1) + 1.)/(n + m))


    # now apply Schmidt normalization
    for n in range(1, nmax + 1):
        for m in range(0, min([n + 1, mmax + 1])):
            P[n, m]  *= S[n, m]
            dP[n, m] *= S[n, m]

    Pmat  = np.hstack(tuple(P[key] for key in keys))
    dPmat = np.hstack(tuple(dP[key] for key in keys)) 
    return Pmat, dPmat    

def nm2ind(n, m):
    """Calculate the 1d index of ALFs with degree n and order m

    Args:
        n (int)     : degree
        m (int)     : order
        
    Return:
        inm (int)   : the 1d index
    """
    inm = int(n*(n+1)/2 + m);
    return inm;

def PnmSchmidt_SFCR(nmax, theta, flag_dP = False, flag_mP = False):
    """Calculate Schmidt semi-normalized associated Legendre functions (SNALFs)
    (optinal) First-Order Derivative dPnm/dtheta & (optinal) First Kind Spherical Harmonic Polynomial m*Pnm/sin(theta)
    Calculations based on standard forward column recursive algorithm found in 
    - Liu, X., Xu, T., Sun, B., Li, X., Wang, Y., & Guan, B. (2021). Non-singular calculation of geomagnetic vectors and geomagnetic gradient tensors. Journal of Geophysical Research: Solid Earth, 126, e2021JB023397. https://doi.org/10.1029/2021JB023397
    
    Args:
        nmax (int)              : maximum degree
        theta (numpy array)     : colatitudes in degrees
        flag_dP (bool, optinal) : flag to compute First-Order Derivative dPnm/dtheta, default False
        flag_mP (bool, optinal) : flag to compute First Kind Spherical Harmonic Polynomial m*Pnm/sin(theta), default False
        
    Returns:
        Pnm (numpy array)   : Array of SNALFs, with shape (theta.size, (nmax+1)*(nmax+2)/2)
        dPnm (numpy array)  : Array of dP/dtheta, with shape (theta.size, len(keys))
    """
    N = int((nmax+1)*(nmax+2)/2);
    theta = theta.flatten();
    theta_rad = np.deg2rad(theta);
    sinth = np.sin(theta_rad);
    costh = np.cos(theta_rad);
    Pnm = np.zeros((theta.size, N));
    Pnm[:,0] = np.ones_like(theta);
    for m in range(1, nmax+1):
        inm1 = nm2ind(m-1, m-1);
        inm2 = nm2ind(m, m);
        if m==1:
            sec_nm = np.sqrt((2*m-1)/m);
        else:
            sec_nm = np.sqrt((2*m-1)/(2*m));
        Pnm[:,inm2] = sec_nm * sinth * Pnm[:,inm1];
    for m in range(0, nmax):
        inm1 = nm2ind(m, m);
        inm2 = nm2ind(m+1, m);
        Pnm[:,inm2] = np.sqrt(2*m+1) * costh * Pnm[:,inm1];
    for m in range(0, nmax-1):
        for n in range(m+2, nmax+1):
            inm1 = nm2ind(n-2, m);
            inm2 = nm2ind(n-1, m);
            inm3 = nm2ind(n, m);
            anm = (2*n-1)/np.sqrt((n-m)*(n+m));
            bnm = np.sqrt((n-m-1)*(n+m-1))/np.sqrt((n-m)*(n+m));
            Pnm[:,inm3] = anm * costh * Pnm[:,inm2] - bnm * Pnm[:,inm1];
    if (not flag_dP) and (not flag_mP):
        return Pnm;

def read_shc(filename):
    """ 
    Read .shc (spherical harmonic coefficient) file

    The function produces data frames that have time as index and 
    spherical harmonic degree and order as columns. In the case of IGRF,
    the times will correspond to the different models 5 years apart

    Parameters
    ----------
    filename : string
        filename of .shc file

    Returns
    -------
    g : DataFrame
        pandas DataFrame of gauss coefficients for cos terms.
    h : DataFrame
        pandas DataFrame of gauss coefficients for sin terms.

    Note
    ----
    This code has no special treatment of "secular variation" coefficients. 
    Instead, the SV coefficients of IGRF should be used to make gauss 
    coefficients. This must be done prior to this code (when making the 
    .shc file).
    """

    header = 2
    coeffdict = {}
    with open(filename, 'r') as f:
        for line in f.readlines():

            if line.startswith('#'): # this is a header that we don't read
                continue

            if header == 2: # read parameters (could be skipped...)
                 N_MIN, N_MAX, NTIMES, SP_ORDER, N_STEPS = list(map(int, line.split()[:5]))
                 header -= 1
                 continue

            if header == 1: # read years
                 times = yearfrac_to_datetime(list(map(float, line.split())) )
                 header -= 1
                 continue

            key = tuple(map(int, line.split()[:2]))
            coeffdict[key] = np.array(list(map(float, line.split()[2:])))

    g = {key:coeffdict[key] for key in coeffdict.keys() if key[1] >= 0}
    h = {(key[0], -key[1]):coeffdict[key] for key in coeffdict.keys() if key[1] < 0 }
    for key in [k for k in g.keys() if k[1] == 0]: # add zero coefficients for m = 0 in h dictionary
        h[key] = 0

    # this must be true:
    assert len(g.keys()) == len(h.keys())

    gdf = pd.DataFrame(g, index = times)
    hdf = pd.DataFrame(h, index = times)

    # make sure that the column keys of h are in same order as in g:
    hdf = hdf[gdf.columns]

    return gdf, hdf

def SHS_gc(Lamda, Theta, Radius, date, shc_file, min_degree=1, max_degree=13):
    """Spherical Harmonic Systhesis of magnetic models e.g. IGRF
    Input & Output in geocentric coordinates

    Args:
    <lamda, theta, r> : numpy array
        longitude [deg], colatitude [deg], radius [km]
    date : 
        one date to evaluate IGRF coefficients
    shc_file : string, optional
        filename of .shc file. Default is latest IGRF
    min_degree : int, optional
        lowest degree of expansion  min_degree >= 1
    max_degree : int, optional
        highest degree of expansion, 1 <= max_degree <= 13

    Returns:
    <Blamda, Btheta, Br> : numpy array
        Magnetic field [nT] in eastward, theta(south), radial directions
    """

    RE = 6371.2 # km
    # read coefficient file:
    g, h = read_shc(shc_file)

    date = np.array([date])

    if min_degree > max_degree:
        print('Warning: Highest degree of expansion must be larger or equal to lowest degree.')
        print('Reset to original range.')
        min_degree, max_degree = 1, 13

    inshp = Lamda.shape;
    Lamda, Theta, Radius = map(lambda x: x.flatten().reshape((-1 ,1)), [Lamda, Theta, Radius]) # column vectors

    # make row vectors of wavenumbers n and m:
    n, m = np.array([k for k in g.columns]).T
    n, m = n.reshape((1, -1)), m.reshape((1, -1))

    # get maximum N and maximum M:
    N, M = np.max(n), np.max(m)

    # get the legendre functions
    P, dP = get_legendre(Theta, g.keys())

    # Append coefficients at desired times (skip if index is already in coefficient data frame):
    index = g.index.union(date)

    g = g.reindex(index).groupby(index).first() # reindex and skip duplicates
    h = h.reindex(index).groupby(index).first() # reindex and skip duplicates

    # interpolate and collect the coefficients at desired times:
    g = g.interpolate(method = 'time').loc[date, :]
    h = h.interpolate(method = 'time').loc[date, :]

    # compute cosmlon and sinmlon:
    lamda_rad = np.radians(Lamda)
    cosmlamda = np.cos(lamda_rad * m) # shape (n_coords x n_model_params/2)
    sinmlamda = np.sin(lamda_rad * m)

    # make versions of n and m that are repeated twice
    nn, mm = np.tile(n, 2), np.tile(m, 2)

    N_map = ((nn >= min_degree) & (nn <= max_degree)).astype(int)

    # calculate Br:
    G  = N_map * (RE / Radius) ** (nn + 2) * (nn + 1) * np.hstack((P * cosmlamda, P * sinmlamda))
    Br = G.dot(np.hstack((g.values, h.values)).T).T # shape (n_times, n_coords)

    # calculate Btheta:
    G  = - N_map * (RE / Radius) ** (nn + 1) * np.hstack((dP * cosmlamda, dP * sinmlamda)) \
         * RE / Radius
    Btheta = G.dot(np.hstack((g.values, h.values)).T).T # shape (n_times, n_coords)

    # calculate Bphi:
    G  = - N_map * (RE / Radius) ** (nn + 1) * mm * np.hstack((-P * sinmlamda, P * cosmlamda)) \
         * RE / Radius / np.sin(np.radians(Theta))
    Blamda = G.dot(np.hstack((g.values, h.values)).T).T # shape (n_times, n_coords)

    # reshape and return
    return Blamda.reshape(inshp), Btheta.reshape(inshp), Br.reshape(inshp);

def SHS_gd(Lon, Lat, Hgt, date, shc_file, ellipsoid='wgs84', min_degree=1, max_degree=13):
    """Spherical Harmonic Systhesis of magnetic models e.g. IGRF
    Input & Output in geodetic coordinates

    Args:
        Lon (numpy array)           : longitude in degrees
        Lat (numpy array)           : geodetic latitude in degrees
        Hgt (numpy array)           : height above ellipsoid in km
        date (datetime)             : datetime to evaluate SHCs
        shc_file (string)           : filename of .shc file
        ellipsoid (str, optional)   : ellipsoid. Default to 'wgs84'
        min_degree (int, optional)  : lowest degree of expansion. Default to 1
        max_degree (int, optional)  : highest degree of expansion. Default to 13

    Returns:
        Bn, Be, Bd (numpy arrays)   : Magnetic field [nT] in <northward, eastward, downward> directions relative to ellipsoid
    """

    if min_degree > max_degree:
        print('Warning: Highest degree of expansion must be larger or equal to lowest degree.')
        print('Reset to original range.')
        min_degree, max_degree = 1, 13

    # Geodetic to Geocentric
    a, _, _, ecc = GetEllipsoid(ellipsoid);
    X, Y, Z = ell2ecef(Lon, Lat, Hgt, a, ecc, degrees=True);
    Lamda, Theta, Radius = ecef2sph(X, Y, Z, colat=True, degrees=True);
    
    # SHS to calculate geocentric components 
    Blamda, Btheta, Br = SHS_gc(Lamda, Theta, Radius, date, shc_file, min_degree=min_degree, max_degree=max_degree);
    
    # Geocentric to Geodetic
    Bx, By, Bz = ned2ecef_v(Lamda, Theta, -Btheta, Blamda, -Br, colat=True, degrees=True);
    Bn, Be, Bd = ecef2ned_v(X, Y, Z, Bx, By, Bz, a, ecc);

    return Bn, Be, Bd;

def ComputeMagB7(Bn, Be, Bd, degrees=True):
    """Compute the 7 magnetic field components

    Args:
        Bn (numpy array): North component
        Be (numpy array): East component
        Bd (numpy array): Downward component
        degrees (bool, optional): Defaults to True.

    Returns:
        BX, BY, BZ, BD, BI, BH, BF: The 7 magnetic field components
    """
    BX = Bn; BY = Be; BZ = Bd;
    BH = np.sqrt(BX**2 + BY**2);
    BF = np.sqrt(BX**2 + BY**2 +BZ**2);
    BI = np.arctan(BZ/BH);
    BD = np.arcsin(BY/BH);  
    BI[BH == 0] = np.sign(BZ[BH == 0]) * np.pi /2;
    BD[BH == 0] = 0;
    if degrees:
        BI = np.rad2deg(BI);
        BD = np.rad2deg(BD);
    return BX, BY, BZ, BD, BI, BH, BF

def SaveMagB7File(file_MagB7, LON, LAT, HGT, BX, BY, BZ, BD, BI, BH, BF, coord='ell', showinfo=True):
    # ***************** Format ******************
    # LON, LAT, HGT, BX, BY, BZ, BD, BI, BH, BF
    # *******************************************
    [LON, LAT, HGT, BX, BY, BZ, BD, BI, BH, BF] = map(lambda x: x.flatten(),\
        [LON, LAT, HGT, BX, BY, BZ, BD, BI, BH, BF]);
    fid = open(file_MagB7, 'w');
    wrtArray=np.column_stack((LON,LAT,HGT,BX,BY,BZ,BD,BI,BH,BF));
    wrtFormat=['%12.6f','%12.6f','%12.4f','%12.1f','%12.1f','%12.1f','%12.6f','%12.6f','%12.1f','%12.1f'];
    if coord=='ell':
        colHead='%10s %12s %12s %12s %12s %12s %12s %12s %12s %12s'%('Lon','Lat','H','B_X','B_Y','B_Z','B_Dec','B_Inc','B_H','B_F');
    elif coord=='sph':
        colHead='%10s %12s %12s %12s %12s %12s %12s %12s %12s %12s'%('Lamda','Psi','R','B_X','B_Y','B_Z','B_Dec','B_Inc','B_H','B_F');
    np.savetxt(fid,wrtArray,fmt=wrtFormat,header=colHead);
    fid.close();
    print('\n*****************************************\n',flush=True);
    print('Saving estimation of B7 components to file: %s'%(file_MagB7),flush=True);
    if showinfo:
        print('BX in [%12.1f,%12.1f] nT;'%(BX.min(),BX.max()),flush=True);
        print('BY in [%12.1f,%12.1f] nT;'%(BY.min(),BY.max()),flush=True);
        print('BZ in [%12.1f,%12.1f] nT;'%(BZ.min(),BZ.max()),flush=True);
        print('BD in [%12.6f,%12.6f] deg;'%(BD.min(),BD.max()),flush=True);
        print('BI in [%12.6f,%12.6f] deg;'%(BI.min(),BI.max()),flush=True);
        print('BH in [%12.1f,%12.1f] nT;'%(BH.min(),BH.max()),flush=True);
        print('BF in [%12.1f,%12.1f] nT;'%(BF.min(),BF.max()),flush=True);
    print('\n*****************************************\n',flush=True);
    return None;
