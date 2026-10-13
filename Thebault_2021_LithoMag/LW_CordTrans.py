##################################################################
# # ! Leyuan Wu
# Coordinate transformation
# Transform 3-D coordinates between 
#               Earth-centered Earth-fixed (ECEF), 
#               geodetic, 
#               east-north-up (ENU), 
#               north-east-down (NED), 
#               azimuth-elevation-range (AER) 
# systems. Transform 3-D vector components between ECEF, ENU, and NED systems.
##################################################################

import numpy as np;

def GetEllipsoid(Name):
    """Get reference ellipsoid parameters

    Args:
        Name (str): name of ellipsoid

    Raises:
        NotImplementedError: name not find
    
    Units:
        Distance in km

    Returns:
        a, b, f, ecc (float): semimajor, semiminor, flattening & eccentricity 
    """
    if Name == "wgs84":
        # * https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        a = 6378137.0/1e3
        b = 6356752.31424518/1e3
    elif Name == "grs80":
        # * https://en.wikipedia.org/wiki/GRS_80
        a = 6378137.0/1e3
        b = 6356752.31414036/1e3
    elif Name == "mars":
        # * https://tharsis.gsfc.nasa.gov/geodesy.html
        a = 3396900.0/1e3
        b = 3376097.80585952/1e3
    elif Name == "moon":
        a = 1738000.0/1e3
        b = a
    elif Name == "venus":
        a = 6051000.0/1e3
        b = a
    elif Name == "jupiter":
        a = 71492000.0/1e3
        b = 66770054.3475922/1e3
    elif Name == "io":
        # * https://doi.org/10.1006/icar.1998.5987
        a = 1829.7/1e3
        b = 1815.8/1e3
    elif Name == "pluto":
        a = 1187000.0/1e3
        b = a
    else:
        raise NotImplementedError(
            "{} not find".format(Name)
        )
    f = (a - b) / a;
    ecc = np.sqrt(2 * f - f ** 2);
    return a, b, f, ecc;

def ell2ecef(Lon, Lat, H, a, ecc, colat=False, degrees=False):
    """Ellipsoidal to ECEF Cartesian coordinates
    Args:
        Lon (numpy array)       : longitudes in radians
        Lat (numpy array)       : latitudes in radians
        H (numpy array)         : ellipsoidal height in km
        a (float)               : major axis of spheroid in km
        ecc (float)             : eccentricity
        colat (bool, optional)  : co-Latitude. defaults to False.
        degrees (bool, optional): degree or radian, defaults to False (radian).
        
    Returns:
        <X, Y, Z> (numpy array) : ECEF Cartesian coordinates in km
    """
    if degrees:
        Lon = np.deg2rad(Lon);
        Lat = np.deg2rad(Lat);
    if colat:
        Lat = np.pi/2-Lat;
    e2 = ecc**2;
    N = a/np.sqrt(1-e2*np.sin(Lat)**2);
    X = (N+H)*(np.cos(Lat)*np.cos(Lon));
    Y = (N+H)*(np.cos(Lat)*np.sin(Lon));
    Z = (N*(1-e2)+H)*np.sin(Lat);
    return X, Y, Z;

def ecef2ell(X, Y, Z, a, ecc, colat=False, degrees=False):
    """From Earth-Centered-Earth-Fixed to 
        Ellipsoidal (geodetic) coordinates using Bowring Algorithm
    Reference:
    Bowring, B.R., 1976. Transformation from spatial to geographical coordinates.
    Surv. Rev. 23 (181), 323-327.

    Args:
        <X, Y, Z> (numpy array)     : ECEF Coordinates in km
        a (float)                   : major axis of spheroid in km
        ecc (float)                 : eccentricity
        colat (bool, optional)      : co-Latitude. defaults to False.
        degrees (bool, optional)    : degree or radian, defaults to False (radian).

    Returns:
        <Lon, Lat, H> (numpy array) : Geodetic Coordinates
    """
    ecc2 = ecc**2;              # eccentricity squared
    b = a*np.sqrt(1-ecc2);      # semi-minor axis   
    Lon = np.arctan2(Y,X);      # Geographic longitude
    # Auxiliary values
    p = np.sqrt(X**2 + Y**2);
    mu = np.arctan(a*Z/(b*p));
    ea2 = ecc2/(1-ecc2);
    # Geographic latitude
    Lat = np.arctan2(Z + ea2*b*np.sin(mu)**3, p - ecc2*a*np.cos(mu)**3);
    # Height above ellipsoid
    H = p/np.cos(Lat) - a/np.sqrt(1-ecc2*np.sin(Lat)**2);
    # Lat to co-Lat
    if colat:
        Lat = np.pi/2-Lat;
    # Radians to degrees
    if degrees:
        Lon = np.rad2deg(Lon);
        Lat = np.rad2deg(Lat);
    return Lon, Lat, H;

def sph2ecef(Lon, Lat, R, colat=False, degrees=False):
    """Spherical to ECEF Cartesian coordinates

    Args:
        Lon (numpy array)       : longitudes in radians
        Lat (numpy array)       : latitudes in radians
        R (numpy array)         : radius in km
        colat (bool, optional)  : co-Latitude. defaults to False.
        degrees (bool, optional): degree or radian, defaults to False (radian).

    Returns:
        <X, Y, Z> (numpy array) : ECEF Cartesian coordinates in km
    """
    if degrees:
        Lon = np.deg2rad(Lon);
        Lat = np.deg2rad(Lat);
    if colat:
        Lat = np.pi/2-Lat;
    Z = R * np.sin(Lat);
    temp = R * np.cos(Lat);
    X = temp * np.cos(Lon);
    Y = temp * np.sin(Lon);
    return X, Y, Z;

def ecef2sph(X, Y, Z, colat=False, degrees=False):
    """ECEF to Spherical coordinates 
    Args:
        <X, Y, Z> (numpy array)     : ECEF Coordinates in km
        colat (bool, optional)      : co-Latitude. defaults to False.
        degrees (bool, optional)    : degree or radian, defaults to False (radian).

    Returns:
        <Lon, Lat, R> (numpy array) : Spherical coordinates
    """
    R=np.sqrt(X**2+Y**2+Z**2);
    Lon=np.arctan2(Y,X);
    Lat=np.arccos(Z/R); # ! Co-latitude in [0,pi]
    Lat = np.pi/2-Lat;  # ! latitude in [-pi/2,pi/2]
    if colat:
        Lat  = np.pi/2-Lat;
    if degrees:
        Lon = np.rad2deg(Lon);
        Lat = np.rad2deg(Lat);
    return Lon, Lat, R;

def ell2ned(Lon, Lat, H, lon0, lat0, h0, a, ecc, colat=False, degrees=False):
    """Transforming geodetic to North-East-Down (NED) coordinates
    Step 1: Convert Geodetic Coordinates to ECEF
    Step 2: Calculate the Relative ECEF Coordinates
    Step 3: Transform to NED Coordinates

    Args:
        <Lon, Lat, H> (numpy array)     : Geodetic Coordinates (H in km)
        <lon0, lat0, h0> (float)        : Reference point
        a (float)                       : semi-major in km
        ecc (float)                     : eccentricity
        colat (bool, optional)          : co-Latitude. defaults to False.
        degrees (bool, optional)        : degree or radian, defaults to False (radian).

    Returns:
        <N, E, D> (numpy array)         : North-East-Down (NED) coordinates
    """
    if degrees:
        Lon = np.deg2rad(Lon);
        Lat = np.deg2rad(Lat);
        lon0 = np.deg2rad(lon0);
        lat0 = np.deg2rad(lat0); 
    if colat:
        lat0 = np.pi/2-lat0;
        Lat  = np.pi/2-Lat;
    x0, y0, z0 = ell2ecef(lon0, lat0, h0, a, ecc);
    X, Y, Z = ell2ecef(Lon, Lat, H, a, ecc);
    Dx = X - x0; Dy = Y - y0; Dz = Z - z0;
    uN = np.array([-np.sin(lat0)*np.cos(lon0), -np.sin(lat0)*np.sin(lon0),  np.cos(lat0)]);
    uE = np.array([-np.sin(lon0)             ,  np.cos(lon0)             ,  0]);
    uD = np.array([-np.cos(lat0)*np.cos(lon0), -np.cos(lat0)*np.sin(lon0), -np.sin(lat0)]);
    ######## * Old algorithm
    # P = np.hstack((Dx.reshape(-1,1), Dy.reshape(-1,1),Dz.reshape(-1,1)))
    # R = np.vstack((uN, uE, uD));
    # NED = P @ np.transpose(R);
    # return np.reshape(NED[:,0],Lon.shape), np.reshape(NED[:,1],Lon.shape), np.reshape(NED[:,2],Lon.shape);
    ######## * New algorithm
    N = uN[0]*Dx + uN[1]*Dy + uN[2]*Dz;
    E = uE[0]*Dx + uE[1]*Dy + uE[2]*Dz;
    D = uD[0]*Dx + uD[1]*Dy + uD[2]*Dz;
    return N, E, D;

def ned2ell(N, E, D, lon0, lat0, h0, a, ecc, colat=False, degrees=False):
    """Transforming North-East-Down (NED) to geodetic coordinates
    Step 1: Calculate the ECEF Coordinates of the Reference (Origin) point
    Step 2: Add to the reference point NED=>ECEF coponents
    Step 3: Transform to geodetic coordinates

    Args:
        <N, E, D> (numpy array)         : North-East-Down (NED) coordinates in km
        <lon0, lat0, h0> (float)        : Reference point
        a (float)                       : semi-major in km
        ecc (float)                     : eccentricity
        colat (bool, optional)          : co-Latitude. defaults to False.
        degrees (bool, optional)        : degree or radian, defaults to False (radian).

    Returns:
        <Lon, Lat, H> (numpy array)     : Geodetic Coordinates (H in km)
    """
    if degrees:
        lon0 = np.deg2rad(lon0);
        lat0 = np.deg2rad(lat0); 
    if colat:
        lat0 = np.pi/2-lat0;
    x0, y0, z0 = ell2ecef(lon0, lat0, h0, a, ecc);
    uN = np.array([-np.sin(lat0)*np.cos(lon0), -np.sin(lat0)*np.sin(lon0),  np.cos(lat0)]);
    uE = np.array([-np.sin(lon0)             ,  np.cos(lon0)             ,  0]);
    uD = np.array([-np.cos(lat0)*np.cos(lon0), -np.cos(lat0)*np.sin(lon0), -np.sin(lat0)]);
    Dx = uN[0]*N + uE[0]*E + uD[0]*D;
    Dy = uN[1]*N + uE[1]*E + uD[1]*D;
    Dz = uN[2]*N + uE[2]*E + uD[2]*D;
    X = Dx + x0; Y = Dy + y0; Z = Dz + z0;
    Lon, Lat, H = ecef2ell(X, Y, Z, a, ecc, colat=colat, degrees=degrees);
    return Lon, Lat, H;

def sph2planar(Lon, Lat, RE, lon0, lat0, colat=False, degrees=False):
    """Spherical to Planar coordinates

    Args:
        <Lon, Lat> (numpy array)    : longitudes/latitudes in radians
        RE (float)                  : Earth radius in km
        <lon0, lat0> (float)        : Origin of the planar approximation
        colat (bool, optional)      : co-Latitude. defaults to False.
        degrees (bool, optional)    : degree or radian, defaults to False (radian).

    Returns:
        Est, Nor                    : Local Cartesian East (Y) and North (X) coordinates
    """
    if degrees:
        Lon = np.deg2rad(Lon);
        Lat = np.deg2rad(Lat);
        lon0 = np.deg2rad(lon0);
        lat0 = np.deg2rad(lat0); 
    if colat:
        lat0 = np.pi/2-lat0;
        Lat  = np.pi/2-Lat;
    dEst_rad = RE*np.cos(lat0);
    dNor_rad = RE;
    Est = (Lon-lon0)*dEst_rad;
    Nor = (Lat-lat0)*dNor_rad;
    return Est, Nor;

def ecef2ned_v(X, Y, Z, U, V, W, a, ecc):
    """Transform vector from Earth-Centered-Earth-Fixed XYZ-components
    to Ellipsoidal (geodetic) North-East-Downward (NED) components 

    Args:
        X,Y,Z (numpy array): ECEF Coordinates
        U,V,W (numpy array): ECEF Vector components
        a (float): semi-major
        ecc (float): eccentricity

    Returns:
        Gn,Ge,Gd: North-East-Downward (NED) components
    """
    assert(X.shape==Y.shape==Z.shape==U.shape==V.shape==W.shape);
    ######## * Old algorithm
    # shpComp=X.shape;
    # X=X.reshape(-1,1);Y=Y.reshape(-1,1);Z=Z.reshape(-1,1);
    # U=U.reshape(-1,1);V=V.reshape(-1,1);W=W.reshape(-1,1);
    # Lon,Lat,_ = ecef2ell(X,Y,Z,a,ecc,colat=False);

    # unitx=np.array([1,0,0]).reshape(-1,1);
    # unity=np.array([0,1,0]).reshape(-1,1);
    # unitz=np.array([0,0,1]).reshape(-1,1);

    # N=np.hstack((-np.sin(Lat)*np.cos(Lon),-np.sin(Lat)*np.sin(Lon),np.cos(Lat)));
    # E=np.hstack((-np.sin(Lon),np.cos(Lon),np.zeros((Lon.shape))));
    # D=np.hstack((-np.cos(Lat)*np.cos(Lon),-np.cos(Lat)*np.sin(Lon),-np.sin(Lat)));

    # Gn=U*(N@unitx)+V*(N@unity)+W*(N@unitz);
    # Ge=U*(E@unitx)+V*(E@unity)+W*(E@unitz);
    # Gd=U*(D@unitx)+V*(D@unity)+W*(D@unitz);
    
    # Gn=np.reshape(Gn,shpComp);
    # Ge=np.reshape(Ge,shpComp);
    # Gd=np.reshape(Gd,shpComp);
    ######## * New algorithm
    Lon, Lat, _ = ecef2ell(X,Y,Z,a,ecc);
    N = [-np.sin(Lat)*np.cos(Lon), -np.sin(Lat)*np.sin(Lon),  np.cos(Lat)];
    E = [-np.sin(Lon),              np.cos(Lon),              np.zeros((Lon.shape))];
    D = [-np.cos(Lat)*np.cos(Lon), -np.cos(Lat)*np.sin(Lon), -np.sin(Lat)];

    Gn = N[0]*U + N[1]*V + N[2]*W;
    Ge = E[0]*U + E[1]*V + E[2]*W;
    Gd = D[0]*U + D[1]*V + D[2]*W;
    
    return Gn, Ge, Gd;

def ned2ecef_v(Lon, Lat, Gn, Ge, Gd, colat=False, degrees=False):
    """Transform vector from 
    Ellipsoidal (geodetic) North-East-Downward (NED) components 
    to Earth-Centered-Earth-Fixed XYZ components

    Args:
        <Lat, Lon> (numpy array)        : Lat and Lon coordinates
        <Gn, Ge, Gd> (numpy array)      : (NED) components
        colat (bool, optional)          : co-Latitude. defaults to False.
        degrees (bool, optional)        : degree or radian, defaults to False (radian).
    Returns:
        <Gx, Gy, Gz>                    : ECEF Vector components
    """
    assert(Lon.shape==Lat.shape==Gn.shape==Ge.shape==Gd.shape);
    if degrees:
        Lon = np.deg2rad(Lon);
        Lat = np.deg2rad(Lat);
    if colat:
        Lat = np.pi/2-Lat;
    ######## * Old algorithm
    # shpComp=Lat.shape;
    # Lat=Lat.reshape(-1,1);Lon=Lon.reshape(-1,1);
    # Gn=Gn.reshape(-1,1);Ge=Ge.reshape(-1,1);Gd=Gd.reshape(-1,1);
    # unitx=np.array([1,0,0]).reshape(-1,1);
    # unity=np.array([0,1,0]).reshape(-1,1);
    # unitz=np.array([0,0,1]).reshape(-1,1);

    # N=np.hstack((-np.sin(Lat)*np.cos(Lon),-np.sin(Lat)*np.sin(Lon),np.cos(Lat)));
    # E=np.hstack((-np.sin(Lon),np.cos(Lon),np.zeros((Lon.shape))));
    # D=np.hstack((-np.cos(Lat)*np.cos(Lon),-np.cos(Lat)*np.sin(Lon),-np.sin(Lat)));

    # Gx=Gn*(N@unitx)+Ge*(E@unitx)+Gd*(D@unitx);
    # Gy=Gn*(N@unity)+Ge*(E@unity)+Gd*(D@unity);
    # Gz=Gn*(N@unitz)+Ge*(E@unitz)+Gd*(D@unitz);
    
    # Gx=np.reshape(Gx,shpComp);
    # Gy=np.reshape(Gy,shpComp);
    # Gz=np.reshape(Gz,shpComp);
    ######## * New algorithm
    N = [-np.sin(Lat)*np.cos(Lon), -np.sin(Lat)*np.sin(Lon),  np.cos(Lat)];
    E = [-np.sin(Lon),              np.cos(Lon),              np.zeros((Lon.shape))];
    D = [-np.cos(Lat)*np.cos(Lon), -np.cos(Lat)*np.sin(Lon), -np.sin(Lat)];

    Gx = N[0]*Gn + E[0]*Ge + D[0]*Gd;
    Gy = N[1]*Gn + E[1]*Ge + D[1]*Gd;
    Gz = N[2]*Gn + E[2]*Ge + D[2]*Gd;
    return Gx, Gy, Gz;

def MD2ecef(MDx,MDy,MDz,md_lon,md_lat,systemIII=True):
    """Magnetic Dipole to ECEF Cartesian coordinates
    Args:
        MDx,MDy,MDz (numpy array):  Magnetic Dipole coordinates
        md_lon: Magnetic dipole Longitude
        md_lat: Magnetic dipole Latitude
        systemIII (bool, optional): . Defaults to True.

    Returns:
        x,y,z: (numpy array):  ECEF coordinates
    """
    if(systemIII):
        md_lon = 360.0-md_lon;
    md_lon = np.deg2rad(md_lon);
    md_lat = np.deg2rad(md_lat);
    unitMDx = np.array([np.cos(md_lat)*np.cos(md_lon),np.cos(md_lat)*np.sin(md_lon),-np.sin(md_lat)]).reshape(1,-1);
    unitMDy = np.array([-np.sin(md_lon),np.cos(md_lon),0]).reshape(1,-1);
    unitMDz = np.array([np.sin(md_lat)*np.cos(md_lon),np.sin(md_lat)*np.sin(md_lon),np.cos(md_lat)]).reshape(1,-1);
    MD_Pts = np.hstack((MDx,MDy,MDz));
    Pts = MD_Pts@np.vstack((unitMDx, unitMDy, unitMDz));
    return Pts[:,0:1], Pts[:,1:2], Pts[:,2:3];

def MD2ecef_v(MDu,MDv,MDw,md_lon,md_lat,systemIII=True):
    """Magnetic Dipole to ECEF Cartesian vectors
    Args:
        MDu,MDv,MDw (numpy array):  Magnetic Dipole vector components
        md_lon: Magnetic dipole Longitude
        md_lat: Magnetic dipole Latitude
        systemIII (bool, optional): . Defaults to True.
    Returns:
        u,v,w,norm: (numpy array):  ECEF vector components and magnitude
    """
    if(systemIII):
        md_lon = 360.0-md_lon;
    md_lon = np.deg2rad(md_lon);
    md_lat = np.deg2rad(md_lat);
    unitMDx = np.array([np.cos(md_lat)*np.cos(md_lon),np.cos(md_lat)*np.sin(md_lon),-np.sin(md_lat)]).reshape(1,-1);
    unitMDy = np.array([-np.sin(md_lon),np.cos(md_lon),0]).reshape(1,-1);
    unitMDz = np.array([np.sin(md_lat)*np.cos(md_lon),np.sin(md_lat)*np.sin(md_lon),np.cos(md_lat)]).reshape(1,-1);
    MD_VECs = np.hstack((MDu,MDv,MDw));
    VECs = MD_VECs@np.vstack((unitMDx, unitMDy, unitMDz));
    Norms = np.sqrt(np.sum(VECs**2,axis=1,keepdims=True));
    return VECs[:,0:1], VECs[:,1:2], VECs[:,2:3], Norms;

def ecef2MD(x,y,z,md_lon,md_lat,systemIII=True):
    """ECEF to Magnetic Dipole Cartesian coordinates
    Args:
        x,y,z (numpy array):  ECEF coordinates
        md_lon: Magnetic dipole Longitude
        md_lat: Magnetic dipole Latitude
        systemIII (bool, optional): . Defaults to True.
    Returns:
        MDx,MDy,MDz: (numpy array):  Magnetic Dipole coordinates
    """
    if(systemIII):
        md_lon = 360.0-md_lon;
    md_lon = np.deg2rad(md_lon);
    md_lat = np.deg2rad(md_lat);
    unitMDx = np.array([np.cos(md_lat)*np.cos(md_lon),np.cos(md_lat)*np.sin(md_lon),-np.sin(md_lat)]).reshape(-1,1);
    unitMDy = np.array([-np.sin(md_lon),np.cos(md_lon),0]).reshape(-1,1);
    unitMDz = np.array([np.sin(md_lat)*np.cos(md_lon),np.sin(md_lat)*np.sin(md_lon),np.cos(md_lat)]).reshape(-1,1);
    Pts = np.hstack((x,y,z));
    MD_Pts = Pts@np.hstack((unitMDx, unitMDy, unitMDz));
    return MD_Pts[:,0:1], MD_Pts[:,1:2], MD_Pts[:,2:3];

def ecef2MD_v(u,v,w,md_lon,md_lat,systemIII=True):
    """ECEF to Magnetic Dipole Cartesian vectors
    Args:
        (u,v,w) (numpy array):  ECEF vector components
        md_lon: Magnetic dipole Longitude
        md_lat: Magnetic dipole Latitude
        systemIII (bool, optional): . Defaults to True.
    Returns:
        MDu,MDv,MDw,MD_Norm: (numpy array):  Magnetic Dipole vector components and magnitude
    """
    if(systemIII):
        md_lon = 360.0-md_lon;
    md_lon = np.deg2rad(md_lon);
    md_lat = np.deg2rad(md_lat);
    unitMDx = np.array([np.cos(md_lat)*np.cos(md_lon),np.cos(md_lat)*np.sin(md_lon),-np.sin(md_lat)]).reshape(-1,1);
    unitMDy = np.array([-np.sin(md_lon),np.cos(md_lon),0]).reshape(-1,1);
    unitMDz = np.array([np.sin(md_lat)*np.cos(md_lon),np.sin(md_lat)*np.sin(md_lon),np.cos(md_lat)]).reshape(-1,1);
    VECs = np.hstack((u,v,w));
    MD_VECs = VECs@np.hstack((unitMDx, unitMDy, unitMDz));
    MD_Norm = np.sqrt(np.sum(MD_VECs**2,axis=1,keepdims=True));
    return MD_VECs[:,0:1], MD_VECs[:,1:2], MD_VECs[:,2:3], MD_Norm;