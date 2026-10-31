"""
magnetic_forward_numba.py
=========================
3D Forward modeling of magnetic fields.

Table of Contents
-----------------
1. Analytical Solution (Magnetic Polyhedron)
2. Numerical Quadrature (Gauss-Legendre Quadrature)
3. FFT-Based Forward Modeling (Fast Multipole Method)
4. Coordinate Transformation
5. Miscellaneous Utilities
"""

import numpy as np
from numba import njit, prange


# =============================================================================
# 1. ANALYTICAL SOLUTION
# =============================================================================

# -----------------------------------------------
# 1.1 Dipole
# -----------------------------------------------

@njit(parallel=True)
def dipole_field(obs, m, r0, mu0_over_4pi=1e-7):
    """
    Compute the magnetic field B of a single dipole at observation points.

    The field is evaluated using the analytical dipole expression:

    .. math::
        \mathbf{B} = \frac{\mu_0}{4\pi}
        \left[ \frac{3(\mathbf{m} \cdot \hat{\mathbf{r}})\hat{\mathbf{r}} - \mathbf{m}}{r^3} \right]

    Parameters
    ----------
    obs : (N, 3) ndarray
        Array of observation point coordinates.
    m : (3,) ndarray
        Dipole moment vector [A·m²].
    r0 : (3,) ndarray
        Dipole source position.
    mu0_over_4pi : float, optional
        Value of :math:`\mu_0 / (4\pi)` in H/m. Default is 1e-7.

    Returns
    -------
    Bx, By, Bz : (N,) ndarray
        Magnetic field components at the observation points [nT].
    """
    N = obs.shape[0]
    B = np.empty((N, 3), dtype=np.float64)

    for i in prange(N):
        rx = obs[i, 0] - r0[0]
        ry = obs[i, 1] - r0[1]
        rz = obs[i, 2] - r0[2]

        r2 = rx * rx + ry * ry + rz * rz
        r = np.sqrt(r2)
        r3 = r2 * r

        # Avoid singularity at the dipole location
        if r3 < 1e-15:
            B[i, 0] = 0.0
            B[i, 1] = 0.0
            B[i, 2] = 0.0
            continue

        # Unit vector r_hat
        rx_hat = rx / r
        ry_hat = ry / r
        rz_hat = rz / r

        # Dot product m · r_hat
        mdotr = m[0] * rx_hat + m[1] * ry_hat + m[2] * rz_hat

        # B = mu0/(4pi) * (3*(m·r_hat)*r_hat - m) / r^3
        coeff = 1.e9 * mu0_over_4pi / r3
        B[i, 0] = coeff * (3.0 * mdotr * rx_hat - m[0])
        B[i, 1] = coeff * (3.0 * mdotr * ry_hat - m[1])
        B[i, 2] = coeff * (3.0 * mdotr * rz_hat - m[2])

    return B[:, 0], B[:, 1], B[:, 2]


@njit(parallel=True)
def dipole_potential(obs, m, r0, mu0_over_4pi=1e-7):
    """
    Compute the scalar magnetic potential V of a single dipole.

    The potential is given by:

    .. math::
        V = \frac{\mu_0}{4\pi} \frac{\mathbf{m} \cdot \hat{\mathbf{r}}}{r^2}
        = \frac{\mu_0}{4\pi} \frac{\mathbf{m} \cdot \mathbf{r}}{r^3}

    Parameters
    ----------
    obs : (N, 3) ndarray
        Array of observation point coordinates.
    m : (3,) ndarray
        Dipole moment vector [A·m²].
    r0 : (3,) ndarray
        Dipole source position.
    mu0_over_4pi : float, optional
        Value of :math:`\mu_0 / (4\pi)` in H/m. Default is 1e-7.

    Returns
    -------
    V : (N,) ndarray
        Scalar magnetic potential at the observation points [A].
    """
    N = obs.shape[0]
    V = np.empty(N, dtype=np.float64)

    for i in prange(N):
        rx = obs[i, 0] - r0[0]
        ry = obs[i, 1] - r0[1]
        rz = obs[i, 2] - r0[2]

        r2 = rx * rx + ry * ry + rz * rz
        r = np.sqrt(r2)
        r3 = r2 * r

        if r3 < 1e-15:
            V[i] = 0.0
            continue

        mdotr = m[0] * rx + m[1] * ry + m[2] * rz
        V[i] = mu0_over_4pi * mdotr / r3

    return V


@njit(parallel=True)
def dipole_tmi(obs, M, inc_m, dec_m, r0, inc_bg, dec_bg,
               mu0_over_4pi=1e-7, degrees=True):
    """
    Compute the total magnetic intensity (TMI) anomaly of a single dipole.

    The TMI anomaly is the projection of the magnetic anomaly vector onto
    the direction of the background (main) magnetic field:

    .. math::
        \Delta T = \mathbf{B}_{\text{anomaly}} \cdot \hat{\mathbf{B}}_0

    Parameters
    ----------
    obs : (N, 3) ndarray
        Array of observation point coordinates (NED: North, East, Down).
    M : float
        Magnitude of the dipole moment [A·m²].
    inc_m : float
        Inclination of the dipole moment, positive downward from horizontal.
    dec_m : float
        Declination of the dipole moment, clockwise from North.
    r0 : (3,) ndarray
        Dipole source position.
    inc_bg : float
        Inclination of the background field, positive downward from horizontal.
    dec_bg : float
        Declination of the background field, clockwise from North.
    mu0_over_4pi : float, optional
        Value of :math:`\\mu_0 / (4\\pi)` in H/m. Default is 1e-7.
    degrees : bool, optional
        If True (default), *inc_m*, *dec_m*, *inc_bg* and *dec_bg* are
        interpreted as degrees. If False, they are interpreted as radians.

    Returns
    -------
    TMI : (N,) ndarray
        Total magnetic intensity anomaly at the observation points [nT].
    """
    # Convert angles to radians if necessary
    if degrees:
        inc_m_rad = np.deg2rad(inc_m)
        dec_m_rad = np.deg2rad(dec_m)
        inc_bg_rad = np.deg2rad(inc_bg)
        dec_bg_rad = np.deg2rad(dec_bg)
    else:
        inc_m_rad = inc_m
        dec_m_rad = dec_m
        inc_bg_rad = inc_bg
        dec_bg_rad = dec_bg

    # Dipole moment vector in NED (North, East, Down)
    mx = M * np.cos(inc_m_rad) * np.cos(dec_m_rad)  
    my = M * np.cos(inc_m_rad) * np.sin(dec_m_rad) 
    mz = M * np.sin(inc_m_rad)        

    # Background field unit vector in NED
    b0x = np.cos(inc_bg_rad) * np.cos(dec_bg_rad)
    b0y = np.cos(inc_bg_rad) * np.sin(dec_bg_rad)
    b0z = np.sin(inc_bg_rad)

    N = obs.shape[0]
    TMI = np.empty(N, dtype=np.float64)

    for i in prange(N):
        rx = obs[i, 0] - r0[0]
        ry = obs[i, 1] - r0[1]
        rz = obs[i, 2] - r0[2]

        r2 = rx * rx + ry * ry + rz * rz
        r = np.sqrt(r2)
        r3 = r2 * r

        if r3 < 1e-15:
            TMI[i] = 0.0
            continue

        rx_hat = rx / r
        ry_hat = ry / r
        rz_hat = rz / r

        mdotr = mx * rx_hat + my * ry_hat + mz * rz_hat

        coeff = 1.e9 * mu0_over_4pi / r3

        Bx = coeff * (3.0 * mdotr * rx_hat - mx)
        By = coeff * (3.0 * mdotr * ry_hat - my)
        Bz = coeff * (3.0 * mdotr * rz_hat - mz)

        TMI[i] = Bx * b0x + By * b0y + Bz * b0z

    return TMI


# =============================================================================
# 2. NUMERICAL QUADRATURE
# =============================================================================



# =============================================================================
# 3. FFT-based Methods
# =============================================================================

# =============================================================================
# 4. Coordinate Transformations
# =============================================================================
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

# =============================================================================
# 5. Miscellaneous Utilities
# =============================================================================