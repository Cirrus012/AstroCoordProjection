# gnomonic_projection.py
"""
Gnomonic projection and its inverse.
This module provides functions for projecting spherical coordinates onto a tangent plane,
and recovering spherical coordinates from plane projections.

All angles (longitude/latitude, projection center) should be given in degrees.
Projection center and points to be projected must be in the same spherical coordinate system.
"""

import numpy as np

def gnomonic_projection(center_lon, center_lat, lon, lat):
    """
    Project spherical coordinates onto a tangent plane using gnomonic projection.

    Parameters
    ----------
    center_lon : float
        Longitude of projection center (degrees)
    center_lat : float
        Latitude of projection center (degrees)
    lon : array-like
        Longitudes of points to project (degrees)
    lat : array-like
        Latitudes of points to project (degrees)

    Returns
    -------
    xi : ndarray
        X coordinates on tangent plane (unitless, can be interpreted as "radians" for small angles)
    eta : ndarray
        Y coordinates on tangent plane (unitless, can be interpreted as "radians" for small angles)
    """
    lon = np.asarray(lon)
    lat = np.asarray(lat)

    dlon = np.radians((lon - center_lon + 180) % 360 - 180)
    sin_lat0 = np.sin(np.radians(center_lat))
    cos_lat0 = np.cos(np.radians(center_lat))
    sin_lat = np.sin(np.radians(lat))
    cos_lat = np.cos(np.radians(lat))

    denom = sin_lat * sin_lat0 + cos_lat * cos_lat0 * np.cos(dlon)

    xi = cos_lat * np.sin(dlon) / denom
    eta = (sin_lat * cos_lat0 - cos_lat * sin_lat0 * np.cos(dlon)) / denom
    return xi, eta

def inverse_gnomonic_projection(center_lon, center_lat, xi, eta):
    """
    Inverse gnomonic projection: recover spherical coordinates from plane coordinates.

    Parameters
    ----------
    center_lon : float
        Longitude of projection center (degrees)
    center_lat : float
        Latitude of projection center (degrees)
    xi : array-like
        X coordinates on tangent plane (unitless)
    eta : array-like
        Y coordinates on tangent plane (unitless)

    Returns
    -------
    lon : ndarray
        Longitudes of recovered points (degrees)
    lat : ndarray
        Latitudes of recovered points (degrees)
    """
    xi = np.asarray(xi, dtype=float)
    eta = np.asarray(eta, dtype=float)
    lat0_rad = np.radians(center_lat)

    sin_lat0 = np.sin(lat0_rad)
    cos_lat0 = np.cos(lat0_rad)

    # Compute denominator for both longitude and latitude
    denom = cos_lat0 - eta * sin_lat0

    # Longitude (in degrees, normalized to [0, 360))
    lam = np.arctan2(xi, denom)
    lon = (np.degrees(lam) + center_lon) % 360

    # Latitude (in degrees)
    numerator = (eta * cos_lat0 + sin_lat0) * np.cos(lam)
    denominator = cos_lat0 - eta * sin_lat0
    lat = np.degrees(np.arctan(numerator / denominator))

    return lon, lat

# Example usage
if __name__ == "__main__":
    # Projection center
    center_lon = 0
    center_lat = 90

    # Spherical coordinates (e.g., in the same system as center)
    lon = [120.1, 121.0]
    lat = [45.0, 46.0]

    xi, eta = gnomonic_projection(center_lon, center_lat, lon, lat)
    print("Projected plane coordinates:", xi, eta)

    lon_rec, lat_rec = inverse_gnomonic_projection(center_lon, center_lat, xi, eta)
    print("Recovered spherical coordinates:", lon_rec, lat_rec)
