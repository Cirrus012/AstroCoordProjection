# gnomonic_projection.py
"""
Gnomonic projection and its inverse.
Numerically stable, handles edge cases such as projection near poles and opposite points.
All angles (longitude/latitude, projection center) should be in degrees.
Projection center and points to be projected must be in the same spherical coordinate system.
"""

import numpy as np
import warnings

def gnomonic_projection(center_lon, center_lat, lon, lat, epsilon=1e-12):
    """
    Project spherical coordinates onto a tangent plane using gnomonic projection.
    Numerically safe for all spherical coordinates; returns nan for undefined cases (e.g., points at/near anti-center).

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
    epsilon : float
        Small value to guard against division by zero (default: 1e-12)

    Returns
    -------
    xi : ndarray
        X coordinates on tangent plane (nan where projection is undefined)
    eta : ndarray
        Y coordinates on tangent plane (nan where projection is undefined)
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)

    dlon = np.radians((lon - center_lon + 180) % 360 - 180)
    sin_lat0 = np.sin(np.radians(center_lat))
    cos_lat0 = np.cos(np.radians(center_lat))
    sin_lat = np.sin(np.radians(lat))
    cos_lat = np.cos(np.radians(lat))

    denom = sin_lat * sin_lat0 + cos_lat * cos_lat0 * np.cos(dlon)

    # Protect against projection singularities (denom ~ 0)
    mask = np.abs(denom) < epsilon
    if np.any(mask):
        warnings.warn("Some points cannot be projected (denominator near zero). Output set to nan.", RuntimeWarning)

    # Standard projection
    xi = np.full_like(denom, np.nan)
    eta = np.full_like(denom, np.nan)
    # Well-defined cases
    good = ~mask

    xi[good] = cos_lat[good] * np.sin(dlon[good]) / denom[good]
    eta[good] = (sin_lat[good] * cos_lat0 - cos_lat[good] * sin_lat0 * np.cos(dlon[good])) / denom[good]

    # Points exactly at projection center (should output 0, 0)
    same_point = (np.abs(lon - center_lon) < epsilon) & (np.abs(lat - center_lat) < epsilon)
    xi[same_point] = 0.0
    eta[same_point] = 0.0

    return xi, eta

def inverse_gnomonic_projection(center_lon, center_lat, xi, eta, epsilon=1e-12):
    """
    Inverse gnomonic projection: recover spherical coordinates from plane coordinates.
    Numerically stable for all input; returns nan for undefined cases (e.g., projection of anti-center).

    Parameters
    ----------
    center_lon : float
        Longitude of projection center (degrees)
    center_lat : float
        Latitude of projection center (degrees)
    xi : array-like
        X coordinates on tangent plane
    eta : array-like
        Y coordinates on tangent plane
    epsilon : float
        Small value to guard against division by zero (default: 1e-12)

    Returns
    -------
    lon : ndarray
        Longitudes of recovered points (degrees, nan if undefined)
    lat : ndarray
        Latitudes of recovered points (degrees, nan if undefined)
    """
    xi = np.asarray(xi, dtype=float)
    eta = np.asarray(eta, dtype=float)
    lat0_rad = np.radians(center_lat)

    sin_lat0 = np.sin(lat0_rad)
    cos_lat0 = np.cos(lat0_rad)

    denom = cos_lat0 - eta * sin_lat0

    # Protect against division by zero
    mask = np.abs(denom) < epsilon
    if np.any(mask):
        warnings.warn("Some points cannot be inverse-projected (denominator near zero). Output set to nan.", RuntimeWarning)

    lon = np.full_like(denom, np.nan)
    lat = np.full_like(denom, np.nan)
    good = ~mask

    # Longitude (in degrees, normalized to [0, 360))
    lam = np.arctan2(xi[good], denom[good])
    lon[good] = (np.degrees(lam) + center_lon) % 360

    # Latitude (in degrees)
    numerator = (eta[good] * cos_lat0 + sin_lat0) * np.cos(lam)
    denominator = cos_lat0 - eta[good] * sin_lat0
    lat[good] = np.degrees(np.arctan(numerator / denominator))

    # Points at projection center (xi=0, eta=0)
    center_mask = (np.abs(xi) < epsilon) & (np.abs(eta) < epsilon)
    lon[center_mask] = center_lon
    lat[center_mask] = center_lat

    return lon, lat

# Example usage
if __name__ == "__main__":
    center_lon = 0
    center_lat = 90  # North pole
    # Try points at the pole, at equator, and at anti-pole
    lon = [0, 90, 180, 0]
    lat = [90, 0, -90, 89.999999]
    xi, eta = gnomonic_projection(center_lon, center_lat, lon, lat)
    print("Projected plane coordinates:", xi, eta)

    lon_rec, lat_rec = inverse_gnomonic_projection(center_lon, center_lat, xi, eta)
    print("Recovered spherical coordinates:", lon_rec, lat_rec)
