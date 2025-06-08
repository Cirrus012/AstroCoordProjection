"""projections.py
Combined implementations of gnomonic, azimuthal equidistant and
orthographic projections with their inverse transforms.
All angles are expected in degrees.
"""

import numpy as np

# ----------------------------------------------------------------------
# Gnomonic projection
# ----------------------------------------------------------------------

def gnomonic_projection(center_lon, center_lat, lon, lat):
    """Project lon/lat to tangent plane using gnomonic projection.

    Parameters
    ----------
    center_lon, center_lat : float
        Longitude and latitude of the projection center in degrees.
    lon, lat : array_like
        Coordinates of points to project in degrees.

    Returns
    -------
    x, y : ndarray
        Coordinates on the tangent plane.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)

    dlon = np.radians((lon - center_lon + 180) % 360 - 180)
    sin_lat0 = np.sin(np.radians(center_lat))
    cos_lat0 = np.cos(np.radians(center_lat))
    sin_lat = np.sin(np.radians(lat))
    cos_lat = np.cos(np.radians(lat))

    denom = sin_lat * sin_lat0 + cos_lat * cos_lat0 * np.cos(dlon)
    x = cos_lat * np.sin(dlon) / denom
    y = (sin_lat * cos_lat0 - cos_lat * sin_lat0 * np.cos(dlon)) / denom
    return x, y


def inverse_gnomonic_projection(center_lon, center_lat, x, y):
    """Recover lon/lat from gnomonic projection coordinates."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    lat0_rad = np.radians(center_lat)
    sin_lat0 = np.sin(lat0_rad)
    cos_lat0 = np.cos(lat0_rad)

    denom = cos_lat0 - y * sin_lat0
    lam = np.arctan2(x, denom)
    lon = (np.degrees(lam) + center_lon + 360) % 360
    numerator = (y * cos_lat0 + sin_lat0) * np.cos(lam)
    lat = np.degrees(np.arctan(numerator / denom))
    return lon, lat

# ----------------------------------------------------------------------
# Azimuthal equidistant projection
# ----------------------------------------------------------------------

def azimuthal_equidistant_projection(center_lon, center_lat, lon, lat):
    """Azimuthal equidistant projection."""
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    dlon = np.radians((lon - center_lon + 180) % 360 - 180)
    sin_lat0 = np.sin(np.radians(center_lat))
    cos_lat0 = np.cos(np.radians(center_lat))
    sin_lat = np.sin(np.radians(lat))
    cos_lat = np.cos(np.radians(lat))
    cos_c = sin_lat0 * sin_lat + cos_lat0 * cos_lat * np.cos(dlon)
    cos_c = np.clip(cos_c, -1.0, 1.0)
    c = np.arccos(cos_c)
    k = np.where(c != 0, c / np.sin(c), 1.0)
    x = k * cos_lat * np.sin(dlon)
    y = k * (cos_lat0 * sin_lat - sin_lat0 * cos_lat * np.cos(dlon))
    return x, y


def inverse_azimuthal_equidistant_projection(center_lon, center_lat, x, y):
    """Inverse azimuthal equidistant projection."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rho = np.hypot(x, y)
    c = rho
    sin_c = np.sin(c)
    cos_c = np.cos(c)
    sin_lat0 = np.sin(np.radians(center_lat))
    cos_lat0 = np.cos(np.radians(center_lat))
    lat = np.arcsin(cos_c * sin_lat0 + (y * sin_c * cos_lat0) / np.where(rho != 0, rho, 1))
    lon = np.radians(center_lon) + np.arctan2(x * sin_c, rho * cos_lat0 * cos_c - y * sin_lat0 * sin_c)
    return (np.degrees(lon) + 360) % 360, np.degrees(lat)

# ----------------------------------------------------------------------
# Orthographic projection
# ----------------------------------------------------------------------

def orthographic_projection(center_lon, center_lat, lon, lat):
    """Orthographic projection."""
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    dlon = np.radians((lon - center_lon + 180) % 360 - 180)
    sin_lat0 = np.sin(np.radians(center_lat))
    cos_lat0 = np.cos(np.radians(center_lat))
    sin_lat = np.sin(np.radians(lat))
    cos_lat = np.cos(np.radians(lat))
    x = cos_lat * np.sin(dlon)
    y = cos_lat0 * sin_lat - sin_lat0 * cos_lat * np.cos(dlon)
    return x, y


def inverse_orthographic_projection(center_lon, center_lat, x, y):
    """Inverse orthographic projection."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rho = np.hypot(x, y)
    if np.any(rho > 1):
        # Points outside the visible hemisphere cannot be projected back
        return np.full_like(x, np.nan), np.full_like(y, np.nan)
    c = np.arcsin(np.clip(rho, -1.0, 1.0))
    sin_c = np.sin(c)
    cos_c = np.cos(c)
    sin_lat0 = np.sin(np.radians(center_lat))
    cos_lat0 = np.cos(np.radians(center_lat))
    lat = np.arcsin(cos_c * sin_lat0 + (y * sin_c * cos_lat0) / np.where(rho != 0, rho, 1))
    lon = np.radians(center_lon) + np.arctan2(x * sin_c, rho * cos_lat0 * cos_c - y * sin_lat0 * sin_c)
    return (np.degrees(lon) + 360) % 360, np.degrees(lat)

