"""
Azimuthal Equidistant and Orthographic projections.
All angles in degrees.
"""
import numpy as np


def azimuthal_equidistant_projection(center_lon, center_lat, lon, lat):
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    dlon = np.radians((lon - center_lon + 180) % 360 - 180)
    sin_lat0 = np.sin(np.radians(center_lat))
    cos_lat0 = np.cos(np.radians(center_lat))
    sin_lat = np.sin(np.radians(lat))
    cos_lat = np.cos(np.radians(lat))
    cos_c = sin_lat0 * sin_lat + cos_lat0 * cos_lat * np.cos(dlon)
    c = np.arccos(np.clip(cos_c, -1.0, 1.0))
    k = np.where(c != 0, c / np.sin(c), 1.0)
    x = k * cos_lat * np.sin(dlon)
    y = k * (cos_lat0 * sin_lat - sin_lat0 * cos_lat * np.cos(dlon))
    return x, y

def inverse_azimuthal_equidistant_projection(center_lon, center_lat, x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rho = np.sqrt(x**2 + y**2)
    c = rho
    sin_c = np.sin(c)
    cos_c = np.cos(c)
    sin_lat0 = np.sin(np.radians(center_lat))
    cos_lat0 = np.cos(np.radians(center_lat))
    lat = np.arcsin(cos_c * sin_lat0 + (y * sin_c * cos_lat0) / np.where(rho!=0, rho, 1))
    lon = np.radians(center_lon) + np.arctan2(x * sin_c, rho * cos_lat0 * cos_c - y * sin_lat0 * sin_c)
    return (np.degrees(lon) + 360) % 360, np.degrees(lat)

def orthographic_projection(center_lon, center_lat, lon, lat):
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
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rho = np.sqrt(x**2 + y**2)
    c = np.arcsin(np.clip(rho, -1.0, 1.0))
    sin_c = np.sin(c)
    cos_c = np.cos(c)
    sin_lat0 = np.sin(np.radians(center_lat))
    cos_lat0 = np.cos(np.radians(center_lat))
    lat = np.arcsin(cos_c * sin_lat0 + (y * sin_c * cos_lat0) / np.where(rho!=0, rho, 1))
    lon = np.radians(center_lon) + np.arctan2(x * sin_c, rho * cos_lat0 * cos_c - y * sin_lat0 * sin_c)
    return (np.degrees(lon) + 360) % 360, np.degrees(lat)
