"""
projections_safe.py
-------------------
Domain‐checking wrappers around projections_core.

status = 0  → 正常
        = 1  → 触发 ε 微调   (near singularity)
        = 2  → 超出有效域但仍返回结果
"""
import numpy as np
import projections_core as pc   # ← 你的原公式脚本

EPS_THETA = 1e-10  # 球面距离角阈值（弧度）
EPS_SHIFT = 1e-12  # 经度/纬度微调量（度）

# ---------- 通用工具 ------------------------------------------------- #
def _spherical_distance(A, D, lon, lat):
    """返回球面距离 θ（弧度），向量化。"""
    lon  = np.asarray(lon,  dtype=float)
    lat  = np.asarray(lat,  dtype=float)
    A, D = np.broadcast_arrays(np.asarray(A, float), np.asarray(D, float))
    dlon = np.radians((lon - A + 180) % 360 - 180)
    sinD, cosD = np.sin(np.radians(D)), np.cos(np.radians(D))
    sinφ, cosφ = np.sin(np.radians(lat)), np.cos(np.radians(lat))
    cosθ = sinφ * sinD + cosφ * cosD * np.cos(dlon)
    return np.arccos(np.clip(cosθ, -1.0, 1.0))

# ---------- Gnomonic ------------------------------------------------ #
def gnomonic_projection(A, D, lon, lat, *, verbose=True):
    θ = _spherical_distance(A, D, lon, lat)
    status = 0
    if np.any(np.isclose(θ, np.pi/2, atol=EPS_THETA)):
        status = 1
        if verbose:
            print("[Gnomonic] 点距中心 θ≈90°，已微调 ±ε。")
        lon = np.asarray(lon, float) + EPS_SHIFT
    elif np.any(θ > np.pi/2):
        status = 2
        if verbose:
            print("[Gnomonic] 有点 θ>90°，理论上发散；仍计算并返回。")
    x, y = pc.gnomonic_projection(A, D, lon, lat)
    return x, y, status

def inverse_gnomonic_projection(A, D, x, y, *, verbose=True):
    x, y = np.asarray(x, float), np.asarray(y, float)
    status = 0
    # 无显式域限；让核心函数处理
    lon, lat = pc.inverse_gnomonic_projection(A, D, x, y)
    return lon, lat, status

# ---------- Azimuthal Equidistant ---------------------------------- #
def azimuthal_equidistant_projection(A, D, lon, lat, *, verbose=True):
    θ = _spherical_distance(A, D, lon, lat)
    status = 0
    if np.any(np.isclose(θ, np.pi, atol=EPS_THETA)):
        status = 1
        if verbose:
            print("[AEQ] 点在对跖位置 θ≈180°，已微调 ±ε。")
        lon = np.asarray(lon, float) + EPS_SHIFT
    elif np.any(θ > np.pi):
        status = 2
        if verbose:
            print("[AEQ] θ>180° 不存在；仍计算并返回（结果可能失真）。")
    x, y = pc.azimuthal_equidistant_projection(A, D, lon, lat)
    return x, y, status

def inverse_azimuthal_equidistant_projection(A, D, x, y, *, verbose=True):
    ρ = np.hypot(x, y)
    status = 0
    if np.any(ρ > np.pi):
        status = 2
        if verbose:
            print("[AEQ‐inv] ρ>π，超出逆变换稳定域；仍计算并返回。")
    lon, lat = pc.inverse_azimuthal_equidistant_projection(A, D, x, y)
    return lon, lat, status

# ---------- Orthographic ------------------------------------------ #
def orthographic_projection(A, D, lon, lat, *, verbose=True):
    θ = _spherical_distance(A, D, lon, lat)
    status = 0
    if np.any(θ > np.pi/2):
        status = 2
        if verbose:
            print("[Ortho] θ>90°，投影点在背半球；结果仅供参考。")
    x, y = pc.orthographic_projection(A, D, lon, lat)
    return x, y, status

def inverse_orthographic_projection(A, D, x, y, *, verbose=True):
    ρ = np.hypot(x, y)
    status = 0
    if np.any(ρ > 1):
        status = 2
        if verbose:
            print("[Ortho‐inv] ρ>1，点在可视半球外；仍返回 NaN。")
    lon, lat = pc.inverse_orthographic_projection(A, D, x, y)
    return lon, lat, status

