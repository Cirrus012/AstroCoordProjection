import numpy as np

# ---------------------------------------------------
# Gnomonic projection (eq.6–9)
# ---------------------------------------------------

def gnomonic_projection(A, D, α, δ):
    """
    Gnomonic projection, Dick 1991 eq.6–7
    """
    α = np.asarray(α, dtype=float)
    δ = np.asarray(δ, dtype=float)
    Δα = np.radians((α - A + 180) % 360 - 180)
    D_rad = np.radians(D)
    δ_rad = np.radians(δ)

    sinD, cosD = np.sin(D_rad), np.cos(D_rad)
    sinδ, cosδ = np.sin(δ_rad), np.cos(δ_rad)

    cosθ = sinδ * sinD + cosδ * cosD * np.cos(Δα)   # eq.1
    ξ_G = cosδ * np.sin(Δα) / cosθ                 # eq.6
    η_G = (sinδ * cosD - cosδ * sinD * np.cos(Δα)) / cosθ  # eq.7
    return ξ_G, η_G

def inverse_gnomonic_projection(A, D, ξ_G, η_G):
    """
    Inverse gnomonic projection, Dick 1991 eq.8–9
    """
    ξ_G = np.asarray(ξ_G, dtype=float)
    η_G = np.asarray(η_G, dtype=float)
    D_rad = np.radians(D)
    sinD, cosD = np.sin(D_rad), np.cos(D_rad)

    Δα = np.arctan2(ξ_G, cosD - η_G * sinD)        # eq.8
    α = (np.degrees(Δα) + A + 360) % 360

    num = η_G * cosD + sinD
    den = np.sqrt(ξ_G**2 + (cosD - η_G * sinD)**2)
    δ = np.degrees(np.arctan(num/den))             # eq.9
    return α, δ

# ---------------------------------------------------
# Azimuthal equidistant projection (eq.11–15)
# ---------------------------------------------------

def azimuthal_equidistant_projection(A, D, α, δ):
    """
    Azimuthal equidistant projection, Dick 1991 eq.11–12
    """
    α = np.asarray(α, dtype=float)
    δ = np.asarray(δ, dtype=float)
    Δα = np.radians((α - A + 180) % 360 - 180)
    D_rad = np.radians(D)
    δ_rad = np.radians(δ)

    sinD, cosD = np.sin(D_rad), np.cos(D_rad)
    sinδ, cosδ = np.sin(δ_rad), np.cos(δ_rad)

    # 球面距离和补偿因子
    cosθ = sinδ * sinD + cosδ * cosD * np.cos(Δα)  
    θ    = np.arccos(np.clip(cosθ, -1.0, 1.0))     # eq.1
    sinθ = np.sin(θ)
    k    = np.where(θ != 0, θ/sinθ, 1.0)

    # 正向公式
    ξ_E = k * cosδ * np.sin(Δα)                    # eq.11
    η_E = k * (sinδ * cosD - cosδ * sinD * np.cos(Δα))  # eq.12
    return ξ_E, η_E

def inverse_azimuthal_equidistant_projection(A, D, ξ_E, η_E):
    """
    Inverse azimuthal equidistant projection, Dick 1991 eq.13–15
    """
    ξ_E = np.asarray(ξ_E, dtype=float)
    η_E = np.asarray(η_E, dtype=float)
    ρ    = np.hypot(ξ_E, η_E)
    θ    = ρ
    sinθ = np.sin(θ)
    cosθ = np.cos(θ)

    D_rad = np.radians(D)
    sinD, cosD = np.sin(D_rad), np.cos(D_rad)

    # eq.14
    sinδ = cosθ * sinD + η_E * sinθ * cosD / np.where(ρ != 0, ρ, 1)
    δ    = np.degrees(np.arcsin(np.clip(sinδ, -1.0, 1.0)))

    # eq.15
    num  = ξ_E * sinθ
    den  = ρ * cosD * cosθ - η_E * sinD * sinθ
    Δα   = np.arctan2(num, den)
    α    = (np.degrees(Δα) + A + 360) % 360
    return α, δ

# ---------------------------------------------------
# Orthographic projection (eq.19–23)
# ---------------------------------------------------

def orthographic_projection(A, D, α, δ):
    """
    Orthographic projection, Dick 1991 eq.19–20
    """
    α = np.asarray(α, dtype=float)
    δ = np.asarray(δ, dtype=float)
    Δα   = np.radians((α - A + 180) % 360 - 180)
    D_rad = np.radians(D)
    δ_rad = np.radians(δ)

    sinD, cosD = np.sin(D_rad), np.cos(D_rad)
    sinδ, cosδ = np.sin(δ_rad), np.cos(δ_rad)

    ξ_O = cosδ * np.sin(Δα)                        # eq.19
    η_O = sinδ * cosD - cosδ * sinD * np.cos(Δα)   # eq.20
    return ξ_O, η_O

def inverse_orthographic_projection(A, D, ξ_O, η_O):
    """
    Inverse orthographic projection, Dick 1991 eq.21–23
    """
    ξ_O = np.asarray(ξ_O, dtype=float)
    η_O = np.asarray(η_O, dtype=float)
    D_rad = np.radians(D)
    sinD, cosD = np.sin(D_rad), np.cos(D_rad)

    sinθ2 = ξ_O**2 + η_O**2
    θ     = np.arcsin(np.clip(np.sqrt(sinθ2), -1.0, 1.0))
    sinθ, cosθ = np.sin(θ), np.cos(θ)

    sinδ = η_O * cosD + cosθ * sinD                # eq.22
    δ    = np.degrees(np.arcsin(np.clip(sinδ, -1.0, 1.0)))

    cosδ = np.cos(np.radians(δ))
    Δα   = np.arcsin(np.clip(ξ_O/np.where(cosδ!=0, cosδ,1), -1.0, 1.0))  # eq.23
    α    = (np.degrees(Δα) + A + 360) % 360
    return α, δ
