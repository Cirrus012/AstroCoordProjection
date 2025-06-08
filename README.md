# AstroCoordProjection

This repository demonstrates several azimuthal map projections used in astronomy.  Implementations are provided in Python with reference code in C, C++ and Fortran.

## Implemented projections / 已实现投影

- **Gnomonic projection** (`gnomonic_projection`, `inverse_gnomonic_projection`)
- **Azimuthal equidistant projection** (`azimuthal_equidistant_projection`, `inverse_azimuthal_equidistant_projection`)
- **Orthographic projection** (`orthographic_projection`, `inverse_orthographic_projection`)

## Usage / 用法

All angles are in degrees.  The forward functions convert longitude and latitude into planar \(x, y\) coordinates relative to a chosen projection center.  The inverse functions recover spherical coordinates from those plane values.

```
from projections import (
    gnomonic_projection,
    azimuthal_equidistant_projection,
    orthographic_projection,
)

center_lon, center_lat = 0.0, 90.0
lon, lat = [120.0], [45.0]

# Project using the gnomonic projection
xi, eta = gnomonic_projection(center_lon, center_lat, lon, lat)
```

See `show_usage.py` for a complete demonstration that also
exercises the inverse transforms.

### C/C++/Fortran

Simple implementations of the same algorithms are provided in `c/`, `cpp/` and `fortran/` directories for reference.
