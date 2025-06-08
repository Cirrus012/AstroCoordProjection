"""Example script demonstrating projections usage."""

from projections import (
    gnomonic_projection,
    inverse_gnomonic_projection,
    azimuthal_equidistant_projection,
    inverse_azimuthal_equidistant_projection,
    orthographic_projection,
    inverse_orthographic_projection,
)


def main():
    # Projection center at the north pole
    center_lon, center_lat = 0.0, 90.0

    # Example spherical coordinates to project
    lon = [120.0, 40.0]
    lat = [45.0, 10.0]

    # Gnomonic
    x, y = gnomonic_projection(center_lon, center_lat, lon, lat)
    print("Gnomonic:", x, y)
    lon2, lat2 = inverse_gnomonic_projection(center_lon, center_lat, x, y)
    print("Recovered:", lon2, lat2)

    # Azimuthal equidistant
    x, y = azimuthal_equidistant_projection(center_lon, center_lat, lon, lat)
    print("Azimuthal Equidistant:", x, y)
    lon2, lat2 = inverse_azimuthal_equidistant_projection(center_lon, center_lat, x, y)
    print("Recovered:", lon2, lat2)

    # Orthographic
    x, y = orthographic_projection(center_lon, center_lat, lon, lat)
    print("Orthographic:", x, y)
    lon2, lat2 = inverse_orthographic_projection(center_lon, center_lat, x, y)
    print("Recovered:", lon2, lat2)


if __name__ == "__main__":
    main()
