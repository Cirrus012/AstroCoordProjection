module projections
contains
  real*8 function rad(d) result(r)
    real*8, intent(in) :: d
    r = d * acos(-1.0d0) / 180.0d0
  end function rad

  real*8 function deg(r)
    real*8, intent(in) :: r
    deg = r * 180.0d0 / acos(-1.0d0)
  end function deg

  subroutine gnomonic_projection(clon, clat, lon, lat, x, y)
    real*8, intent(in) :: clon, clat, lon, lat
    real*8, intent(out) :: x, y
    real*8 dlon, sinlat0, coslat0, sinlat, coslat, denom
    dlon = rad(lon - clon)
    sinlat0 = sin(rad(clat))
    coslat0 = cos(rad(clat))
    sinlat = sin(rad(lat))
    coslat = cos(rad(lat))
    denom = sinlat*sinlat0 + coslat*coslat0*cos(dlon)
    x = coslat*sin(dlon)/denom
    y = (sinlat*coslat0 - coslat*sinlat0*cos(dlon))/denom
  end subroutine gnomonic_projection

  subroutine inverse_gnomonic_projection(clon, clat, x, y, lon, lat)
    real*8, intent(in) :: clon, clat, x, y
    real*8, intent(out) :: lon, lat
    real*8 lat0, sinlat0, coslat0, denom, lam, num
    lat0 = rad(clat)
    sinlat0 = sin(lat0)
    coslat0 = cos(lat0)
    denom = coslat0 - y*sinlat0
    lam = atan2(x, denom)
    lon = mod(deg(lam)+clon+360.d0,360.d0)
    num = (y*coslat0 + sinlat0)*cos(lam)
    lat = deg(atan(num/denom))
  end subroutine inverse_gnomonic_projection

  subroutine azimuthal_equidistant_projection(clon, clat, lon, lat, x, y)
    real*8, intent(in) :: clon, clat, lon, lat
    real*8, intent(out) :: x, y
    real*8 dlon, sinlat0, coslat0, sinlat, coslat, cosc, c, k
    dlon = rad(lon - clon)
    sinlat0 = sin(rad(clat))
    coslat0 = cos(rad(clat))
    sinlat = sin(rad(lat))
    coslat = cos(rad(lat))
    cosc = sinlat0*sinlat + coslat0*coslat*cos(dlon)
    if(cosc > 1.d0) cosc=1.d0
    if(cosc < -1.d0) cosc=-1.d0
    c = acos(cosc)
    if (c /= 0.d0) then
       k = c / sin(c)
    else
       k = 1.d0
    end if
    x = k*coslat*sin(dlon)
    y = k*(coslat0*sinlat - sinlat0*coslat*cos(dlon))
  end subroutine azimuthal_equidistant_projection

  subroutine inverse_azimuthal_equidistant_projection(clon, clat, x, y, lon, lat)
    real*8, intent(in) :: clon, clat, x, y
    real*8, intent(out) :: lon, lat
    real*8 rho, c, sinc, cosc, sinlat0, coslat0, lam
    rho = sqrt(x*x + y*y)
    c = rho
    sinc = sin(c)
    cosc = cos(c)
    sinlat0 = sin(rad(clat))
    coslat0 = cos(rad(clat))
    lat = deg(asin(cosc*sinlat0 + (y*sinc*coslat0)/(merge(rho,1.d0,rho==0.d0))))
    lam = atan2(x*sinc, rho*coslat0*cosc - y*sinlat0*sinc)
    lon = mod(deg(lam)+clon+360.d0,360.d0)
  end subroutine inverse_azimuthal_equidistant_projection

  subroutine orthographic_projection(clon, clat, lon, lat, x, y)
    real*8, intent(in) :: clon, clat, lon, lat
    real*8, intent(out) :: x, y
    real*8 dlon, sinlat0, coslat0, sinlat, coslat
    dlon = rad(lon - clon)
    sinlat0 = sin(rad(clat))
    coslat0 = cos(rad(clat))
    sinlat = sin(rad(lat))
    coslat = cos(rad(lat))
    x = coslat*sin(dlon)
    y = coslat0*sinlat - sinlat0*coslat*cos(dlon)
  end subroutine orthographic_projection

  subroutine inverse_orthographic_projection(clon, clat, x, y, lon, lat)
    real*8, intent(in) :: clon, clat, x, y
    real*8, intent(out) :: lon, lat
    real*8 rho, c, sinc, cosc, sinlat0, coslat0, lam
    rho = sqrt(x*x + y*y)
    if (rho > 1.d0) then
       lon = 0.d0; lat = 0.d0; return
    end if
    c = asin(rho)
    sinc = sin(c)
    cosc = cos(c)
    sinlat0 = sin(rad(clat))
    coslat0 = cos(rad(clat))
    lat = deg(asin(cosc*sinlat0 + (y*sinc*coslat0)/(merge(rho,1.d0,rho==0.d0))))
    lam = atan2(x*sinc, rho*coslat0*cosc - y*sinlat0*sinc)
    lon = mod(deg(lam)+clon+360.d0,360.d0)
  end subroutine inverse_orthographic_projection

end module projections

#ifdef TEST
program test
  use projections
  real*8 x,y,lon,lat
  call gnomonic_projection(0.d0,90.d0,120.d0,45.d0,x,y)
  print *, 'gno', x, y
  call inverse_gnomonic_projection(0.d0,90.d0,x,y,lon,lat)
  print *, 'back', lon, lat
end program test
#endif
