#include <math.h>
#include <stdio.h>

/* Convert degrees to radians */
static double rad(double d){ return d * M_PI / 180.0; }
/* Convert radians to degrees */
static double deg(double r){ return r * 180.0 / M_PI; }

/* Difference of longitudes in radians normalized to [-pi, pi] */
static double dlon_rad(double lon, double clon){
    double d = fmod(lon - clon + 180.0, 360.0);
    if(d < 0) d += 360.0;
    d -= 180.0;
    return rad(d);
}

void gnomonic_projection(double clon, double clat, double lon, double lat, double *x, double *y){
    double dlon = dlon_rad(lon, clon);
    double sinlat0 = sin(rad(clat));
    double coslat0 = cos(rad(clat));
    double sinlat = sin(rad(lat));
    double coslat = cos(rad(lat));
    double denom = sinlat * sinlat0 + coslat * coslat0 * cos(dlon);
    *x = coslat * sin(dlon) / denom;
    *y = (sinlat * coslat0 - coslat * sinlat0 * cos(dlon)) / denom;
}

void inverse_gnomonic_projection(double clon, double clat, double x, double y, double *lon, double *lat){
    double lat0 = rad(clat);
    double sinlat0 = sin(lat0);
    double coslat0 = cos(lat0);
    double denom = coslat0 - y * sinlat0;
    double lam = atan2(x, denom);
    *lon = fmod(deg(lam) + clon + 360.0, 360.0);
    double numerator = (y * coslat0 + sinlat0) * cos(lam);
    *lat = deg(atan(numerator / denom));
}

void azimuthal_equidistant_projection(double clon, double clat, double lon, double lat, double *x, double *y){
    double dlon = dlon_rad(lon, clon);
    double sinlat0 = sin(rad(clat));
    double coslat0 = cos(rad(clat));
    double sinlat = sin(rad(lat));
    double coslat = cos(rad(lat));
    double cosc = sinlat0*sinlat + coslat0*coslat*cos(dlon);
    if(cosc>1) cosc=1; if(cosc<-1) cosc=-1;
    double c = acos(cosc);
    double k = (c!=0.0)? c/sin(c):1.0;
    *x = k * coslat * sin(dlon);
    *y = k * (coslat0*sinlat - sinlat0*coslat*cos(dlon));
}

void inverse_azimuthal_equidistant_projection(double clon, double clat, double x, double y, double *lon, double *lat){
    double rho = hypot(x,y);
    double c = rho;
    double sinc = sin(c); double cosc = cos(c);
    double sinlat0 = sin(rad(clat));
    double coslat0 = cos(rad(clat));
    double lat_rad = asin(cosc*sinlat0 + (y*sinc*coslat0)/(rho==0?1:rho));
    double lam = atan2(x*sinc, rho*coslat0*cosc - y*sinlat0*sinc);
    *lon = fmod(deg(lam)+clon+360.0,360.0);
    *lat = deg(lat_rad);
}

void orthographic_projection(double clon, double clat, double lon, double lat, double *x, double *y){
    double dlon = dlon_rad(lon, clon);
    double sinlat0 = sin(rad(clat));
    double coslat0 = cos(rad(clat));
    double sinlat = sin(rad(lat));
    double coslat = cos(rad(lat));
    *x = coslat * sin(dlon);
    *y = coslat0 * sinlat - sinlat0 * coslat * cos(dlon);
}

void inverse_orthographic_projection(double clon, double clat, double x, double y, double *lon, double *lat){
    double rho = hypot(x,y);
    if(rho>1) { *lon=NAN; *lat=NAN; return; }
    double c = asin(rho);
    double sinc = sin(c); double cosc = cos(c);
    double sinlat0 = sin(rad(clat));
    double coslat0 = cos(rad(clat));
    double lat_rad = asin(cosc*sinlat0 + (y*sinc*coslat0)/(rho==0?1:rho));
    double lam = atan2(x*sinc, rho*coslat0*cosc - y*sinlat0*sinc);
    *lon = fmod(deg(lam)+clon+360.0,360.0);
    *lat = deg(lat_rad);
}

#ifdef TEST
int main(){
    double x,y,lon,lat;
    /* Gnomonic */
    gnomonic_projection(0,90,120,45,&x,&y);
    printf("Gnomonic: %f %f\n",x,y);
    inverse_gnomonic_projection(0,90,x,y,&lon,&lat);
    printf("Recovered: %f %f\n",lon,lat);

    /* Azimuthal Equidistant */
    azimuthal_equidistant_projection(0,90,120,45,&x,&y);
    printf("Azimuthal Equidistant: %f %f\n",x,y);
    inverse_azimuthal_equidistant_projection(0,90,x,y,&lon,&lat);
    printf("Recovered: %f %f\n",lon,lat);

    /* Orthographic */
    orthographic_projection(0,90,120,45,&x,&y);
    printf("Orthographic: %f %f\n",x,y);
    inverse_orthographic_projection(0,90,x,y,&lon,&lat);
    printf("Recovered: %f %f\n",lon,lat);
    return 0;
}
#endif
