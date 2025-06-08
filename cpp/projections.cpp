#include <cmath>
#include <iostream>

static double rad(double d){ return d * M_PI / 180.0; }
static double deg(double r){ return r * 180.0 / M_PI; }

// Difference of longitudes in radians normalized to [-pi, pi]
static double dlon_rad(double lon, double clon){
    double d = fmod(lon - clon + 180.0, 360.0);
    if(d < 0) d += 360.0;
    d -= 180.0;
    return rad(d);
}

struct Point { double x, y; };
struct Spherical { double lon, lat; };

Point gnomonic_projection(double clon, double clat, double lon, double lat){
    double dlon = dlon_rad(lon, clon);
    double sinlat0 = sin(rad(clat));
    double coslat0 = cos(rad(clat));
    double sinlat = sin(rad(lat));
    double coslat = cos(rad(lat));
    double denom = sinlat*sinlat0 + coslat*coslat0*cos(dlon);
    Point p;
    p.x = coslat * sin(dlon) / denom;
    p.y = (sinlat*coslat0 - coslat*sinlat0*cos(dlon)) / denom;
    return p;
}

Spherical inverse_gnomonic_projection(double clon, double clat, double x, double y){
    double lat0 = rad(clat);
    double sinlat0 = sin(lat0);
    double coslat0 = cos(lat0);
    double denom = coslat0 - y*sinlat0;
    double lam = atan2(x, denom);
    Spherical s;
    s.lon = fmod(deg(lam)+clon+360.0,360.0);
    double numerator = (y*coslat0 + sinlat0)*cos(lam);
    s.lat = deg(atan(numerator/denom));
    return s;
}

Point azimuthal_equidistant_projection(double clon, double clat, double lon, double lat){
    double dlon = dlon_rad(lon, clon);
    double sinlat0 = sin(rad(clat));
    double coslat0 = cos(rad(clat));
    double sinlat = sin(rad(lat));
    double coslat = cos(rad(lat));
    double cosc = sinlat0*sinlat + coslat0*coslat*cos(dlon);
    cosc = std::min(1.0, std::max(-1.0, cosc));
    double c = acos(cosc);
    double k = (c!=0.0)? c/sin(c):1.0;
    Point p; p.x = k*coslat*sin(dlon); p.y = k*(coslat0*sinlat - sinlat0*coslat*cos(dlon));
    return p;
}

Spherical inverse_azimuthal_equidistant_projection(double clon, double clat, double x, double y){
    double rho = hypot(x,y); double c = rho;
    double sinc = sin(c), cosc = cos(c);
    double sinlat0 = sin(rad(clat));
    double coslat0 = cos(rad(clat));
    Spherical s;
    s.lat = deg(asin(cosc*sinlat0 + (y*sinc*coslat0)/(rho==0?1:rho)));
    double lam = atan2(x*sinc, rho*coslat0*cosc - y*sinlat0*sinc);
    s.lon = fmod(deg(lam)+clon+360.0,360.0);
    return s;
}

Point orthographic_projection(double clon, double clat, double lon, double lat){
    double dlon = dlon_rad(lon, clon);
    double sinlat0 = sin(rad(clat));
    double coslat0 = cos(rad(clat));
    double sinlat = sin(rad(lat));
    double coslat = cos(rad(lat));
    Point p; p.x = coslat*sin(dlon); p.y = coslat0*sinlat - sinlat0*coslat*cos(dlon); return p; }

Spherical inverse_orthographic_projection(double clon, double clat, double x, double y){
    double rho = hypot(x,y);
    if(rho>1) return {NAN,NAN};
    double c = asin(rho); double sinc=sin(c), cosc=cos(c);
    double sinlat0=sin(rad(clat)), coslat0=cos(rad(clat));
    Spherical s;
    s.lat = deg(asin(cosc*sinlat0 + (y*sinc*coslat0)/(rho==0?1:rho)));
    double lam = atan2(x*sinc, rho*coslat0*cosc - y*sinlat0*sinc);
    s.lon = fmod(deg(lam)+clon+360.0,360.0);
    return s;
}

#ifdef TEST
int main(){
    Point p; Spherical s;

    p = gnomonic_projection(0,90,120,45);
    std::cout << "Gnomonic: "<<p.x<<" "<<p.y<<"\n";
    s = inverse_gnomonic_projection(0,90,p.x,p.y);
    std::cout << "Recovered: "<<s.lon<<" "<<s.lat<<"\n";

    p = azimuthal_equidistant_projection(0,90,120,45);
    std::cout << "Azimuthal Equidistant: "<<p.x<<" "<<p.y<<"\n";
    s = inverse_azimuthal_equidistant_projection(0,90,p.x,p.y);
    std::cout << "Recovered: "<<s.lon<<" "<<s.lat<<"\n";

    p = orthographic_projection(0,90,120,45);
    std::cout << "Orthographic: "<<p.x<<" "<<p.y<<"\n";
    s = inverse_orthographic_projection(0,90,p.x,p.y);
    std::cout << "Recovered: "<<s.lon<<" "<<s.lat<<"\n";
}
#endif
