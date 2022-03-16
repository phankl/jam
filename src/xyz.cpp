#include "xyz.h"

XYZ::XYZ() :
  x(0.0),
  y(0.0),
  z(0.0)
{
}

XYZ::XYZ(double xNew, double yNew, double zNew) :
  x(xNew),
  y(yNew),
  z(zNew)
{
}


XYZ rotate(const XYZ& r, const XYZ& n, double theta) {
  double c = cos(theta);
  double s = sin(theta);
  return (1.0-c)*(n*r)*n + c*r - s*cross(n,r);
}

bool pbc(const XYZ& xyzMin, const XYZ& xyzMax, const XYZ& r1, XYZ& r2) {
  double dx = xyzMax.x - xyzMin.x;
  double dy = xyzMax.y - xyzMin.y;
  double dz = xyzMax.z - xyzMin.z;
  double dxHalf = 0.5 * dx;
  double dyHalf = 0.5 * dy;
  double dzHalf = 0.5 * dz;

  bool wrapped = false;

  if (r2.x - r1.x > dxHalf) {
    r2.x -= dx;
    wrapped = true;
  }
  else if (r2.x - r1.x < -dxHalf) {
    r2.x += dx;
    wrapped = true;
  }
  if (r2.y - r1.y > dyHalf) {
    r2.y -= dy;
    wrapped = true;
  }
  else if (r2.y - r1.y < -dyHalf) {
    r2.y += dy;
    wrapped = true;
  }
  if (r2.z - r1.z > dzHalf) {
    r2.z -= dz;
    wrapped = true;
  }
  else if (r2.z - r1.z < -dzHalf) {
    r2.z += dz;
    wrapped = true;
  }

  return wrapped;
}
