#ifndef XYZ_HEADER
#define XYZ_HEADER

#include <cmath>

struct XYZ {
  
  public:
    
    double x;
    double y;
    double z;
  
    XYZ();
    XYZ(double, double, double);
    
    inline double length() {
      return sqrt(x*x + y*y + z*z);
    };
    
    inline void normalise() {
      double length = sqrt(x*x + y*y + z*z);
      x /= length;
      y /= length;
      z /= length;
    };

    double& operator [] (int i) {
      if (i == 0) return x;
      else if (i == 1) return y;
      else return z;
    };
};

inline XYZ operator + (const XYZ& a, const XYZ& b) {
  return XYZ(a.x + b.x,a.y + b.y,a.z + b.z);
}

inline XYZ operator - (const XYZ& a, const XYZ& b) {
  return XYZ(a.x - b.x,a.y - b.y,a.z - b.z);
}

inline XYZ operator * (double a, const XYZ& b) {
  return XYZ(a*b.x,a*b.y,a*b.z);
}

inline double operator * (const XYZ& a, const XYZ& b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline XYZ cross(const XYZ& a, const XYZ& b) {
  return XYZ(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
};

XYZ rotate(const XYZ&, const XYZ&, double);

bool pbc(const XYZ&, const XYZ&, const XYZ&, XYZ&);

#endif
