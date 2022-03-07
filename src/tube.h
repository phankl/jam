#ifndef TUBE_HEADER
#define TUBE_HEADER

#include "xyz.h"

struct Tube{

  public:

    double r;
    double l;

    XYZ s;
    XYZ t;

    Tube(double, double, XYZ, XYZ);
};

#endif
