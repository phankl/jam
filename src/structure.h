#ifndef STRUCTURE_HEADER
#define STRUCTURE_HEADER

#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <numbers>

#include <complex>
#include <fftw3.h>

#include "analysis.h"
#include "tube.h"
#include "xyz.h"

using namespace std;

class Structure{
  
  public:
   
    Structure(int, int, double, double, XYZ, XYZ, XYZ);
    
    void segment(int, string);
    vector<double> odf(int);
    
  private:

    int nTube;
    double rTube;
    double lTube;

    XYZ boxSize;

    XYZ mean;
    XYZ std;

    vector<Tube> tubes;
    vector<Tube> ghostTubes;
    
    double distance(Tube, Tube, double, double);
    bool checkOverlap(Tube);

    void generate(int);
};

#endif
