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

#include "tube.h"
#include "xyz.h"

using namespace std;

class Structure {
  
  public:
   
    Structure(int, int, int, double, double, XYZ, XYZ, XYZ);
    
    vector<double> odf(int);

    void printDataFile(double, string);
    void printInputFile(double, double, string);

  private:

    int nTube;
    int nSegment;
    double rTube;
    double lTube;

    XYZ boxSize;

    XYZ mean;
    XYZ std;

    vector<Tube> tubes;
    vector<Tube> ghostTubes;

    vector<XYZ> atoms;
    
    double distance(Tube, Tube, double, double);
    bool checkOverlap(Tube);

    void generate(int);
    void segment(int);
};

#endif
