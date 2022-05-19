#ifndef STRUCTURE_HEADER
#define STRUCTURE_HEADER

#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "tube.h"
#include "xyz.h"

using namespace std;

class Structure {
  
  public:
   
    Structure(int, int, int, double, double, double, double, XYZ, XYZ, XYZ);
    
    vector<double> odf(int);

    virtual void printDataFile(string) = 0;
    virtual void printInputFile(int, double, string) = 0;
    
    XYZ boxSize;

    XYZ mean;
    XYZ std;

  protected:
    
    int nTube;
    int nSegment;
    double rTube;
    double lTube;
    double mTube;
    double mAtom;
    double lSegment;

    double skin;

    vector<Tube> tubes;
    vector<Tube> ghostTubes;

    vector<XYZ> atoms;
    
    double distance(Tube, Tube, double, double);
    bool checkOverlap(Tube);

    virtual void generate(int);
    void segment(int);
};

#endif
