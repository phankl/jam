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
   
    Structure(int, int, int, double, double, double, XYZ, XYZ, XYZ);
    
    vector<double> odf(int);

    void printDataFile(string);
    void printInertiaFile(string);
    void printInputFile(double, double, string);

  private:

    int nTube;
    int nSegment;
    double rTube;
    double lTube;
    double mTube;
    double mAtom;

    XYZ boxSize;

    XYZ mean;
    XYZ std;

    vector<Tube> tubes;
    vector<Tube> ghostTubes;

    vector<XYZ> atoms;
    vector<vector<vector<double>>> inertiaTensors;
    vector<XYZ> coms;
    
    double distance(Tube, Tube, double, double);
    bool checkOverlap(Tube);

    void generate(int);
    void segment(int);
    void calculateInertia();
    void calculateCOMs();
};

#endif
