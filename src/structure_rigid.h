#ifndef STRUCTURE_RIGID_HEADER
#define STRUCTURE_RIGID_HEADER

#include "structure.h"

using namespace std;

class StructureRigid : public Structure {

  public:

    StructureRigid(int, int, int, double, double, double, double, XYZ, XYZ, XYZ);

    void printDataFile(string) override;
    void printInputFile(int, double, string) override;
    void printInertiaFile(string);

  protected:

    vector<vector<vector<double>>> inertiaTensors;
    
    vector<XYZ> coms;
    
    void calculateCOMs();
    void calculateInertia();
};

#endif
