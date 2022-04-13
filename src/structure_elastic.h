#ifndef STRUCTURE_ELASTIC_HEADER
#define STRUCTURE_ELASTIC_HEADER

#include "structure.h"

using namespace std;

class StructureElastic : public Structure {

  public:

    StructureElastic(int, int, int, double, double, double, XYZ, XYZ, XYZ);

    void printDataFile(string) override;
    void printInputFile(int, double, string) override;

  protected:

    double stiffnessBond;

    double stiffnessHarmonic;
    double stiffnessBuckling;
    double bucklingAngle;
};

#endif
