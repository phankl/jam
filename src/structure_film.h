#ifndef STRUCTURE_FILM_HEADER
#define STRUCTURE_FILM_HEADER

#include "structure.h"

using namespace std;

class StructureFilm : public Structure {

  public:

    StructureFilm(int, int, int, double, double, double, double, XYZ, XYZ, XYZ);
    StructureFilm(int, int, int, double, double, double, double, XYZ, double, double);

    void printDataFile(string) override;
    void printInputFile(int, double, string) override;

  protected:

    double stiffnessBond;

    double stiffnessHarmonic;
    double stiffnessBuckling;
    double bucklingAngle;

    void generate(int) override;
    void generateGaussian(int);
};

#endif
