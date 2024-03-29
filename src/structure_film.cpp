#include "structure_film.h"

StructureFilm::StructureFilm(int seed, int nTube, int nSegment, double rTube, double lTube, double mTube, double skin, XYZ boxSize, XYZ mean, XYZ std):
  Structure(seed, 0, nSegment, rTube, lTube, mTube, skin, boxSize, mean, std),
  stiffnessBond(0.5 * (86.64 + 100.56*rTube) / lSegment),
  stiffnessHarmonic(0.5 * 63.80 * pow(rTube, 2.93) / lSegment),
  stiffnessBuckling(0.7 * 63.80 * pow(rTube, 2.93) / 275.0),
  bucklingAngle(45.0 * atan(lSegment/275.0) / atan(1))
{
  this->nTube = nTube;
  generate(seed);
  segment(nSegment);
}

StructureFilm::StructureFilm(int seed, int nTube, int nSegment, double rTube, double lTube, double mTube, double skin, XYZ boxSize, double xySigma, double zRange):
  Structure(seed, 0, nSegment, rTube, lTube, mTube, skin, boxSize, XYZ(), XYZ(xySigma, xySigma, zRange)),
  stiffnessBond(0.5 * (86.64 + 100.56*rTube) / lSegment),
  stiffnessHarmonic(0.5 * 63.80 * pow(rTube, 2.93) / lSegment),
  stiffnessBuckling(0.7 * 63.80 * pow(rTube, 2.93) / 275.0),
  bucklingAngle(45.0 * atan(lSegment/275.0) / atan(1))
{
  this->nTube = nTube;
  generateGaussian(seed);
  segment(nSegment);
}

void StructureFilm::printDataFile(string fileName) {
  ofstream file(fileName);

  // find min and max z values

  double zMin = 0.0;
  double zMax = 0.0;
  for (int i = 0; i < atoms.size(); i++) {
    if (atoms[i].z < zMin) zMin = atoms[i].z;
    if (atoms[i].z > zMax) zMax = atoms[i].z;
  }

  // write file header
  
  file << "LAMMPS CNT Data File\n\n";
  
  file << "\t" << atoms.size() << " atoms\n";
  file << "\t" << nSegment*nTube << " bonds\n";
  file << "\t" << (nSegment-1)*nTube << " angles\n";
  file << "\t0 dihedrals\n";
  file << "\t0 impropers\n\n";
  
  file << "\t2 atom types\n";
  file << "\t1 bond types\n";
  file << "\t1 angle types\n";
  file << "\t0 " << boxSize.x << " xlo xhi\n";
  file << "\t0 " << boxSize.y << " ylo yhi\n";
  file << "\t" << zMin - 1.0 << " " << zMax + 1.0 << " zlo zhi\n\n";

  file << "Masses\n\n";

  file << "\t1 " << 0.5 * mAtom << "\n";
  file << "\t2 " << mAtom << "\n\n";

  // write atom data

  file << "Atoms\n\n";

  for (int i = 0; i < atoms.size(); i++) {
    file << i+1 << " " << i / (nSegment + 1) + 1;
    if (i % (nSegment+1) == 0 || i % (nSegment+1) == nSegment)
      file << " 1 0 ";
    else
      file << " 2 0 ";
    file << atoms[i].x << " " << atoms[i].y << " " << atoms[i].z << "\n";
  }

  // write bond data

  file << "\nBonds\n\n";

  int index = 1;
  for (int i = 0; i < atoms.size(); i++) {
    if (i % (nSegment+1) == 0) continue;
    file << index++ << " 1 " << i << " " << i+1 << "\n"; 
  }

  file << "\nAngles\n\n";

  index = 1;
  for (int i = 0; i < atoms.size(); i++) {
    if (i % (nSegment+1) == 0 || i % (nSegment+1) == nSegment) continue;
    file << index++ << " 1 " << i << " " << i+1 << " " << i+2 << "\n";
  }

  file.close();
}

void StructureFilm::printInputFile(int steps, double temp, string fileName) {
  ofstream file(fileName);

  string tab = "\t\t\t";

  // general simulation setup

  file << "#Initialisation\n\n";

  file << "units metal\n";
  file << "dimension 3\n";
  file << "boundary p p s\n";
  file << "atom_style full\n";
  file << "neighbor " << 4.0*lSegment - 23.77 << " bin\n";
  file << "neigh_modify every 5 delay 0 check yes\n";
  file << "newton on\n\n";
  
  file << "#Read data\n\n";
  
  file << "read_data cnt.data\n\n";

  file << "#Force field\n\n";

  file << "bond_style harmonic\n";
  file << "bond_coeff 1 " << stiffnessBond << " " << lSegment << "\n";
  file << "angle_style mesocnt\n";
  file << "angle_coeff 1 " << stiffnessHarmonic << " " << stiffnessBuckling << " " << bucklingAngle << "\n";
  file << "pair_style mesocnt\n";
  file << "pair_coeff * * C_10_10.mesocnt 1\n\n";

  file << "#Output\n\n";

  file << "thermo 100\n";
  // file << "compute ebond all pe bond\n";
  file << "compute epair all pe/atom pair\n";
  file << "dump custom all custom 1000 cnt.lmp mol x y z c_epair\n";
  file << "dump_modify custom sort id\n\n";

  file << "#Simulation setup\n\n";

  file << "velocity all create 600.0 2022\n";
  file << "timestep 1.0e-2\n";
  file << "fix nvt all nvt temp " << temp << " " << temp << " 100\n";
  file << "run " << steps << "\n";

  file.close();
}

void StructureFilm::generate(int seed) {

  mt19937_64 generator(seed);
  generator.discard(10000);
  double pi = 4 * atan(1);
  double deg2rad = pi / 180.0;

  uniform_real_distribution<double> sDistribution(0.0, 1.0);
  normal_distribution<double> tDistributionX(mean.x, std.x);
  normal_distribution<double> tDistributionY(mean.y, std.y);
  uniform_real_distribution<double> angleDistribution(mean.z - std.z, mean.z + std.z);

  int attempts = 0;
  int tubeIndex = 1;

  while (tubes.size() < nTube) {
    
    // generate new tube

    double sX = sDistribution(generator) * boxSize.x;
    double sY = sDistribution(generator) * boxSize.y;
    double sZ = sDistribution(generator) * boxSize.z;
  
    double tX = tDistributionX(generator);
    double tY = tDistributionY(generator);

    double angleRadians = angleDistribution(generator) * deg2rad;
    double tZ = tan(angleRadians) * sqrt(tX*tX + tY*tY);
 
    XYZ s(sX, sY, sZ);
    XYZ t(tX, tY, tZ);

    t.normalise();
    
    Tube tube(rTube, lTube, s, t);
    
    // check periodic boundary conditions and create image tubes if necessary

    XYZ e = s + lTube*t;
    XYZ sTemp = s;

    // image flags

    vector<int> ixs(1, 0);
    vector<int> iys(1, 0);

    if (e.x < rTube+skin) ixs = {0, 1};
    else if (e.x > boxSize.x-rTube-skin) ixs = {0, -1};
    if (e.y < rTube+skin) iys = {0, 1};
    else if (e.y > boxSize.y-rTube-skin) iys = {0, -1};

    // generate ghost tubes

    vector<Tube> ghosts;
    for (auto ix: ixs) {
      sTemp.x = s.x + ix*boxSize.x;
      for (auto iy: iys) {
        sTemp.y = s.y + iy*boxSize.y;
        ghosts.push_back(Tube(rTube, lTube, sTemp, t));
      }  
    }

    attempts++;

    // add new tube if no overlap is detected

    bool overlap = false;
    for (auto ghost: ghosts) {
      overlap = checkOverlap(ghost);
      if (overlap) break;
    }

    if (!overlap) {
      cout << "Tube " << tubeIndex << "/" << nTube << " generated after " << attempts << " attempts." << endl;
      attempts = 0;
      tubes.push_back(tube);
      for (auto ghost: ghosts)
        ghostTubes.push_back(ghost);
      tubeIndex++;
    }
  }
}

void StructureFilm::generateGaussian(int seed) {

  mt19937_64 generator(seed);
  generator.discard(10000);
  double pi = 4 * atan(1);
  double deg2rad = pi / 180.0;

  uniform_real_distribution<double> sDistribution(0.0, 1.0);
  normal_distribution<double> tDistributionXY(mean.x, std.x);
  discrete_distribution<int> signDistribution(-1, 1);
  uniform_real_distribution<double> angleDistribution(mean.z - std.z, mean.z + std.z);

  int attempts = 0;
  int tubeIndex = 1;

  while (tubes.size() < nTube) {
    
    // generate new tube

    double sX = sDistribution(generator) * boxSize.x;
    double sY = sDistribution(generator) * boxSize.y;
    double sZ = sDistribution(generator) * boxSize.z;

    double xyAngle;
    do {
      xyAngle = tDistributionXY(generator);
    } while (xyAngle < 0.0 || xyAngle > 0.5 * pi);

    int signX = signDistribution(generator);
    int signY = signDistribution(generator);

    double tX = signX * cos(xyAngle);
    double tY = signY * sin(xyAngle);

    double angleRadians = angleDistribution(generator) * deg2rad;
    double tZ = tan(angleRadians) * sqrt(tX*tX + tY*tY);
 
    XYZ s(sX, sY, sZ);
    XYZ t(tX, tY, tZ);

    t.normalise();
    
    Tube tube(rTube, lTube, s, t);
    
    // check periodic boundary conditions and create image tubes if necessary

    XYZ e = s + lTube*t;
    XYZ sTemp = s;

    // image flags

    vector<int> ixs(1, 0);
    vector<int> iys(1, 0);

    if (e.x < rTube+skin) ixs = {0, 1};
    else if (e.x > boxSize.x-rTube-skin) ixs = {0, -1};
    if (e.y < rTube+skin) iys = {0, 1};
    else if (e.y > boxSize.y-rTube-skin) iys = {0, -1};

    // generate ghost tubes

    vector<Tube> ghosts;
    for (auto ix: ixs) {
      sTemp.x = s.x + ix*boxSize.x;
      for (auto iy: iys) {
        sTemp.y = s.y + iy*boxSize.y;
        ghosts.push_back(Tube(rTube, lTube, sTemp, t));
      }  
    }

    attempts++;

    // add new tube if no overlap is detected

    bool overlap = false;
    for (auto ghost: ghosts) {
      overlap = checkOverlap(ghost);
      if (overlap) break;
    }

    if (!overlap) {
      cout << "Tube " << tubeIndex << "/" << nTube << " generated after " << attempts << " attempts." << endl;
      attempts = 0;
      tubes.push_back(tube);
      for (auto ghost: ghosts)
        ghostTubes.push_back(ghost);
      tubeIndex++;
    }
  }
}
