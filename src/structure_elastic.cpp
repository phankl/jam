#include "structure_elastic.h"

StructureElastic::StructureElastic(int seed, int nTube, int nSegment, double rTube, double lTube, double mTube, double skin, XYZ boxSize, XYZ mean, XYZ std):
  Structure(seed, nTube, nSegment, rTube, lTube, mTube, skin, boxSize, mean, std),
  stiffnessBond(0.5 * (86.64 + 100.56*rTube) / lSegment),
  stiffnessHarmonic(0.5 * 63.80 * pow(rTube, 2.93) / lSegment),
  stiffnessBuckling(0.7 * 63.80 * pow(rTube, 2.93) / 275.0),
  bucklingAngle(180.0 * atan(lSegment/275.0) / atan(1))
{
}

void StructureElastic::printDataFile(string fileName) {
  ofstream file(fileName);

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
  file << "\t0 " << boxSize.z << " zlo zhi\n\n";

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

void StructureElastic::printInputFile(int steps, double temp, string fileName) {
  ofstream file(fileName);

  string tab = "\t\t\t";

  // general simulation setup

  file << "#Initialisation\n\n";

  file << "units metal\n";
  file << "dimension 3\n";
  file << "boundary p p p\n";
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
