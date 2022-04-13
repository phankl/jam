#include "structure_rigid.h"

StructureRigid::StructureRigid(int seed, int nTube, int nSegment, double rTube, double lTube, double mTube, double skin, XYZ boxSize, XYZ mean, XYZ std):
  Structure(seed, nTube, nSegment, rTube, lTube, mTube, skin, boxSize, mean, std)
{
  calculateInertia();
  calculateCOMs();
}


void StructureRigid::printDataFile(string fileName) {
  ofstream file(fileName);

  // write file header
  
  file << "LAMMPS CNT Data File\n\n";
  
  file << "\t" << atoms.size() << " atoms\n";
  file << "\t" << nSegment*nTube << " bonds\n";
  file << "\t0 angles\n";
  file << "\t0 dihedrals\n";
  file << "\t0 impropers\n\n";
  
  file << "\t2 atom types\n";
  file << "\t1 bond types\n";
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
  
  file.close();
}

void StructureRigid::printInputFile(int steps, double temp, string fileName) {
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

  file << "bond_style zero\n";
  file << "bond_coeff *\n";
  file << "pair_style mesocnt\n";
  file << "pair_coeff * * C_10_10.mesocnt 1\n\n";

  file << "#Output\n\n";

  file << "thermo 100\n";
  file << "dump custom all custom 1000 cnt.lmp mol x y z mol\n";
  file << "dump_modify custom sort id\n\n";

  file << "#Simulation setup\n\n";

  file << "velocity all create " << temp << " 2022\n";
  file << "timestep 1.0e-2\n";
  file << "fix rigid all rigid/nvt molecule temp " << temp << " " << temp << "100 infile cnt.inertia reinit no\n";
  file << "run " << steps << "\n";

  file.close();
}

void StructureRigid::printInertiaFile(string fileName) {
  ofstream file(fileName);

  file << nTube << "\n";

  for (int i = 0; i < nTube; i++) {
    // id, mass, com
    file << i+1 << " " << mTube << " " << coms[i].x << " " << coms[i].y << " " << coms[i].z << " ";
    // inertia diagonal
    file << inertiaTensors[i][0][0] << " " << inertiaTensors[i][1][1] << " " << inertiaTensors[i][2][2] << " ";
    // inertia off-diagonal
    file << inertiaTensors[i][0][1] << " " << inertiaTensors[i][0][2] << " " << inertiaTensors[i][1][2] << " ";
    // velocity, angular momemtum 
    file << "0 0 0 0 0 0 ";
    // image
    file << "0 0 0\n";
  }

  file.close();
}

void StructureRigid::calculateInertia() {

  inertiaTensors = vector<vector<vector<double>>>(nTube, vector<vector<double>>(3, vector<double>(3, 0)));

  // define inertia tensor of tube in principal coordinate system

  vector<vector<double>> inertia(3, vector<double>(3, 0.0));
  inertia[0][0] = mTube * (6.0*rTube*rTube + lTube*lTube) / 12.0;
  inertia[1][1] = inertia[0][0];
  inertia[2][2] = mTube * rTube * rTube;

  vector<vector<double>> rotation(3, vector<double>(3, 0.0));

  for (int i = 0; i < nTube; i++) {
    XYZ t = tubes[i].t;

    // angle between t and z axis

    double theta = acos(t.z);
    double c = pow(cos(0.5*theta), 2);
    double s = pow(sin(0.5*theta), 2);

    // axis of rotation

    XYZ u(-t.y, t.x, 0.0);
    u.normalise();

    // construct rotation matrix

    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        if (j == k)
          rotation[j][k] = c + s*(2*u[j]*u[j] - 1);
        else {
          rotation[j][k] = 2 * u[j] * u[k] * s;
          for (int l = 0; l < 3; l++) {
            vector<int> ind({j, k, l});
            if (ind == vector<int>({0, 1, 2}) || ind == vector<int>({2, 0, 1}) || ind == vector<int>({1, 2, 0}))
              rotation[j][k] -= u[l] * sin(theta);
            else if (ind == vector<int>({1, 0, 2}) || ind == vector<int>({2, 1, 0}) || ind == vector<int>({0, 2, 1}))
              rotation[j][k] += u[l] * sin(theta);
          }
        }
      }
 
    // rotate inertia tensor

    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          for (int m = 0; m < 3; m++)
            inertiaTensors[i][j][k] += rotation[j][l] * inertia[l][m] * rotation[k][m];
  }
}

void StructureRigid::calculateCOMs() {
  coms = vector<XYZ>(nTube, XYZ());

  for (int i = 0; i < nTube; i++)
    coms[i] = tubes[i].s + 0.5*lTube*tubes[i].t;
}
