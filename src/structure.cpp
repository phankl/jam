#include "structure.h"

Structure::Structure(int seed, int nTubeNew, int nSegmentNew, double rTubeNew, double lTubeNew, double mTubeNew, XYZ boxSizeNew, XYZ meanNew, XYZ stdNew):
  nTube(nTubeNew),
  nSegment(nSegmentNew),
  rTube(rTubeNew),
  lTube(lTubeNew),
  mTube(mTubeNew),
  boxSize(boxSizeNew),
  mean(meanNew),
  std(stdNew)
{
  generate(seed);
  segment(nSegment);
  mAtom = mTube / (nSegment+1);
  calculateInertia();
  calculateCOMs();
}

void Structure::printDataFile(string fileName) {
  ofstream file(fileName);

  // write file header
  
  file << "LAMMPS CNT Data File\n\n";
  
  file << "\t" << atoms.size() << " atoms\n";
  file << "\t" << nSegment*nTube << " bonds\n";
  file << "\t0 angles\n";
  file << "\t0 dihedrals\n";
  file << "\t0 impropers\n\n";
  
  file << "\t1 atom types\n";
  file << "\t1 bond types\n";
  file << "\t0 " << boxSize.x << " xlo xhi\n";
  file << "\t0 " << boxSize.y << " ylo yhi\n";
  file << "\t0 " << boxSize.z << " zlo zhi\n\n";

  file << "Masses\n\n";

  file << "\t1 " << mAtom << "\n\n";

  // write atom data

  file << "Atoms\n\n";

  for (int i = 0; i < atoms.size(); i++) {
    file << i+1 << " " << i / (nSegment + 1) + 1;
    file << " 1 0 " << atoms[i].x << " " << atoms[i].y << " " << atoms[i].z << "\n";
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

void Structure::printInertiaFile(string fileName) {
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

void Structure::printInputFile(double neighCutoff, double commCutoff, string fileName) {
  ofstream file(fileName);

  string tab = "\t\t\t";

  // general simulation setup

  file << "#Initialisation\n\n";

  file << "units" << tab << "metal\n";
  file << "dimension" << tab << "3\n";
  file << "boundary" << tab << "p p p\n";
  file << "atom_style" << tab << "full\n";
  file << "comm_modify" << tab << "cutoff " << commCutoff << "\n";
  file << "neighbor" << tab << neighCutoff << " bin\n";
  file << "neigh_modify" << tab << "every 5 delay 0 check yes\n";
  file << "newton" << tab << "on\n\n";
  
  file << "#Read data\n\n";
  
  file << "read_data" << tab << "cnt.data\n\n";

  file << "#Force field\n\n";

  file << "bond_style zero\n";
  file << "bond_coeff *\n";
  file << "pair_style" << tab << "mesocnt\n";
  file << "pair_coeff" << tab << "* * C_10_10.mesocnt\n\n";

  file << "#Output\n\n";

  file << "thermo" << tab << "10\n";
  file << "dump" << tab << "custom all custom 100 cnt.lmp mol x y z mol\n";
  file << "dump_modify" << tab << "custom sort id\n";

  file << "#Simulation setup\n\n";

  file << "velocity" << tab << "all create 600.0 2022\n";
  file << "timestep" << tab << "1.0e-2\n";
  file << "fix rigid all rigid/nvt molecule temp 600 600 100 infile cnt.inertia reinit no\n";
  file << "run" << tab << "10000\n";

  file.close();
}

vector<double> Structure::odf(int nBins) {
  
  vector<double> result(nBins, 0.0);
  double binSize = numbers::pi / nBins;

  for (int i = 0; i < nTube; i++) {
    XYZ t = tubes[i].t;
    double angle = acos(t.z / t.length());
    int bin = floor(angle / binSize);
    result[bin] += 1.0;
  }

  double sum = 0.0;
  for (int i = 0; i < nBins; i++)
    sum += result[i] * binSize;

  for (int i = 0; i < nBins; i++)
    result[i] /= sin((i+0.5) * binSize) * sum;

  return result;
}

double Structure::distance(Tube tube1, Tube tube2, double lambda1, double lambda2) {
  return (tube1.s - tube2.s + lambda1*tube1.t - lambda2*tube2.t).length();
}

bool Structure::checkOverlap(Tube tube1) {
  
  double r1 = tube1.r;
  double l1 = tube1.l;

  XYZ t1 = tube1.t;
  XYZ s1 = tube1.s;

  for (int i = 0; i < ghostTubes.size(); i++) {   
    
    Tube tube2 = ghostTubes[i];
    
    double r2 = tube2.r;
    double l2 = tube2.l;

    XYZ t2 = tube2.t;
    XYZ s2 = tube2.s;

    // check for minimum distance of lines

    XYZ n = cross(t1, t2);
    n.normalise();

    XYZ delta = s2 - s1;

    double d = fabs(delta * n);

    if (d > r1+r2) continue;

    // check if minimum distance is within line segments

    double c = t1 * t2;
    double frac = 1.0 / (1.0 - c*c);
    
    double lambda1 = (t1 - c*t2) * delta * frac;
    double lambda2 = (c*t1 - t2) * delta * frac;
 
    if (lambda1 >= 0.0 && lambda1 <= l1 && lambda2 >= 0.0 && lambda2 <= l2) {
      if (d < r1+r2) return true;
      else continue;
    }

    // check cases where minimum distance is outside the line segments
    // not fully exploiting convexity on subdomains yet, kinda brute force

    // domain edges
    
    lambda1 = 0.0;
    lambda2 = -(delta * t2);
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda2 >= 0.0 && lambda2 <= l2 && d < r1+r2) return true;

    lambda1 = l1;
    lambda2 = c*l1 - delta*t2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda2 >= 0.0 && lambda2 <= l2 && d < r1+r2) return true;

    lambda1 = delta * t1;
    lambda2 = 0.0;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda1 >= 0.0 && lambda1 <= l1 && d < r1+r2) return true;

    lambda1 = c*l2 + delta*t1;
    lambda2 = l2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda1 >= 0.0 && lambda1 <= l1 && d < r1+r2) return true;

    // domain vertices

    lambda1 = 0.0;
    lambda2 = 0.0;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2) return true;

    lambda1 = 0.0;
    lambda2 = l2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2) return true;

    lambda1 = l1;
    lambda2 = 0.0;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2) return true;
    
    lambda1 = l1;
    lambda2 = l2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2) return true;
  }

  return false;
}

void Structure::generate(int seed) {

  mt19937_64 generator(seed);
  generator.discard(10000);

  uniform_real_distribution<double> sDistribution(0.0, 1.0);
  normal_distribution<double> tDistributionX(mean.x, std.x);
  normal_distribution<double> tDistributionY(mean.y, std.y);
  normal_distribution<double> tDistributionZ(mean.z, std.z);

  int attempts = 0;
  int tubeIndex = 1;

  while (tubes.size() < nTube) {
    
    // generate new tube

    double sX = sDistribution(generator) * boxSize.x;
    double sY = sDistribution(generator) * boxSize.y;
    double sZ = sDistribution(generator) * boxSize.z;
  
    double tX = tDistributionX(generator);
    double tY = tDistributionY(generator);
    double tZ = tDistributionZ(generator);
 
    if (sDistribution(generator) < 0.5) tZ *= -1.0;
    
    XYZ s(sX, sY, sZ);
    XYZ t(tX, tY, tZ);

    t.normalise();
    
    // move centre of mass inside simulation box

    XYZ com = s + 0.5*lTube*t;
    if (com.x > boxSize.x) s.x -= boxSize.x;
    else if (com.x < 0) s.x += boxSize.x;
    if (com.y > boxSize.y) s.y -= boxSize.y;
    else if (com.y < 0) s.y += boxSize.y;
    if (com.z > boxSize.z) s.z -= boxSize.z;
    else if (com.z < 0) s.z += boxSize.z;

    Tube tube(rTube, lTube, s, t);
    
    // check periodic boundary conditions and create secondary ghost tube for trial

    XYZ eTemp = s + lTube*t;
    XYZ sTemp = s;
    
    bool ghost = false;

    if (eTemp.x < 0.0) {
      sTemp.x += boxSize.x;
      ghost = true;
    }
    if (eTemp.y < 0.0) {
      sTemp.y += boxSize.y; 
      ghost = true;
    }
    if (eTemp.z < 0.0) {
      sTemp.z += boxSize.z; 
      ghost = true;
    }
    if (eTemp.x > boxSize.x) {
      sTemp.x -= boxSize.x;
      ghost = true;
    }
    if (eTemp.y > boxSize.y) {
      sTemp.y -= boxSize.y;
      ghost = true;
    }
    if (eTemp.z > boxSize.z) {
      sTemp.z -= boxSize.z;
      ghost = true;
    }

    attempts++;

    // add new tube if no overlap is detected

    bool ghostOverlap = (ghost) ? checkOverlap(Tube(rTube, lTube, sTemp, t)) : false;

    if (!checkOverlap(tube) && !ghostOverlap) {
      
      cout << "Tube " << tubeIndex << "/" << nTube << " generated after " << attempts << " attempts." << endl;
      attempts = 0;
      tubes.push_back(tube);
      ghostTubes.push_back(tube);
      tubeIndex++;

      // check periodic boundary conditions and add ghost tube if boundary is crossed
      
      XYZ eTemp = s + lTube*t;
      XYZ sTemp = s;
      
      bool ghost = false;

      if (eTemp.x < 0.0) {
        sTemp.x += boxSize.x;
        ghost = true;
      }
      if (eTemp.y < 0.0) {
        sTemp.y += boxSize.y; 
        ghost = true;
      }
      if (eTemp.z < 0.0) {
        sTemp.z += boxSize.z; 
        ghost = true;
      }
      if (eTemp.x > boxSize.x) {
        sTemp.x -= boxSize.x;
        ghost = true;
      }
      if (eTemp.y > boxSize.y) {
        sTemp.y -= boxSize.y;
        ghost = true;
      }
      if (eTemp.z > boxSize.z) {
        sTemp.z -= boxSize.z;
        ghost = true;
      }

      if (ghost) ghostTubes.push_back(Tube(rTube, lTube, sTemp, t));
    }
  }
}

void Structure::segment(int nSegment) {
  atoms = vector<XYZ>((nSegment+1)*nTube, XYZ());
  int index = 0;
  
  for (int i = 0; i < nTube; i++) {
    Tube tube = tubes[i];
    double lSegment = tube.l / nSegment;
    XYZ delta = lSegment * tube.t;
    
    atoms[index++] = tube.s;
    for (int j = 0; j < nSegment; j++)
      atoms[index++] = tube.s + (j+1)*delta;
  }
}

void Structure::calculateInertia() {

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

void Structure::calculateCOMs() {
  coms = vector<XYZ>(nTube, XYZ());

  for (int i = 0; i < nTube; i++)
    coms[i] = tubes[i].s + 0.5*lTube*tubes[i].t;
}
