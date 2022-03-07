#include "structure.h"

Structure::Structure(int seed, int nTubeNew, double rTubeNew, double lTubeNew, XYZ boxSizeNew, XYZ meanNew, XYZ stdNew):
  nTube(nTubeNew),
  rTube(rTubeNew),
  lTube(lTubeNew),
  boxSize(boxSizeNew),
  mean(meanNew),
  std(stdNew)
{
  generate(seed);
}

void Structure::segment(int nSegment, string fileName) {
  ofstream file(fileName); 
  
  for (int i = 0; i < nTube; i++) {
    Tube tube = tubes[i];
    XYZ t = tube.t;
    XYZ s = tube.s;
  }
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
      
      // cout << "Tube " << tubeIndex << " generated after " << attempts << " attempts." << endl;
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