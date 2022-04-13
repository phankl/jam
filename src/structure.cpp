#include "structure.h"

Structure::Structure(int seed, int nTube, int nSegment, double rTube, double lTube, double mTube, double skin, XYZ boxSize, XYZ mean, XYZ std): 
  nTube(nTube),
  nSegment(nSegment),
  rTube(rTube),
  lTube(lTube),
  mTube(mTube),
  skin(skin),
  mAtom(mTube / nSegment),
  lSegment(lTube / nSegment),
  boxSize(boxSize),
  mean(mean),
  std(std)
{
  generate(seed);
  segment(nSegment);
}

vector<double> Structure::odf(int nBins) {
  
  vector<double> result(nBins, 0.0);
  double pi = 4 * atan(1);
  double binSize = pi / nBins;

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

    if (d > r1+r2+skin) continue;

    // check if minimum distance is within line segments

    double c = t1 * t2;
    double frac = 1.0 / (1.0 - c*c);
    
    double lambda1 = (t1 - c*t2) * delta * frac;
    double lambda2 = (c*t1 - t2) * delta * frac;
 
    if (lambda1 >= 0.0 && lambda1 <= l1 && lambda2 >= 0.0 && lambda2 <= l2) {
      if (d < r1+r2+skin) return true;
      else continue;
    }

    // check cases where minimum distance is outside the line segments
    // not fully exploiting convexity on subdomains yet, kinda brute force

    // domain edges
    
    lambda1 = 0.0;
    lambda2 = -(delta * t2);
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda2 >= 0.0 && lambda2 <= l2 && d < r1+r2+skin) return true;

    lambda1 = l1;
    lambda2 = c*l1 - delta*t2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda2 >= 0.0 && lambda2 <= l2 && d < r1+r2+skin) return true;

    lambda1 = delta * t1;
    lambda2 = 0.0;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda1 >= 0.0 && lambda1 <= l1 && d < r1+r2+skin) return true;

    lambda1 = c*l2 + delta*t1;
    lambda2 = l2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (lambda1 >= 0.0 && lambda1 <= l1 && d < r1+r2+skin) return true;

    // domain vertices

    lambda1 = 0.0;
    lambda2 = 0.0;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2+skin) return true;

    lambda1 = 0.0;
    lambda2 = l2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2+skin) return true;

    lambda1 = l1;
    lambda2 = 0.0;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2+skin) return true;
    
    lambda1 = l1;
    lambda2 = l2;
    d = distance(tube1, tube2, lambda1, lambda2);
    if (d < r1+r2+skin) return true;
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
 
    // if (sDistribution(generator) < 0.5) tZ *= -1.0;
    
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
    vector<int> izs(1, 0);

    if (e.x < 0.0) ixs = {0, 1};
    else if (e.x > boxSize.x) ixs = {0, -1};
    if (e.y < 0.0) iys = {0, 1};
    else if (e.y > boxSize.y) iys = {0, -1};
    if (e.z < 0.0) izs = {0, 1};
    else if (e.z > boxSize.z) izs = {0, -1};

    // generate ghost tubes

    vector<Tube> ghosts;
    for (auto ix: ixs) {
      sTemp.x = s.x + ix*boxSize.x;
      for (auto iy: iys) {
        sTemp.y = s.y + iy*boxSize.y;
        for (auto iz: izs) {
          sTemp.z = s.z + iz*boxSize.z;
          ghosts.push_back(Tube(rTube, lTube, sTemp, t));
        }
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
