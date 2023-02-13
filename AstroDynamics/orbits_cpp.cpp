"""

@author: dgodinez
"""

#include <cmath>
#include <iostream>
#include <array>

using namespace std;

array<double, 2> calc_acceleration(array<double, 2> X, array<double, 2> x, double M, double m) {
  array<double, 2> a, A;
  
  double distance = sqrt(pow(x[0] - X[0], 2) + pow(x[1] - X[1], 2));
  
  a[0] = -M * (x[0] - X[0]) / pow(distance, 3);
  a[1] = -M * (x[1] - X[1]) / pow(distance, 3);
  
  A[0] = -m * (X[0] - x[0]) / pow(distance, 3);
  A[1] = -m * (X[1] - x[1]) / pow(distance, 3);

  return {a, A};
}

double calc_energy(array<double, 2> X, array<double, 2> x, double Vx, double vx, double Vy, double vy, double M, double m) {
  double r = sqrt(pow(x[0] - X[0], 2) + pow(x[1] - X[1], 2));
  double Vtot = sqrt(pow(Vx, 2) + pow(Vy, 2));
  double vtot = sqrt(pow(vx, 2) + pow(vy, 2));
  
  double energy = m * pow(vtot, 2) / 2.0 + M * pow(Vtot, 2) / 2.0 - m * M / r;
  
  return energy;
}

double calc_com(double M, double m, array<double, 2> X, array<double, 2> x) {
  double tot_mass = M + m;
  
  double rp = sqrt(pow(x[0], 2) + pow(x[1], 2));
  double rs = sqrt(pow(X[0], 2) + pow(X[1], 2));
  
  double r = (M * rs + m * rp) / tot_mass;
  
  return r;
}

double calc_comv(double M, double m, double Vx, double vx, double Vy, double vy) {
  double tot_mass = M + m;
  
  double vp = sqrt(pow(vx, 2) + pow(vy, 2));
  double vs = sqrt(pow(Vx, 2) + pow(Vy, 2));
  
  double v = (M * vs + m * vp) / tot_mass;
  
  return v;
}

double calc_momentum(double X[2], double x[2], double Vx, double vx, double Vy, double vy) {
double x_ = x[0] - X[0];
double y_ = x[1] - X[1];
double v_x = vx - Vx;
double v_y = vy - Vy;

double r = sqrt(x_ * x_ + y_ * y_);
double phi = atan2(y_, x_);

double vrad = v_x * cos(phi) + v_y * sin(phi);
double vphi = -v_x * sin(phi) + v_y * cos(phi);

double phidot = vphi / r;
double h = r * r * phidot;

return h;
}



