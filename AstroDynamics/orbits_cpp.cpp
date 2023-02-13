#include <cmath>
#include <iostream>
#include <vector>
#include <ctime>
#include <utility>

//////////////////////
//Helper functions
//////////////////////

//Function to calculate the acceleration of both bodies
std::pair<double, double> calc_acceleration(const std::vector<double> &X, const std::vector<double> &x, double M, double m) {
  double a = -M * (x[0] - X[0]) / pow(pow(x[0] - X[0], 2) + pow(x[1] - X[1], 2), 3.0 / 2.0);
  double A = -m * (X[0] - x[0]) / pow(pow(x[0] - X[0], 2) + pow(x[1] - X[1], 2), 3.0 / 2.0);

  return std::make_pair(a, A);
}

//Function to calculate the energy of the system, should be conserved
double calc_energy(const std::vector<double> &X, const std::vector<double> &x, double Vx, double vx, double Vy, double vy, double M, double m) {
  double r = sqrt(pow(x[0] - X[0], 2) + pow(x[1] - X[1], 2));

  double Vtot = sqrt(Vx * Vx + Vy * Vy);
  double vtot = sqrt(vx * vx + vy * vy);

  return m * vtot * vtot / 2.0 + M * Vtot * Vtot / 2.0 - m * M / r;
}

//Function to calculate the C.O.M. of the system
double calc_com(double M, double m, const std::vector<double> &X, const std::vector<double> &x) {
  double tot_mass = M + m;

  double rp = sqrt(x[0] * x[0] + x[1] * x[1]);
  double rs = sqrt(X[0] * X[0] + X[1] * X[1]);

  return (M * rs + m * rp) / tot_mass;
}

//Function to calculate the C.O.M. velocity of the system
double calc_comv(double M, double m, double Vx, double vx, double Vy, double vy) {
  double tot_mass = M + m;

  double vp = sqrt(vx * vx + vy * vy);
  double vs = sqrt(Vx * Vx + Vy * Vy);

  return (M * vs + m * vp) / tot_mass;
}

//Function to calculate the angular momentum of the system

double calc_momentum(const std::vector<double> &X, const std::vector<double> &x, double Vx, double vx, double Vy, double vy) {
  double x_ = x[0] - X[0];
  double y_ = x[1] - X[1];
  double v_x = vx - Vx;
  double v_y = vy - Vy;

  double r = sqrt(x_ * x_ + y_ * y_);
  double phi = atan2(y_, x_);

  double vrad = v_x * cos(phi) + v_y * sin(phi);
  double vphi = -v_x * sin(phi) + v_y;

  double phidot = vphi / r;
  double h = r * r * phidot;

  return h;
}

//////////////////////
//Below is the function to run the euler integrator
//////////////////////

std::pair<std::vector<double>, std::vector<double> > euler_integrator(std::vector<double> X, std::vector<double> x, std::vector<double> V, std::vector<double> v, double tend, double dt, double M, double m) {
  std::vector<double> X_vec, x_vec, Y_vec, y_vec;
  std::vector<double> energy, h, com_x, com_vx;
  
  for (double t = 0; t < tend; t += dt) {
    auto acceleration = calc_acceleration(X, x, M, m);
    double a = acceleration.first;
    double A = acceleration.second;

    x[0] = x[0] + v[0] * dt;
    x[1] = x[1] + v[1] * dt;
    X[0] = X[0] + V[0] * dt;
    X[1] = X[1] + V[1] * dt;

    v[0] = v[0] + a * dt;
    v[1] = v[1] + a * dt;
    V[0] = V[0] + A * dt;
    V[1] = V[1] + A * dt;

    x_vec.push_back(x[0]);
    y_vec.push_back(x[1]);
    X_vec.push_back(X[0]);
    Y_vec.push_back(X[1]);

    double energy1 = calc_energy(X, x, V[0], v[0], V[1], v[1], M, m);
    double h1 = calc_momentum(X, x, V[0], v[0], V[1], v[1]);
    double comx = calc_com(M, m, X, x);
    double comvx = calc_comv(M, m, V[0], v[0], V[1], v[1]);

    com_x.push_back(comx);
    com_vx.push_back(comvx);
    energy.push_back(energy1);
    h.push_back(h1);
  }

  std::vector<double> energy_error, h_error;
  for (int i = 0; i < energy.size(); i++) {
    energy_error.push_back(std::abs((energy[i] - energy[0]) / energy[0]));
    h_error.push_back(std::abs((h[i] - h[0]) / h[0]));
  }

  return std::make_pair(energy_error, h_error);
}

//Main, which by convention should return an integer value to indicate whether the program ran successfully or not.
int main() {
  std::vector<double> X = {0, 0, 0};
  std::vector<double> x = {1, 0, 0};
  std::vector<double> V = {0, 0, 0};
  std::vector<double> v = {0, 1, 0};
  double M = 1 - 3e-6;
  double m = 3e-6;
  double tend = 10;
  double dt = 1e-2;

  std::pair<std::vector<double>, std::vector<double>> result = euler_integrator(X, x, V, v, tend, dt, M, m);
  std::vector<double> energy = result.first;
  std::vector<double> energy_error = result.second;

  return 0; // return 0 to indicate success
}



