#include <iostream>
#include <vector>
#include <utility>
#include <cmath>


//////////////////////
//Helper functions
//////////////////////

//Function to calculate the acceleration of both bodies
std::pair<double, double> calc_acceleration(std::vector<double> X, std::vector<double> x, double M, double m) {
    double r = std::sqrt((X[0] - x[0]) * (X[0] - x[0]) + (X[1] - x[1]) * (X[1] - x[1]));
    double a = M * (X[0] - x[0]) / (r * r * r) + m * (x[0] - X[0]) / (r * r * r);
    double A = M * (X[1] - x[1]) / (r * r * r) + m * (x[1] - X[1]) / (r * r * r);
    return std::make_pair(a, A);
}

//Function to calculate the energy of the system, should be conserved
double calc_energy(std::vector<double> X, std::vector<double> x, double Vx, double vx, double Vy, double vy, double M, double m) {
    double r = std::sqrt((X[0] - x[0]) * (X[0] - x[0]) + (X[1] - x[1]) * (X[1] - x[1]));
    double T = 0.5 * M * (Vx * Vx + Vy * Vy) + 0.5 * m * (vx * vx + vy * vy);
    double U = -M * m / r;
    return T + U;
}

//Function to calculate the angular momentum of the system
double calc_momentum(std::vector<double> X, std::vector<double> x, double Vx, double vx, double Vy, double vy) {
    return (Vx - vx) * (X[1] - x[1]) - (Vy - vy) * (X[0] - x[0]);
}

//Function to calculate the C.O.M. of the system
double calc_com(double M, double m, std::vector<double> X, std::vector<double> x) {
    return (M * X[0] + m * x[0]) / (M + m);
}

//Function to calculate the C.O.M. velocity of the system
double calc_comv(double M, double m, double Vx, double vx, double Vy, double vy) {
    return (M * Vx + m * vx) / (M + m);
}

//////////////////////
//Below is the function to run the euler integrator
//////////////////////

std::pair<std::vector<double>, std::vector<double> > euler_integrator(std::vector<double> X, std::vector<double> x, std::vector<double> V, std::vector<double> v, double tend, double dt, double M, double m) {
  std::vector<double> X_vec, x_vec, Y_vec, y_vec;
  std::vector<double> energy, h, com_x, com_vx;

  for (double t = 0; t < tend; t += dt) {
    std::pair<double, double> acceleration = calc_acceleration(X, x, M, m);
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
  std::vector<double> X(3);
  X[0] = 0;
  X[1] = 0;
  X[2] = 0;

  std::vector<double> x(3);
  x[0] = 1;
  x[1] = 0;
  x[2] = 0;

  std::vector<double> V(3);
  V[0] = 0;
  V[1] = 0;
  V[2] = 0;

  std::vector<double> v(3);
  v[0] = 0;
  v[1] = 1;
  v[2] = 0;

  double M = 1 - 3e-6;
  double m = 3e-6;
  double tend = 10;
  double dt = 1e-2;

  std::pair<std::vector<double>, std::vector<double> > result = euler_integrator(X, x, V, v, tend, dt, M, m);
  std::vector<double> energy_error = result.first;
  std::vector<double> h_error = result.second;

  // Print the resulting energy error and h error
  std::cout << "Energy error:\n";
  for (double error : energy_error) {
    std::cout << error << "\n";
  }

  std::cout << "h error:\n";
  for (double error : h_error) {
    std::cout << error << "\n";
  }

  return 0;
}







