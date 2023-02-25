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
//Below is the function to run the Runge-Kutta 4 integrator
//////////////////////

std::pair<std::vector<double>, std::vector<double> > RK4_integrator(std::vector<double> X, std::vector<double> x, std::vector<double> V, std::vector<double> v, double tend, double dt, double M, double m) {

  std::vector<double> X_vec, x_vec, Y_vec, y_vec;
  std::vector<double> energy, h, com_x, com_vx;

  for (double t = 0; t < tend; t += dt) {
    std::pair<double, double> acceleration = calc_acceleration(X, x, M, m);
    double a1 = acceleration.first;
    double A1 = acceleration.second;
    std::vector<double> x1(2);
    x1[0] = x[0] + v[0] * dt / 2;
    x1[1] = x[1] + v[1] * dt / 2;
    std::vector<double> X1(2);
    X1[0] = X[0] + V[0] * dt / 2;
    X1[1] = X[1] + V[1] * dt / 2;
    std::pair<double, double> acceleration2 = calc_acceleration(X1, x1, M, m);
    double a2 = acceleration2.first;
    double A2 = acceleration2.second;
    std::vector<double> x2(2);
    x2[0] = x[0] + v[0] * dt / 2;
    x2[1] = x[1] + v[1] * dt / 2;
    std::vector<double> X2(2);
    X2[0] = X[0] + V[0] * dt / 2;
    X2[1] = X[1] + V[1] * dt / 2;
    std::pair<double, double> acceleration3 = calc_acceleration(X2, x2, M, m);
    double a3 = acceleration3.first;
    double A3 = acceleration3.second;
    std::vector<double> x3(2);
    x3[0] = x[0] + v[0] * dt;
    x3[1] = x[1] + v[1] * dt;
    std::vector<double> X3(2);
    X3[0] = X[0] + V[0] * dt;
    X3[1] = X[1] + V[1] * dt;
    std::pair<double, double> acceleration4 = calc_acceleration(X3, x3, M, m);
    double a4 = acceleration4.first;
    double A4 = acceleration4.second;


    // Update positions at each timestamp
    x[0] = x[0] + (v[0] + 2 * (a1 + a2) + a3) * dt / 6;
    x[1] = x[1] + (v[1] + 2 * (A1 + A2) + A3) * dt / 6;
    X[0] = X[0] + (V[0] + 2 * (A1 + A2) + A3) * dt / 6;
    X[1] = X[1] + (V[1] + 2 * (A1 + A2) + A3) * dt / 6;
    
    // Update velocities at each timestamp
    v[0] = v[0] + (a1 + 2 * (a2 + a3) + a4) * dt / 6;
    v[1] = v[1] + (A1 + 2 * (A2 + A3) + A4) * dt / 6;
    V[0] = V[0] + (A1 + 2 * (A2 + A3) + A4) * dt / 6;
    V[1] = V[1] + (a1 + 2 * (a2 + a3) + a4) * dt / 6;

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

  std::pair<std::vector<double>, std::vector<double> > result = RK4_integrator(X, x, V, v, tend, dt, M, m);
  std::vector<double> energy_error = result.first;
  std::vector<double> h_error = result.second;
  
  // Print the resulting energy error and h error
    std::cout << "Time\tTimestep\tX1\tX2\tX3\tV1\tV2\tV3\n";
    for (int i = 0; i < energy_error.size(); i++) {
    double time = i * dt;
        std::cout << time << "\t" << dt << "\t" << X[0] << "\t" << X[1] << "\t" << X[2] << "\t" << V[0] << "\t" << V[1] << "\t" << V[2] << "\n";
  }

return 0;
}
