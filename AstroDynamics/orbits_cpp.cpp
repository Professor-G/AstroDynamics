#include <cmath>
#include <iostream>
#include <vector>
#include <tuple>

// Function to calculate the acceleration of both bodies
std::pair<double, double> calc_acceleration(const std::vector<double>& X, const std::vector<double>& x, double M, double m, bool approx) {
  double r = sqrt(pow(x[0] - X[0], 2) + pow(x[1] - X[1], 2));
  double a = -M * (x[0] - X[0]) / pow(r, 3);
  double A = -m * (X[0] - x[0]) / pow(r, 3);

  if (approx) {
    a = -M * (x[0] - X[0]) / pow(r, 2);
    A = -m * (X[0] - x[0]) / pow(r, 2);
  }

  return std::make_pair(a, A);
}

// Function to calculate the energy of the system, should be conserved
double calc_energy(const std::vector<double>& X, const std::vector<double>& x, double Vx, double vx, double Vy, double vy, double M, double m) {
  double r = sqrt(pow(x[0] - X[0], 2) + pow(x[1] - X[1], 2));

  double Vtot = sqrt(Vx * Vx + Vy * Vy);
  double vtot = sqrt(vx * vx + vy * vy);

  return m * vtot * vtot / 2.0 + M * Vtot * Vtot / 2.0 - m * M / r;
}

// Function to calculate the C.O.M. of the system
double calc_com(double M, double m, const std::vector<double>& X, const std::vector<double>& x) {
  double tot_mass = M + m;

  double rp = sqrt(x[0] * x[0] + x[1] * x[1]);
  double rs = sqrt(X[0] * X[0] + X[1] * X[1]);

  return (M * rs + m * rp) / tot_mass;
}

// Function to calculate the C.O.M. velocity of the system
double calc_comv(double M, double m, double Vx, double vx, double Vy, double vy) {
  double tot_mass = M + m;

  double vp = sqrt(vx * vx + vy * vy);
  double vs = sqrt(Vx * Vx + Vy * Vy);

  return (M * vs + m * vp) / tot_mass;
}

// Function to calculate the angular momentum of the system
double calc_momentum(const std::vector<double>& X, const std::vector<double>& x, double Vx, double vx, double Vy, double vy) {
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


//Below is the function to run the Leapfrog Integrator
std::vector<std::vector<double> > leapfrog(std::vector<double>& X, std::vector<double>& x, std::vector<double>& V, std::vector<double>& v, double tend, double dt, double M, double m) {

    std::vector<std::vector<double> > result; //The vectors defined below will be stored in this result vector for easier access
    std::vector<double> x_vec, y_vec, X_vec, Y_vec, energy, h, com_x, com_vx, timesteps;

    int num_steps = (int)(tend / dt);
    for (int i = 0; i < num_steps; i++) {
        timesteps.push_back(i * dt);

        std::pair<double, double> acc = calc_acceleration(X, x, M, m, true);

        double KICK_x = V[0] + acc.second * dt / 2.0;
        double KICK_y = V[1] + acc.first * dt / 2.0;
        double kick_x = v[0] + acc.first * dt / 2.0;
        double kick_y = v[1] + acc.second * dt / 2.0;

        x[0] += kick_x * dt;
        x[1] += kick_y * dt;
        X[0] += KICK_x * dt;
        X[1] += KICK_y * dt;

        std::pair<double, double> acc2 = calc_acceleration(X, x, M, m, true);

        v[0] = kick_x + acc2.first * dt / 2.0;
        v[1] = kick_y + acc2.second * dt / 2.0;
        V[0] = KICK_x + acc2.second * dt / 2.0;
        V[1] = KICK_y + acc2.first * dt / 2.0;

        x_vec.push_back(x[0]);
        y_vec.push_back(x[1]);
        X_vec.push_back(X[0]);
        Y_vec.push_back(X[1]);

        double energy1 = calc_energy(X, x, V[0], v[0], V[1], v[1], M, m);
        double h1 = calc_momentum(X, x, V[0], v[0], V[1], v[1]);
        double comx = calc_com(M, m, X, x);
        double comvx = calc_comv(M, m, V[0], v[0], V[1], v[1]);

        energy.push_back(energy1);
        h.push_back(h1);
        com_x.push_back(comx);
        com_vx.push_back(comvx);
    }

    //The vectors will be stored in the result vector
    result.push_back(energy);
    result.push_back(x_vec);
    result.push_back(y_vec);
    result.push_back(X_vec);
    result.push_back(Y_vec);
    result.push_back(com_x);
    result.push_back(com_vx);
    result.push_back(timesteps);

    return result;
}


int main() {
  double M = 1.0; // Mass of the larger body
  double m = 0.01; // Mass of the smaller body
  double dt = 0.01; // Timestep for the simulation
  double tend = 10.0; // End time for the simulation

  std::vector<double> X(2);
  X[0] = 1.0;
  X[1] = 0.0;
  std::vector<double> x(2);
  x[0] = 0.0;
  x[1] = 0.0;
  std::vector<double> V(2);
  V[0] = 0.0;
  V[1] = 0.0;
  std::vector<double> v(2);
  v[0] = 0.0;
  v[1] = 1.0;

  std::vector<std::vector<double> > result = leapfrog(X, x, V, v, tend, dt, M, m); // Call the leapfrog function to simulate the system

  // Print out the metrics
  std::cout << "Energy, x_vec, y_vec, X_vec, Y_vec, com_x, com_vx, timestep" << std::endl;
  for (int i = 0; i < result[2].size(); i++) {
    std::cout << result[0][i] << ", " << result[1][i] << ", " << result[2][i] << result[3][i] << ", " << result[4][i] << ", " << result[5][i] << result[6][i] << ", " << result[7][i] << std::endl;
  }

  return 0;
}
