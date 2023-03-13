#include <cmath>
#include <iostream>
#include <vector>
#include <tuple>

// Function to calculate the acceleration of both bodies, including approx option (default is true)
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

//Below is the function to run the Leapfrog Integrator
void leapfrog(std::vector<double>& X, std::vector<double>& x, std::vector<double>& V, std::vector<double>& v, double tend, double dt, double M, double m) {

    std::cout << "time, timestep, x_star, y_star, vx_star, vy_star, x_planet, y_planet, vx_planet, vy_planet" << std::endl;

    int num_steps = (int)(tend / dt);
    for (int i = 0; i < num_steps; i++) {
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

        std::cout << i*dt << ", " << dt << ", " << x[0] << ", " << x[1] << ", " << v[0] << ", " << v[1] << ", " << X[0] << ", " << X[1] << ", " << V[0] << ", " << V[1] << std::endl;
    }
}

int main() {
  double M = 1.0; // Mass of the larger body
  double m = 0.01; // Mass of the smaller body
  double dt = 0.01; // Timestep for the simulation
  double tend = 10.0; // End time for the simulation

  std::vector<double> X(2); // Initial position of the larger body
  X[0] = 1.0;
  X[1] = 0.0;

  std::vector<double> x(2); // Initial position of the smaller body
  x[0] = 0.0;
  x[1] = 0.0;

  std::vector<double> V(2); // Initial velocity of the larger body
  V[0] = 0.0;
  V[1] = 0.0;

  std::vector<double> v(2); // Initial velocity of the smaller body
  v[0] = 0.0;
  v[1] = 1.0;

  // Open output file for writing
  std::ofstream outfile;
  outfile.open("output.txt");
  // outfile.open("/Users/daniel/Desktop/output.txt");

  // Call the leapfrog function to simulate the system and output results to file
  outfile << "time, timestep, x_star, y_star, vx_star, vy_star, x_planet, y_planet, vx_planet, vy_planet" << std::endl;
  leapfrog(X, x, V, v, tend, dt, M, m, outfile); 

  // Close output file
  // outfile.close();

  return 0;
}

// PYTHON SCRIPT TO PLOT THE ENERGIES
/*
import numpy as np
import matplotlib.pyplot as plt

# Read in data from output file
data = np.loadtxt('/Users/daniel/Desktop/output.txt', skiprows=1, delimiter=',')
time = data[:,0]
timestep = data[:,1]
x_star = data[:,2]
y_star = data[:,3]
vx_star = data[:,4]
vy_star = data[:,5]
x_planet = data[:,6]
y_planet = data[:,7]
vx_planet = data[:,8]
vy_planet = data[:,9]

# Calculate total energy of system
G = M = 1.0
m = 0.01

r = np.sqrt((x_star - x_planet)**2 + (y_star - y_planet)**2)
v_star = np.sqrt(vx_star**2 + vy_star**2)
v_planet = np.sqrt(vx_planet**2 + vy_planet**2)

KE_star = 0.5 * M * v_star**2
KE_planet = 0.5 * m * v_planet**2
PE = -G * M * m / r

energy = KE_star + KE_planet + PE

# Plot energy over time
plt.plot(time, energy)
plt.xlabel('Time')
plt.ylabel('Total Energy')
plt.show()
*/
