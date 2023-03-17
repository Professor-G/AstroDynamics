#include <cmath>
#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <iomanip>  // for std::setprecision

//Function to calculate the acceleration of both bodies, including approx option which will use the r squared approximation (default is true) 
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

// Below is the function to run the Leapfrog Integrator
void leapfrog(std::vector<double>& X, std::vector<double>& x, std::vector<double>& V, std::vector<double>& v, double tend, double dt, double M, double m, std::ofstream& outfile) {
  
  //Set output precision to 16 digits
  outfile << std::setprecision(16);

  //Write the header row, avoid spaces!: time, timestep, x_star, y_star, vx_star, vy_star, x_planet, y_planet, vx_planet, vy_planet
  outfile << "time,timestep,x_star,y_star,vx_star,vy_star,x_planet,y_planet,vx_planet,vy_planet" << std::endl;

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
    V[0] = KICK_x + acc.second * dt / 2.0;
    V[1] = KICK_y + acc.first * dt / 2.0;

    outfile << i*dt << "," << i << "," << X[0] << "," << X[1] << "," << V[0] << "," << V[1] << "," << x[0] << "," << x[1] << "," << v[0] << "," << v[1] << std::endl;
  }
}

int main() {
  std::vector<double> X(2); //Initial position of the star
  X[0] = 1.0;
  X[1] = 0.0;

  std::vector<double> x(2); //Initial position of the planet
  x[0] = 0.0;
  x[1] = 0.0;

  std::vector<double> V(2); //Initial velocity of the star
  V[0] = 0.0;
  V[1] = 0.0;

  std::vector<double> v(2); //Initial velocity of the planet
  v[0] = 0.0;
  v[1] = 1.0;  
  
  //Physical parameters
  double M = 1.0-3e-6;
  double m = 3e-6;
  double tend = 100.0;
  double dt = 0.01;

  //Open output file for writing which will be passed to the leapfrog function
  std::ofstream outfile("/Users/daniel/Desktop/output.csv");

  //Run Leapfrog integration and write data to file
  leapfrog(X, x, V, v, tend, dt, M, m, outfile);

  //Close output file
  outfile.close();

  return 0;
}


// PYTHON SCRIPT TO PLOT THE ENERGIES
/*
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

#Read in data from output file
df = pd.read_csv('/Users/daniel/Desktop/output.csv')

# Calculate total energy of system
G = 1.0
M = 1-3e-6
m = 3e-6

r = np.sqrt((df.x_star - df.x_planet)**2 + (df.y_star - df.y_planet)**2)
v_star = np.sqrt(df.vx_star**2 + df.vy_star**2)
v_planet = np.sqrt(df.vx_planet**2 + df.vy_planet**2)

KE_star = 0.5 * M * v_star**2
KE_planet = 0.5 * m * v_planet**2
PE = -G * M * m / r
energy = KE_star + KE_planet + PE

# Plot energy over time
plt.plot(df.time, energy)
plt.xlabel('Time')
plt.ylabel('Total Energy')
plt.show()
*/
