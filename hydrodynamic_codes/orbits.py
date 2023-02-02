import numpy as np  
import matplotlib.pyplot as plt  
import time

class orbit:
  """
  Calculates orbital characteristics

  Note:
    The inputs should be in code units!

  Args:
    M (float): Mass of the star.
    m (float): Mass of the planet.
    X (ndarray): Initial position vector of the star
    x (ndarray): Initial position vector of the planet.
    V (ndarray): Initial velocity vector of the star.
    v (ndarray): Initial velocity vector of the planet.
    dt (float): Timestep, defaults to 0.001.
    tend (float): Number of timesteps, defaults to 1000.
    integrator (str): Integrator to use, defaults to 'Euler'.

  """

  def __init__(self, M, m, X, x, V, v, dt, tend, integrator):
    self.M = M  
    self.m = m  
    self.X = X 
    self.x = x  
    self.V = V 
    self.v = v 
    self.integrator = integrator 
    self.dt = dt  
    self.tend = tend 

    self._run_()

  def _run_(self):
    if isinstance(self.integrator, str) is False:
      raise ValueError('No integrator input!')
    if self.integrator == 'euler':
      self.euler_integrator()

  def euler_integrator(self):
    """
    Excecutes the Euler integration method
    and prints out the time taken to complete

    """

    #Energy quantity and postion vectors to save values
    Energy, ang_momentum, x_vec, y_vec, X_vec, Y_vec = [], [],[],[],[],[]

    start_time = time.time()
    for t in np.arange(0, self.tend, self.dt):

      #Earth's acceleration
      a = -self.M*(self.x-self.X)/np.linalg.norm(self.x-self.X)**3 
      
      #Sun's acceleration
      A = -self.m*(self.X-self.x)/np.linalg.norm(self.X-self.x)**3 

      #update positions and velocities
      self.x = self.x + self.v*self.dt
      self.X = self.X + self.V*self.dt
      self.v = self.v + a*self.dt
      self.V = self.V + A*self.dt

      x_vec.append(self.x[0]), y_vec.append(self.x[1])
      X_vec.append(self.X[0]), Y_vec.append(self.X[1])

      #Energy calculation
      vx, vy = self.v[0], self.v[1]
      Vx, Vy = self.V[0], self.V[1]
      vtot, Vtot = np.sqrt(vx**2 + vy**2), np.sqrt(Vx**2 + Vy**2)
   
      r = np.sqrt((self.x[0] - self.X[0])**2 + (self.x[1] - self.X[1])**2)
      energy = self.m*(vtot**2)/2.0 + self.M*(Vtot**2)/2.0 - self.m*self.M/r
      
      phi = np.arctan((self.x[1] - self.X[1]) / self.x[0] - self.X[0])
      V_phi = -Vx*np.sin(phi) + Vy*np.cos(phi)
      phi_r = V_phi / r 
      h = r**2 * phi_r 

      Energy.append(energy), ang_momentum.append(h)

    self.x_vec, self.y_vec = x_vec, y_vec
    self.X_vec, self.Y_vec = X_vec, Y_vec 
    self.Energy, self.h = Energy, ang_momentum
    self.time = np.arange(0, self.tend, self.dt)
  
    end_time = time.time()
    self.integration_time = end_time - start_time

    print(f'Time to execute: {np.round(self.integration_time, 3)} seconds.')

    return 


  def plot_orbit(self, savefig=False):
    """
    Plots the X & Y position of the planet's orbit

    Returns:
      AxesImages
    """

    plt.figure(figsize =(10,10))
    plt.plot(self.x_vec, self.y_vec, 'blue', label = 'Earth')
    plt.plot(self.X_vec, self.Y_vec, 'orange', label = "Sun")
    plt.xlabel("X Position", fontsize = 17)
    plt.ylabel("Y Position", fontsize = 17)
    plt.legend(fontsize = 12)
    if savefig is False:
      plt.show()
    else:
      plt.savefig('orbit_plot.png', bbox_inches='tight', dpi=300)
      plt.clf()

    return 

  def plot_energy(self, savefig=False):
    """
    Plots the energy of the system as a function of the
    integration timesteps.

    Returns:
      AxesImages
    """

    energy_error = np.abs((np.array(self.Energy)-self.Energy[0])/self.Energy[0])

    plt.figure(figsize =(10,10))
    plt.plot(self.time, energy_error, 'blue', marker = '*')
    plt.ylabel(r'$\Delta \rm E / \rm E$', fontsize = 17)
    plt.xlabel('Time', fontsize = 17)
    plt.yscale('log')
    if savefig is False:
      plt.show()
    else:
      plt.savefig('energy_plot.png', bbox_inches='tight', dpi=300)
      plt.clf()

    return 

  def plot_momentum(self, savefig=False):
    """
    Plots the energy of the system as a function of the
    integration timesteps.

    Returns:
        AxesImages
    """

    h_error = np.abs((np.array(self.h)-self.h[0])/self.h[0])

    plt.figure(figsize =(10,10))
    plt.plot(self.time, h_error, 'blue', marker = '*')
    plt.ylabel(r'$\Delta \rm h / \rm h$', fontsize = 17)
    plt.xlabel('Time', fontsize = 17)
    plt.yscale('log')
    if savefig is False:
      plt.show()
    else:
      plt.savefig('momentum_plot.png', bbox_inches='tight', dpi=300)
      plt.clf()

    return 

"""

import numpy as np 
from orbits import orbit 

M = 1.0 #capital letters = SUN, lower case = EARTH
m = 3.0e-6
X = np.array([0., 0., 0.])
V = np.array([0., 0., 0.])
x = np.array([1., 0., 0.])
v = np.array([0., 1., 0.])
dt = 1e-3
tend = 100.

o = orbit(M=M, m=m, X=X, V=V, x=x, v=v, dt=dt, tend=tend, integrator='euler')

"""
"""
import matplotlib.pyplot as plt
from cycler import cycler

plt.rcParams["xtick.color"] = "323034"
plt.rcParams["ytick.color"] = "323034"
plt.rcParams["text.color"] = "323034"
plt.rcParams["lines.markeredgecolor"] = "black"
plt.rcParams["patch.facecolor"] = "bc80bd"
plt.rcParams["patch.force_edgecolor"] = True
plt.rcParams["patch.linewidth"] = 0.8
plt.rcParams["scatter.edgecolors"] = "black"
plt.rcParams["grid.color"] = "b1afb5"
plt.rcParams["axes.titlesize"] = 16
plt.rcParams["legend.title_fontsize"] = 12
plt.rcParams["xtick.labelsize"] = 16
plt.rcParams["ytick.labelsize"] = 16
plt.rcParams["font.size"] = 15
plt.rcParams["axes.prop_cycle"] = (cycler('color', ['bc80bd' ,'fb8072', 'b3de69','fdb462','fccde5','8dd3c7','ffed6f','bebada','80b1d3', 'ccebc5', 'd9d9d9']))
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["font.family"] = "STIXGeneral"
plt.rcParams["lines.linewidth"] = 2
plt.rcParams["lines.markersize"] = 6
plt.rcParams["legend.frameon"] = True
plt.rcParams["legend.framealpha"] = 0.8
plt.rcParams["legend.fontsize"] = 13
plt.rcParams["legend.edgecolor"] = "black"
plt.rcParams["legend.borderpad"] = 0.2
plt.rcParams["legend.columnspacing"] = 1.5
plt.rcParams["legend.labelspacing"] = 0.4
plt.rcParams["text.usetex"] = False
plt.rcParams["axes.labelsize"] = 17
plt.rcParams["axes.titlelocation"] = "center"
plt.rcParams["axes.formatter.use_mathtext"] = True
plt.rcParams["axes.autolimit_mode"] = "round_numbers"
plt.rcParams["axes.labelpad"] = 3
plt.rcParams["axes.formatter.limits"] = (-4, 4)
plt.rcParams["axes.labelcolor"] = "black"
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 1
plt.rcParams["axes.grid"] = False
plt.rcParams["axes.spines.right"] = True
plt.rcParams["axes.spines.left"] = True
plt.rcParams["axes.spines.top"] = True
plt.rcParams["figure.titlesize"] = 18
plt.rcParams["autolayout"] = True
plt.rcParams["dpi"] = 300

"""
