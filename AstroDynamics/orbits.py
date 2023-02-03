import time
import numpy as np  
import matplotlib.pyplot as plt  
from cycler import cycler
from pathlib import Path
from progress import bar  

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
    integrator (str): Integrator to use.

  Attributes:
    _run_ : Calculates the orbital characteristics, this is run
      once when a class object is initialized, and must be called manually
      if the class parameters are updated. 

    euler_integrator (ndarray): Integrates using the Euler method.

    plot_orbit (AxesImage): Plots the data_x parameter space using t-SNE

    plot_energy (AxesImage): Plots the confusion matrix, assessed with data_x.

    plot_momentum (AxesImage): Plots ROC curve, assessed with data_x

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

    self.path = str(Path.home()) + '/'
    self._run_()

  def _run_(self):
    """
    Function that calls and executes using the assigned integrator.
    This function must be called directly if the class parameters
    are updated after initialization. 
    """
    if isinstance(self.integrator, str) is False:
      raise ValueError('No integrator input!')

    if self.integrator == 'euler':
      self.euler_integrator()

  def euler_integrator(self):
    """
    Excecutes the Euler integration method
    and prints out the time taken to complete. The
    time is saved and can be loaded using .integration_time
    """

    #Energy quantity and postion vectors to save values
    self.energy, self.h, self.X_vec, self.x_vec, self.Y_vec, self.y_vec = [],[],[],[],[],[]

    timesteps = np.arange(0, self.tend, self.dt)
    progess_bar = bar.FillingSquaresBar('Running integrator...', max=len(timesteps))
    start_time = time.time()

    for t in timesteps:

      #Earth's acceleration
      a = -self.M*(self.x-self.X)/np.linalg.norm(self.x-self.X)**3 
      
      #Sun's acceleration
      A = -self.m*(self.X-self.x)/np.linalg.norm(self.X-self.x)**3 

      #Update positions and velocities at each timestamp
      self.x = self.x + self.v*self.dt
      self.X = self.X + self.V*self.dt
      self.v = self.v + a*self.dt
      self.V = self.V + A*self.dt

      self.x_vec.append(self.x[0]), self.y_vec.append(self.x[1])
      self.X_vec.append(self.X[0]), self.Y_vec.append(self.X[1])
      r = np.sqrt((self.x[0] - self.X[0])**2 + (self.x[1] - self.X[1])**2)
      
      #Energy calculation
      energy = self.calc_energy(r=r, Vx=self.V[0], vx=self.v[0], Vy=self.V[1], vy=self.v[1])
      
      #Momentum calculation
      h = self.calc_momentum(r=r, Vx=self.V[0], Vy=self.V[1])

      self.energy.append(energy), self.h.append(h)
      progess_bar.next()

    end_time = time.time()
    progess_bar.finish()
    self.integration_time = end_time - start_time
    print(f'Time to execute: {np.round(self.integration_time, 4)} seconds.')

    return 

  def calc_energy(self, r, Vx, vx, Vy, vy):
    """
    Calculates the energy of the system using the velocity
    and position vectors.

    Args:
      r (float): The length of the vector connecting both bodies.
      Vx (float): The x-component of the star's velocity vector.
      vx (float): The x-component of the planet's velocity vector.
      Vy (float): The y-component of the star's velocity vector.
      vy (float): The y-component of the planet's velocity vector.

    Returns:
      Energy of the two-body system at the input positions. 
    """

    Vtot, vtot = np.sqrt(Vx**2 + Vy**2), np.sqrt(vx**2 + vy**2)
    energy = self.m*(vtot**2)/2.0 + self.M*(Vtot**2)/2.0 - self.m*self.M/r
    
    return energy

  def calc_momentum(self, r, Vx, Vy):
    """
    Calculates the angular momentum given a set of position and velocity vectors.

    Args:
      r (float): The length of the vector connecting both bodies.
      Vx (float): The x-component of the star's velocity vector.
      Vy (float): The y-component of the star's velocity vector.

    Returns:
      The angular momentum at the specified timestamp.
    """

    phi = np.arctan((self.x[1] - self.X[1]) / self.x[0] - self.X[0])
    V_phi = -Vx*np.sin(phi) + Vy*np.cos(phi)
    phi_r = V_phi / r 
    h = r**2 * phi_r 

    return h 

  def plot_orbit(self, savefig=False):
    """
    Plots the X & Y positions of the star and planet's orbit.

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
      _set_style_()
      plt.savefig(self.path+'orbit_plot.png', bbox_inches='tight', dpi=300)
      plt.clf(), plt.style.use('default')
      print('Files saved in: {}'.format(self.path))

    return 

  def plot_energy(self, savefig=False):
    """
    Plots the energy of the system as a function of the
    integration timesteps.

    Returns:
      AxesImages
    """

    energy_error = np.abs((np.array(self.energy)-self.energy[0])/self.energy[0])

    plt.figure(figsize =(10,10))
    plt.plot(np.arange(0, self.tend, self.dt), energy_error, 'blue', marker = '*')
    plt.ylabel(r'$\Delta \rm E / \rm E$', fontsize = 17)
    plt.xlabel('Time', fontsize = 17)
    plt.yscale('log')
    if savefig is False:
      plt.show()
    else:
      _set_style_()
      plt.savefig(self.path+'energy_plot.png', bbox_inches='tight', dpi=300)
      plt.clf(), plt.style.use('default')
      print('Image saved in: {}'.format(self.path))

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
    plt.plot(np.arange(0, self.tend, self.dt), h_error, 'blue', marker = '*')
    plt.ylabel(r'$\Delta \rm h / \rm h$', fontsize = 17)
    plt.xlabel('Time', fontsize = 17)
    plt.yscale('log')
    if savefig is False:
      plt.show()
    else:
      _set_style_()
      plt.savefig(self.path+'momentum_plot.png', bbox_inches='tight', dpi=300)
      plt.clf(), plt.style.use('default')
      print('Image saved in: {}'.format(self.path))

    return 

def _set_style_():
  """
  Function to configure the matplotlib.pyplot style.
  This function is called before any images are saved.
  """
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
  plt.rcParams["figure.autolayout"] = True
  plt.rcParams["figure.dpi"] = 300

