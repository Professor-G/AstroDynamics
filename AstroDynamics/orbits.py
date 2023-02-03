import time
import numpy as np  
import matplotlib.pyplot as plt  
from cycler import cycler
from pathlib import Path
from progress import bar  

class orbit:
  """Numerical techniques for solving the two-body problem.

  The methods calculate the orbital elements of the system, including the energy and angular momentum at each timestep.

  Note: The inputs must be in code units.

  Args:
      M (float): Mass of the star.
      m (float): Mass of the planet.
      X (np.ndarray): Initial position vector of the star.
      x (np.ndarray): Initial position vector of the planet.
      V (np.ndarray): Initial velocity vector of the star.
      v (np.ndarray): Initial velocity vector of the planet.
      dt (float): Timestep to use for the integration.
      tend (int): Number of timesteps.
      integrator (str): Integrator to use, options include 'euler' and ...
      approx (bool): If True then the acceleration will be calculated using the 1/r^3 approximation
        in which the unit vector is abosorbed by the denominator. Defaults to True.

  Attributes:
      integration_time (float): The integration duration in seconds.
      energy (np.ndarray): Energy of the system, which should be conserved at all timesteps.
      energy_error (np.ndarray): Energy error of the system.
      h (np.ndarray): Angular momentum of the system, as a function of time.
      h_error (np.ndarray): Angular momentum error of the system, as a function of time.
      X_vec (np.ndarray): X-position of the star as a function of integrated time.
      Y_vec (np.ndarray): Y-position of the star as a function of integrated time.
      x_vec (np.ndarray): X-position of the planet as a function of integrated time.
      y_vec (np.ndarray): Y-position of the planet as a function of integrated time.
      path (str): Absolute path to the directory where the images should be saved.

  """

  def __init__(self, M, m, X, x, V, v, dt, tend, integrator, approx=True):
    self.M = M  
    self.m = m  
    self.X = X 
    self.x = x  
    self.V = V 
    self.v = v
    self.dt = dt  
    self.tend = tend 
    self.integrator = integrator 
    self.approx = approx

    self.init_params = np.c_[X, x, V, v]
    self.path = str(Path.home()) + '/'
    self._run_()

  def _run_(self):
    """Calculates the orbital characteristics.

    Args:
        None

    Returns:
        None

    Raises:
        ValueError: If `self.integrator` is not a string or is not 'euler'

    This method must be called directly if the class parameters
    are updated after initialization.
    """

    if isinstance(self.integrator, str) is False:
      raise ValueError("integrator parameter must be 'euler' or !")

    if self.integrator == 'euler':
      self.euler_integrator()

  def euler_integrator(self):
    """Executes the Euler integration method and calculates time taken to complete.
    
    Args:  
        None

    Returns:
        None

    Saves time taken to complete in `self.integration_time`.
    """

    #To ensure the initial conditions are reset
    self.X, self.x, self.V, self.v = np.transpose(self.init_params)

    #Postion vectors to save values
    self.X_vec, self.x_vec, self.Y_vec, self.y_vec = [],[],[],[]

    self.energy, self.h, self.com_pos, self.com_vel = [],[],[],[]
    
    self.timesteps = np.arange(0, self.tend, self.dt)
    progess_bar = bar.FillingSquaresBar('Running integrator...', max=len(self.timesteps))
    start_time = time.time()

    for t in self.timesteps:

      #Calculate the acceleration
      a, A = self.calc_acceleration()

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

      #Center of mass calculation
      com = self.calc_com()

      #Velocity of center of mass calculation
      com_vel = self.calc_vcom(Vx=self.V[0], vx=self.v[0], Vy=self.V[1], vy=self.v[1])

      self.energy.append(energy), self.h.append(h), self.com_pos.append(com), self.com_vel.append(com_vel)
      progess_bar.next()

    end_time = time.time()
    progess_bar.finish()
    self.integration_time = end_time - start_time
    print(f'Time to execute: {np.round(self.integration_time, 4)} seconds.')

    self.energy_error = np.abs((np.array(self.energy)-self.energy[0])/self.energy[0])
    self.h_error = np.abs((np.array(self.h)-self.h[0])/self.h[0])

    return 

  def calc_acceleration(self):
    """Calculates the acceleration of both bodies at a given position.

    Args:
        None

    Returns:
        float: Two floats, the acceleration of the planet and the acceleration of the star, respectively.
    """
    if self.approx:     
        #Earth and Sun acceleration
        a = -self.M*(self.x-self.X)/np.sqrt(np.sum(np.square(self.x - self.X)))#np.linalg.norm(self.x-self.X)**3
        A = -self.m*(self.X-self.x)/np.sqrt(np.sum(np.square(self.X - self.x)))#np.linalg.norm(self.X-self.x)**3 
    else: #Need to calculate the unit vectors
        r1, r2 = np.sqrt(np.sum(np.square(self.x - self.X))), np.sqrt(np.sum(np.square(self.X - self.x)))
        u1, u2 = (self.x - self.X) / r1, (self.X - self.x) / r2
        a = -self.M * u1 / (r1**2)
        A = -self.m * u2 / (r2**2)

    return a, A 

  def calc_energy(self, r, Vx, vx, Vy, vy):
    """Calculates the energy of the two-body system, assuming z-components are zero!

    Args:
        r (float): The length of the vector connecting both bodies.
        Vx (float): The x-component of the star's velocity vector.
        vx (float): The x-component of the planet's velocity vector.
        Vy (float): The y-component of the star's velocity vector.
        vy (float): The y-component of the planet's velocity vector.

    Returns:
        float: The energy of the two-body system at the input positions. 
    """

    Vtot, vtot = np.sqrt(Vx**2 + Vy**2), np.sqrt(vx**2 + vy**2)
    energy = self.m*(vtot**2)/2.0 + self.M*(Vtot**2)/2.0 - self.m*self.M/r
    
    return energy

  def calc_com(self):
    """Calculates the norm of the center of mass vector at a given timestep.
    
    Args:
        None

    Returns:
        float: The length of the center of mass position vector.
    """

    tot_mass = (self.M + self.m)

    x = (self.M * self.X[0] + self.m * self.x[0]) / tot_mass
    y = (self.M * self.X[1] + self.m * self.x[1]) / tot_mass
    z = (self.M * self.X[2] + self.m * self.x[2]) / tot_mass

    r = np.array([x,y,z])

    return np.sqrt(np.dot(r, r))  

  def calc_vcom(self, Vx, vx, Vy, vy):
    """Calculates the velocity of the center of mass of the two objects.

    Args:
        Vx (float): The x-component of the star's velocity vector.
        vx (float): The x-component of the planet's velocity vector.
        Vy (float): The y-component of the star's velocity vector.
        vy (float): The y-component of the planet's velocity vector.

    Returns:
        float: The velocity of the center of mass.
    """

    V, v = np.array([Vx, Vy]), np.array([vx, vy])
    total_momentum = self.M * V + self.m * v

    #p=mv
    vel = total_momentum / (self.M + self.m)  
    
    return np.sqrt(np.dot(vel, vel))

  def calc_momentum(self, r, Vx, Vy):
    """Calculates the angular momentum of the system.

    Args:
        r (float): The length of the vector connecting both bodies.
        Vx (float): The x-component of the star's velocity vector.
        Vy (float): The y-component of the star's velocity vector.

    Returns:
        float: The angular momentum at the specified timestamp.
    """

    phi = np.arctan((self.x[1] - self.X[1]) / self.x[0] - self.X[0])
    V_phi = -Vx*np.sin(phi) + Vy*np.cos(phi)
    phi_r = V_phi / r 
    h = r**2 * phi_r 

    return h 

  def plot_com(self, savefig=False):
    """Plots the COM position of the star and the planet's orbit.

    Args:
        savefig (bool, optional): If True, the image will be saved to the home directory, 
            unless a path attribute is set. Defaults to False, which will output the figure.

    Returns:
        AxesImage: The resulting plot.
    """

    plt.plot(self.timesteps, self.com_pos, 'blue', marker = '*')
    plt.ylabel('C.O.M. Position'), plt.xlabel('Time')
    plt.yscale('log')
    if savefig is False:
      plt.show()
    else:
      _set_style_()
      plt.savefig(self.path+'com_plot.png', bbox_inches='tight', dpi=300)
      plt.clf(), plt.style.use('default')
      print('Image saved in: {}'.format(self.path))

    return 

  def plot_vcom(self, savefig=False):
    """Plots the COM velocity of the star and the planet's orbit.

    Args:
        savefig (bool, optional): If True, the image will be saved to the home directory, 
            unless a path attribute is set. Defaults to False, which will output the figure.

    Returns:
        AxesImage: The resulting plot.
    """

    plt.plot(self.timesteps, self.com_vel, 'blue', marker = '*')
    plt.ylabel('C.O.M. Velocity'), plt.xlabel('Time')
    plt.yscale('log')
    if savefig is False:
      plt.show()
    else:
      _set_style_()
      plt.savefig(self.path+'vcom_plot.png', bbox_inches='tight', dpi=300)
      plt.clf(), plt.style.use('default')
      print('Image saved in: {}'.format(self.path))

    return 

  def plot_orbit(self, savefig=False):
    """Plots the XY position of the star and the planet's orbit.

    Args:
        savefig (bool, optional): If True, the image will be saved to the home directory, 
            unless a path attribute is set. Defaults to False, which will output the figure.

    Returns:
        AxesImage: The resulting plot.
    """

    plt.plot(self.x_vec, self.y_vec, 'blue', label = 'Earth')
    plt.plot(self.X_vec, self.Y_vec, 'orange', label = "Sun")
    plt.xlabel("X Position"), plt.ylabel("Y Position")
    plt.legend()
    if savefig is False:
      plt.show()
    else:
      _set_style_()
      plt.savefig(self.path+'orbit_plot.png', bbox_inches='tight', dpi=300)
      plt.clf(), plt.style.use('default')
      print('Files saved in: {}'.format(self.path))

    return 

  def plot_energy(self, savefig=False):
    """Plots the energy of the system as a function of the integration timesteps.

    Args:
        savefig (bool, optional): If True, the image will be saved to the home directory,
            unless a path attribute is set. Defaults to False, which will output the figure.

    Returns:
        AxesImage: The resulting plot.
    """

    plt.plot(self.timesteps, self.energy_error, 'blue', marker = '*')
    plt.ylabel(r'$\Delta \rm E / \rm E$')
    plt.xlabel('Time')
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
    """Plots the momentum of the system as a function of the integration timesteps.

    Args:
        savefig (bool, optional): If True, the image will be saved to the home directory,
            unless a path attribute is set. Defaults to False, which will output the figure.

    Returns:
        AxesImage: The resulting plot.
    """

    plt.plot(self.timesteps, self.h_error, 'blue', marker = '*')
    plt.ylabel(r'$\Delta \rm h / \rm h$')
    plt.xlabel('Time')
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
  """Function to configure the matplotlib.pyplot style.
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

  return

