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
    elif self.integrator == 'runge-kutta':
      self.runge_kutta_integrator()

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

    self.energy, self.h, com_x, com_vx = [],[],[],[]

    self.timesteps = np.arange(0, self.tend, self.dt)
    progess_bar = bar.FillingSquaresBar('Running Euler Integrator...', max=len(self.timesteps))
    start_time = time.time()

    for t in self.timesteps:

      #Calculate the acceleration
      a, A = self.calc_acceleration(x=self.x, X=self.X)

      #Update positions and velocities at each timestamp
      self.x = self.x + self.v*self.dt
      self.X = self.X + self.V*self.dt
      self.v = self.v + a*self.dt
      self.V = self.V + A*self.dt

      self.x_vec.append(self.x[0]), self.y_vec.append(self.x[1])
      self.X_vec.append(self.X[0]), self.Y_vec.append(self.X[1])

      #Energy calculation
      energy = self.calc_energy(Vx=self.V[0], vx=self.v[0], Vy=self.V[1], vy=self.v[1])
      
      #Momentum calculation
      h = self.calc_momentum(Vx=self.V[0], vx=self.v[0], Vy=self.V[1], vy=self.v[1])

      #Center of mass calculation
      comx = self.calc_com()

      #Velocity of center of mass calculation
      comvx = self.calc_comv(Vx=self.V[0], vx=self.v[0], Vy=self.V[1], vy=self.v[1])

      com_x.append(comx), com_vx.append(comvx), self.energy.append(energy), self.h.append(h)

      progess_bar.next()

    end_time = time.time()
    progess_bar.finish()
    self.integration_time = end_time - start_time
    print(f'Time to execute: {np.round(self.integration_time, 4)} seconds.')

    self.energy_error = np.abs((np.array(self.energy)-self.energy[0])/self.energy[0])
    self.h_error = np.abs((np.array(self.h)-self.h[0])/self.h[0])
    self.com, self.comv = com_x, com_vx

    return 

  def runge_kutta_integrator(self):
    """Executes the Runge-Kutta integration method and calculates time taken to complete.
    
    Args:  
        None

    Returns:
        None

    Saves time taken to complete in `self.integration_time`.
    """

    #To ensure the initial conditions are reset
    self.X, self.x, self.V, self.v = np.transpose(self.init_params)
    # Postion vectors to save values
    self.X_vec, self.x_vec, self.Y_vec, self.y_vec = [],[],[],[]

    self.energy, self.h, com_x, com_vx = [],[],[],[]

    self.timesteps = np.arange(0, self.tend, self.dt)
    progess_bar = bar.FillingSquaresBar('Running Runge-Kutta Integrator...', max=len(self.timesteps))
    start_time = time.time()

    for t in self.timesteps:
      
      #Calculate the acceleration
      a, A = self.calc_acceleration(x=self.x, X=self.X)

      # Intermediate velocities and positions using k1 values
      x1 = self.x + self.v*self.dt/2
      X1 = self.X + self.V*self.dt/2
      V1 = self.V + A*self.dt/2
      v1 = self.v + a*self.dt/2
      
      # Recalculate acceleration using intermediate positions and velocities
      a1, A1 = self.calc_acceleration(x=x1, X=X1)

      # Intermediate velocities and positions using k2 values
      x2 = self.x + v1*self.dt/2
      X2 = self.X + V1*self.dt/2
      v2 = self.v + a1*self.dt/2
      V2 = self.V + A1*self.dt/2
      
      # Recalculate acceleration using intermediate positions and velocities
      a2, A2 = self.calc_acceleration(x=x2, X=X2)

      # Intermediate velocities and positions using k3 values
      x3 = self.x + v2*self.dt
      X3 = self.X + V2*self.dt
      v3 = self.v + a2*self.dt
      V3 = self.V + A2*self.dt
    
      # Recalculate acceleration using intermediate positions and velocities
      a3, A3 = self.calc_acceleration(x=x3, X=X3)

      self.x = self.x + (self.v + 2*v1 + 2*v2 + v3)*self.dt / 6.
      self.X = self.X + (self.V + 2*V1 + 2*V2 + V3)*self.dt / 6.
      self.v = self.v + (a + 2*a1 + 2*a2 + a3)*self.dt / 6.
      self.V = self.V + (A + 2*A1 + 2*A2 + A3)*self.dt / 6.
      
      self.x_vec.append(self.x[0]), self.y_vec.append(self.x[1])
      self.X_vec.append(self.X[0]), self.Y_vec.append(self.X[1])

      # Energy calculation
      energy = self.calc_energy(Vx=self.V[0], vx=self.v[0], Vy=self.V[1], vy=self.v[1])
      
      # Momentum calculation
      h = self.calc_momentum(Vx=self.V[0], vx=self.v[0], Vy=self.V[1], vy=self.v[1])

      # Center of mass calculation
      comx = self.calc_com()

      #Velocity of center of mass calculation
      comvx = self.calc_comv(Vx=self.V[0], vx=self.v[0], Vy=self.V[1], vy=self.v[1])

      com_x.append(comx), com_vx.append(comvx), self.energy.append(energy), self.h.append(h)

      progess_bar.next()

    end_time = time.time()
    progess_bar.finish()
    self.integration_time = end_time - start_time
    print(f'Time to execute: {np.round(self.integration_time, 4)} seconds.')

    self.energy_error = np.abs((np.array(self.energy)-self.energy[0])/self.energy[0])
    self.h_error = np.abs((np.array(self.h)-self.h[0])/self.h[0])
    self.com, self.comv = com_x, com_vx

    return 

  def calc_acceleration(self, x, X):
    """Calculates the acceleration of both bodies at a given position.

    Args:
        None

    Returns:
        float: Two floats, the acceleration of the planet and the acceleration of the star, respectively.
    """
    if self.approx:     
        #Earth and Sun acceleration
        a = -self.M*(x-X)/np.sqrt(np.sum(np.square(x-X)))#np.linalg.norm(self.x-self.X)**3
        A = -self.m*(X-x)/np.sqrt(np.sum(np.square(X-x)))#np.linalg.norm(self.X-self.x)**3 
    else: #Need to calculate the unit vectors
        r1, r2 = np.sqrt(np.sum(np.square(x-X))), np.sqrt(np.sum(np.square(X-x)))
        u1, u2 = (x-X) / r1, (X-x) / r2
        a = -self.M * u1 / (r1**2)
        A = -self.m * u2 / (r2**2)

    return a, A 

  def calc_energy(self, Vx, vx, Vy, vy):
    """Calculates the energy of the two-body system, assuming z-components are zero!

    Args:
        Vx (float): The x-component of the star's velocity vector.
        vx (float): The x-component of the planet's velocity vector.
        Vy (float): The y-component of the star's velocity vector.
        vy (float): The y-component of the planet's velocity vector.

    Returns:
        float: The energy of the two-body system at the input positions. 
    """
    
    #The length of the vector connecting both bodies.
    r = np.sqrt((self.x[0] - self.X[0])**2 + (self.x[1] - self.X[1])**2)

    Vtot, vtot = np.sqrt(Vx**2 + Vy**2), np.sqrt(vx**2 + vy**2)
    energy = self.m*(vtot**2)/2.0 + self.M*(Vtot**2)/2.0 - self.m*self.M/r
    
    return energy

  def calc_com(self):
    """Calculates the position vector of the center of mass vector at a given timestep.
    
    Args:
        None

    Returns:
        float: The length of the center of mass position vector.
    """

    tot_mass = self.M + self.m

    rp = np.sqrt(self.x[0]**2 + self.x[1]**2)
    rs = np.sqrt(self.X[0]**2 + self.X[1]**2)

    r = (self.M * rs + self.m * rp) / tot_mass

    return r

  def calc_comv(self, Vx, vx, Vy, vy):
    """Calculates the velocity vector of the center of mass of the two objects.

    Args:
        Vx (float): The x-component of the star's velocity vector.
        vx (float): The x-component of the planet's velocity vector.
        Vy (float): The y-component of the star's velocity vector.
        vy (float): The y-component of the planet's velocity vector.

    Returns:
        float: The velocity of the center of mass.
    """

    tot_mass = self.M + self.m

    vp = np.sqrt(vx**2 + vy**2)
    vs = np.sqrt(Vx**2 + Vy**2)

    v = (self.M * vs + self.m * vp) / tot_mass
         
    return v

  def calc_momentum(self, Vx, vx, Vy, vy):
    """Calculates the angular momentum of the system.

    Args:
        vx (float): The x-component of the planet's velocity vector.
        vy (float): The y-component of the planet's velocity vector.

    Returns:
        float: The angular momentum at the specified timestamp.
    """

    x, y = self.x[0] - self.X[0], self.x[1] - self.X[1]
    v_x, v_y = vx - Vx, vy - Vy

    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)

    vrad =  v_x*np.cos(phi) + v_y*np.sin(phi)
    vphi = -v_x*np.sin(phi) + v_y*np.cos(phi)

    phidot = vphi / r
    h = r**2 * phidot 

    return h 

  def plot_com(self, savefig=False):
    """Plots the COM position of the star and the planet's orbit.

    Args:
        savefig (bool, optional): If True, the image will be saved to the home directory, 
            unless a path attribute is set. Defaults to False, which will output the figure.

    Returns:
        AxesImage: The resulting plot.
    """

    plt.plot(self.timesteps, self.com, 'blue', marker = '*')
    plt.xlabel('Time'), plt.ylabel(r'$r_{COM}$')
    plt.title('Center of Mass')
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

    plt.plot(self.timesteps, self.comv, 'blue', marker = '*')
    plt.xlabel('Time'), plt.ylabel(r'$V_{COM}$')
    plt.yscale('log')
    plt.title('Center of Mass')
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
    plt.xlabel("X"), plt.ylabel("Y")
    plt.legend()
    plt.title('Orbits')
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
    plt.xlabel('Time'), plt.ylabel(r'|$\Delta \rm E / \rm E_0|$')
    plt.yscale('log')
    plt.title('Error Growth')
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
    plt.ylabel(r'|$\Delta \rm h / \rm h_0$|')
    plt.xlabel(r'$\rm Time \ (\Omega = 1)$')
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
  """Function to configure the matplotlib.pyplot style. This function is called before any images are saved.
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

