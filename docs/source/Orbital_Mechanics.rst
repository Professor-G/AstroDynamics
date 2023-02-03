.. _Orbital_Mechanics:

Orbital Mechanics
==================

In this excercise the two-body system of the Sun and Earth is used to compare the results of different integrators. The code is written in code units and therefore the mass of the Sun, :math:`M`, is assumed to be 1, the mass of the Earth, :math:`m`, is :math:`3 \times 10^{-6}`, the gravitational constant :math:`G` is assumed to be 1, and the initial conditions for the position and velocity are given as 3-dimensional vectors on an XY plane. The necessary modules and pamereters can be set as follows:

.. code-block:: python
	
    from AstroDynamics import orbits
    import numpy as np

    #capital letters = SUN, lower case = EARTH
    M, m = 1.0, 3.0e-6

    #Position vectors
    X, x = np.array([0., 0., 0.]), x = np.array([1., 0., 0.])

    #Velocity vectors
    V, v = np.array([0., 0., 0.]), np.array([0., 1., 0.])
    
    #Integration parameters
    dt, tend = 1e-3, 100


where :math:`dt` and :math:`tend` correspond to the timestep and the number of steps to take during the integration, respectively. The orbit object can then be created using the `orbit class<https://astrodynamics.readthedocs.io/en/latest/autoapi/astrodynamics/orbits/index.html>`_. 

Different integrators are useful for different types of systems. For example, the fourth-order Hermite integrator is optimal for star cluster dynamics while symplectic integrators such as the Wisdom-Holman integrator are often used for secular evolution of planetary systems. The simplest integrator is the forward Euler integration scheme, which calculates the expected position and velocity of an object by assuming it drifts at a uniform velocity during a small time interval. 


Euler Integrator
------------------


.. code-block:: python
	
    orbit = orbits.orbit(M=M, m=m, X=X, V=V, x=x, v=v, dt=dt, tend=tend, integrator='euler')


