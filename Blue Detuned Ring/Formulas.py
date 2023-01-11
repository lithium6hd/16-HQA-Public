# ------------------------------- load class from separate file -------------------------------

# from Gaussian_ODT import Gauss_beam
import numpy as np
import scipy.constants
from scipy.special import erf
from scipy import integrate
import matplotlib.pyplot as plt


# ------------------------------- Formulas --------------------------------

class ThermalFacts:
    """
    TODO
    """

    @staticmethod
    def boltzmann(v, temperature, mass):
        """ Boltzmann distribution in 3D v = velocity; In 3D -> dv = dvx + dvy + dvz
        :param v: thermal velocity
        :param temperature: of gas
        :param mass: atomic mass
        :return: probability of the certain state with v
        """
        a = mass / (2 * scipy.constants.k * temperature)
        return 4 * np.pi * (a / np.pi)**(3 / 2) * v**2 * np.exp(-1 * v**2 * a)

    @staticmethod
    def boltzmann_energy(e, temperature):
        """ Boltzmann distribution in 3D v = velocity; In 3D -> dv = dvx + dvy + dvz
        :param e: thermal energy
        :param temperature: of gas
        :return: probability of the certain state with e
        """
        # a = mass / (2 * scipy.constants.k * temperature)
        return 2 * np.sqrt(e / np.pi) * (1 / (scipy.constants.k * temperature)) ** (3 / 2) \
               * np.exp(- e / (scipy.constants.k * temperature))

    @staticmethod
    def boltzmann_integral(temperature, mass, lower, upper):
        """
        Boltzmann velocity distribution
        :param temperature: of gas [K]
        :param mass: particle mass [kG]
        :param lower: lower bound of integration [m/s]
        :param upper: upper bound of integration [m/s]
        :return: Integral of boltzmann velocity distribution [value, error]
        """
        return integrate.quad(lambda v: ThermalFacts.boltzmann(v, temperature, mass), lower, upper)

    @staticmethod
    def boltzmann_newtemp(temperature, temperature_cut):
        """
        New thermalized temperature of truncated boltzmann distribution
        :param temperature: of gas [K]
        :param temperature_cut: truncated boltzmann distribution, trap depth in [K]
        :return: new temperature after evaporation [K]
        """
        root = np.sqrt(temperature_cut / temperature)
        return temperature / 3 * \
               ((3 * np.sqrt(np.pi) * erf(root) - (6 * root + 4 * (temperature_cut / temperature) ** (3 / 2)) *
                 np.exp(-temperature_cut / temperature)) /
                (np.sqrt(np.pi) * erf(root) - 2 * root * np.exp(-temperature_cut / temperature)))

    @staticmethod
    def v_therm(temperature, mass):
        """
        Thermal velocity in 3D
        :param temperature: of gas
        :param mass: atomic mass
        :return: thermal velocity [m/s]
        """
        return np.sqrt(2 * scipy.constants.k * temperature / mass)  # in 3 dimensions

    @staticmethod
    def lambda_dB(temperature, mass):
        """
        De Broglie wavelength
        :param temperature: of gas [K]
        :param atomic: mass [kG]
        :return: De Broglie wavelength [m]
        """
        return scipy.constants.h / np.sqrt(2 * np.pi * mass * scipy.constants.k * temperature)

    @staticmethod
    def cross_section(temperature, mass):
        """
        Cross section: Assumed to only depend on the De Broglie wavelength!
        :param temperature: of gas [K]
        :param mass: atomic mass [kG]
        :return: Scattering cross section [m^2]
        """
        return 4 * np.pi * ThermalFacts.lambda_dB(temperature, mass) ** 2

    @staticmethod
    def scattering_rate(temperature, mass, density):
        """
        Collision rate of particles
        :param temperature: of gas [K]
        :param mass: atomic mass [kg]
        :param density: particle number density [1/m^3]
        :return: scattering rate of thermal movement [1/s]
        """
        return ThermalFacts.v_therm(temperature, mass) * ThermalFacts.cross_section(temperature, mass) * density
