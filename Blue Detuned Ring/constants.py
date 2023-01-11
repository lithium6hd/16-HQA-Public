from scipy.constants import k, c, h, hbar, g, epsilon_0, physical_constants, tera, giga, mega, kilo, centi, milli, \
    micro, nano
from numpy import pi

def load_constants():
    const_dict = dict()

    # general constants
    const_dict['kB'] = [k, 'J/K', 'Boltzmann constant']
    const_dict['c'] = [c, 'm/s', 'speed of light in vacuum']
    const_dict['h'] = [h, 'J*s', 'Planck constant']
    const_dict['hbar'] = [hbar, 'J*s', 'reduced Planck constant']
    const_dict['g'] = [g, 'm/s^2', 'gravitational acceleration']
    const_dict['epsilon_0'] = [epsilon_0, 'F/m', 'vacuum permittivity']

    const_dict['u0'] = list(physical_constants['atomic mass constant'][:2]) + ['atomic mass constant']
    const_dict['a0'] = list(physical_constants['Bohr radius'][:2]) + ['Bohr radius']
    const_dict['uB'] = list(physical_constants['Bohr magneton'][:2]) + ['Bohr magneton']

    # 6Li D1 transition (2 2S_1/2 -> 2 2P_1/2)
    const_dict['omegaD1_6Li'] = [2 * pi * 446.789634 * tera, '1/s',
                                 'angular frequency of D1 transition of 6Li {Properties of 6Li}']
    const_dict['linewidthD1_6Li'] = [2 * pi * 5.8724 * mega, '1/s',
                                     'Natural Line Width (FWHM)/decay rate of D1 transition in angular frequency (6Li) {Properties of 6Li}']

    # 6Li D2 transition (2 2S_1/2 -> 2 2P_3/2)
    const_dict['omegaD2_6Li'] = [2 * pi * 446.799677 * tera, '1/s',
                                 'angular frequency of D2 transition of 6Li {Properties of 6Li}']
    const_dict['linewidthD2_6Li'] = [2 * pi * 5.8724 * mega, '1/s',
                                     'Natural Line Width (FWHM)/decay rate of D2 transition in angular frequency (6Li) {Properties of 6Li}']

    return const_dict
