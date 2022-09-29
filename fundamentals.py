import numpy as np
from constants import load_constants
cd = load_constants()  # set constants dictionary as global


class PolarizationFcts:
    @ staticmethod
    def real_alpha(d1_line_strength, d2_line_strength, linewidth_d1, linewidth_d2, omega_d1, omega_d2, omega):
        d1_part = d1_line_strength * linewidth_d1 / omega_d1 ** 3 \
                  * (1 / (omega_d1 - omega)
                     + 1 / (omega_d1 + omega))

        d2_part = d2_line_strength * linewidth_d2 / omega_d2 ** 3 \
                  * (1 / (omega_d2 - omega)
                     + 1 / (omega_d2 + omega))

        return 3 * np.pi * cd['epsilon_0'][0] * cd['c'][0]**3 * (d1_part + d2_part)

    @staticmethod
    def imag_alpha(d1_line_strength, d2_line_strength, linewidth_d1, linewidth_d2, omega_d1, omega_d2, omega):
        d1_part = d1_line_strength * linewidth_d1**2 / omega_d1**6 * (1 / (omega_d1 - omega) + 1 / (omega_d1 + omega))**2

        d2_part = d2_line_strength * linewidth_d2**2 / omega_d2**6 * (1 / (omega_d2 - omega) + 1 / (omega_d2 + omega))**2

        return 3 * np.pi * cd['epsilon_0'][0] * cd['c'][0]**3 * omega**3 / 2 * (d1_part + d2_part)

    @staticmethod
    def dipole_potential(d1_line_strength, d2_line_strength, linewidth_d1, linewidth_d2, omega_d1, omega_d2, omega, I):
        real_alpha = PolarizationFcts.real_alpha(d1_line_strength, d2_line_strength, linewidth_d1, linewidth_d2,
                                                 omega_d1, omega_d2, omega)
        return - 1 / (2 * cd['epsilon_0'][0] * cd['c'][0]) * real_alpha * I
