import numpy as np
from Calculations.ParCa.Constants.constants import load_constants
cd = load_constants()  # set constants dictionary as global

'''
Units convention:
stay in units which are used in textbook formulas as long as possible!
'''


class ParCaError(Exception):
    """
    create error messages
    """
    def __init__(self, message):
        message = "ParCa ERROR: %s"%str(message)
        Exception.__init__(self, message)

class HyperfineFcts:
    @staticmethod  #means you can can use it via HyerfineFct.state_check_F(...) or as usual
    def state_check_F(I, J, F, mF):
        """
        Check if the given state is a physical valid on. If not: ParCaError appear.
        :param I: nuclear angular momentum [1]
        :param F: total angular momentum [1]
        :param mF: z projection of the total angular momentum [1]
        """
        state_I = np.array([1.0 * I, 2.0 * I])
        state_J = np.array([1.0 * J, 2.0 * J])
        state_F = np.array([1.0 * F, 1.0 * mF, 2.0 * F, 2.0 * mF])
        diff = (F - abs(I - J)) * 1.0

        if not (I >= 0 and state_I[1].is_integer() and J >= 0 and state_J[1].is_integer()):
            raise ParCaError("The input state of: 'I + J + |F,mF>' is not valid")

        if not (I + J >= F >= abs(I - J) and diff.is_integer()):
            raise ParCaError("The input state of: 'I + J + |F,mF>' is not valid")

        if state_F[0].is_integer():
            if not (state_F[1].is_integer() and abs(state_F[1]) <= state_F[0]):
                raise ParCaError("The input state of: 'I + J + |F,mF>' is not valid")
        else:
            if not (state_F[2].is_integer() and state_F[3].is_integer() and abs(state_F[1]) <= state_F[0]):
                raise ParCaError("The input state of: 'I + J + |F,mF>' is not valid")

    @staticmethod
    def state_check_JI(I, J, mI, mJ):
        """
        Check if the given state is a physical valid one. If not: ParCaError appear.
        :param I: nuclear angular momentum [1]
        :param mI: z projection of the total angular momentum [1]
        :param mJ: z projection of the total electronic angular momentum [1]
        """
        state_I = np.array([1.0 * I, 1.0 * mI, 2.0 * I, 2.0 * mI])
        state_J = np.array([1.0 * J, 1.0 * mJ, 2.0 * J, 2.0 * mJ])

        if not (I >= 0 and state_I[2].is_integer() and J >= 0 and state_J[2].is_integer()):
            raise ParCaError("The input state of: '|I,mI,J,mJ>' is not valid")

        for state in [state_I, state_J]:
            if state[0].is_integer():
                if not (state[1].is_integer() and abs(state[1]) <= state[0]):
                    raise ParCaError("The input state of: '|I,mI,J,mJ>' is not valid")
            else:
                if not (state[2].is_integer() and state[3].is_integer() and abs(state[1]) <= state[0]):
                    raise ParCaError("The input state of: '|I,mI,J,mJ>' is not valid")

    @staticmethod
    def hf_lande(F, I, J, gJ, gI):
        """
        calculate the hyperfine Lande g-factor. (23) in {Cesium D Line Data}
        :param F: total angular momentum [1]
        :param I: nuclear angular momentum [1]
        :param J: total electronic angular momentum [1]
        :param gJ: Lande-factor (g-factor) of the fine structure [1]
        :param gI: Lande-factor (g-factor) of the nuclear [1]
        :return: hyperfine Lande-factor (g-factor) [1]
        """
        gF = 1 / (2 * F * (F+1)) \
             * (gJ * (F * (F+1) - I * (I+1) + J * (J+1)) + gI * (F * (F+1) + I * (I+1) - J * (J+1)))
        return gF

class BreitRabi:
    @staticmethod
    def state_mapping_F_to_IJ(I, F, mF):
        """
        State mapping in the sense of the Breit-Rabi shift.
        Please note that this formula only is valid in case: L=0, S=1/2 => J=1/2
        A representation of the hyperfine state with |F,mF> make only sense for 0 magnetic field or - in first order
        perturbation theory - for small magnetic fields (where the magnetic field splitting is small compared to the
        hyperfine splitting). In some textbooks labeled as the Zeeman regime.
        A representation via |I,mI;J,mJ> make  only sense for very strong fields (compared to the hyperfine splitting).
        In some textbooks labeled as Paschen–Back regime.
        However, since energy levels which can be calculated via the Breit-Rabi formula are non-crossing
        continuous functions, it is possible to map the labeling of the line on the right hand side of the diagram
        (low B) to the labeling of the same line on the left hand side (high B).
        But please note, beyond the corresponding regimes the labeling |F,mF> or |I,mI;J,mJ> have no physical meaning.
        :param I: nuclear angular momentum [1]
        :param F: total angular momentum [1]
        :param mF: z projection of the total angular momentum [1]
        :return: np array with two floats;
                mI: z projection of the nuclear angular momentum [1]
                mJ: z projection of the total electronic angular momentum [1]
        """
        # check if the input is valid
        HyperfineFcts.state_check_F(I, 1 / 2, F, mF)  # J = 1/2

        if mF == I + 1/2:
            mJ = 1/2

        elif mF == - (I + 1/2):
            mJ = -1/2

        elif F == I + 1/2:
            mJ = 1/2

        elif F == I - 1/2:
            mJ = -1/2

        mI = mF - mJ

        return np.array([mI, mJ])

    @staticmethod
    def state_mapping_IJ_to_F(I, mI, mJ):
        """
        See state_mapping_F_to_JI
        :param I: nuclear angular momentum [1]
        :param mI: z projection of the total angular momentum [1]
        :param mJ: z projection of the total electronic angular momentum [1]
        :return: np array with two floats;
                F: total angular momentum [1]
                mF: z projection of the total angular momentum [1]
        """
        # check if the input is valid
        HyperfineFcts.state_check_JI(I, 1 / 2, mI, mJ)  # J = 1/2

        mF = mI + mJ

        if mF == I + 1 / 2:
            F = I + 1 / 2

        elif mF == -(I + 1 / 2):
            F = I + 1 / 2

        elif mJ == 1 / 2:
            F = I + 1 / 2

        elif mJ == -1 / 2:
            F = I - 1 / 2

        return np.array([F, mF])

    @staticmethod
    def hbf_splitting(A, I):
        """
        calculate the hyperfine splitting out of the magnetic dipole constant
        {Cesium D Line Data}
        Note: I + 1/2 = F in this case
        :param A: magnetic dipole constant, ordinary frequency [J]
        :param I: nuclear spin [1]
        :return: hyperfine splitting, ordinary frequency [J]
        """
        return A * (I + 1/2)

    @staticmethod
    def breit_rabi(F, I, mF, deltaEhfs, gJ, gI, B):
        """
        The Breit-Rabi-Formula is valid for alkali metals (H-atom like) where the valence electron need to be in the
        s-shell (L=0, since just one e- in the outer shell => S=1/2, thus J=1/2).
        It calculated the hyperfine splitting with external magnetic field. Where the field strange can be low (Zeeman-
        Regime), very high (Paschen-Back-Regime) or intermediate.
        To understand the labeling of the states, please read the comment of state_mapping_F_to_JI!
        Formula/definitions like in {Cesium D Line Data}, for some background understanding: {wiki Breit-Rabi-Formel}
        :param F: total angular momentum [1]
        :param I: nuclear angular momentum [1]
        :param mF: z projection of the total angular momentum [1]
        :param deltaEhfs: hyperfine splitting [J]
        :param gJ: Lande-factor (g-factor) of the fine structure [1]
        :param gI: Lande-factor (g-factor) of the nuclear [1]
        :param B: magnetic field in Tesla [T]
        :return: energy of the state |F,mf> in Joule [J]
        """
        # check if the input is valid
        HyperfineFcts.state_check_F(I, 1/2, F, mF)  # J = 1/2

        if abs(mF) == I + 1/2:
            sign = mF/abs(mF)  # the sign is given by the sign of mF

            energy = deltaEhfs * I / (2 * I + 1) \
                     + sign * 1/2 * (gJ + 2 * I * gI) * cd['uB'][0] * B

        else:
            # in this case the sign is given either by the sign of mJ or by F
            if F == I + 1/2:
                sign = +1
            elif F == I - 1/2:
                sign = -1

            x = (gJ - gI) * cd['uB'][0] * B / deltaEhfs

            energy = - deltaEhfs / (2 * (2 * I + 1)) \
                     + gI * mF * cd['uB'][0] * B \
                     + sign * deltaEhfs / 2 * np.sqrt(1 + 4 * mF * x / (2 * I + 1) + x**2)

        return energy


class StateEnergy:
    """
    Class to calculate the energy of hyperfine state in an alkali atom (in B-Field) using the Breit-Rabi equation. Here
    the specific values are inserted.
    """
    def __init__(self, species, state_kind, state):
        """
        Init al necessary spect to calculate the energy. Without the magnetic field
        :param species: char, 'Cs' or '6Li'
        :param state_kind: char, 'F,mF' or 'mI,mJ' (the way how to label a certain state)
        :param state: np.array with two numbers
        """
        self.J = 1 / 2  # as always for Breit-Rabi

        if species == 'Cs':
            self.I = cd['I_Cs'][0]  # [1]
            self.gI = cd['gI_Cs'][0]  # [1]
            self.gJ = cd['gJ_S12_Cs'][0]  # Breit-Rabi just for ground state valid S12 {Cesium D Line Data}
            self.deltaEhfs = BreitRabi.hbf_splitting(cd['A_S12_Cs'][0], self.I)  # [J]
        elif species == '6Li':
            self.I = cd['I_6Li'][0]  # [1]
            self.gI = cd['gI_6Li'][0]  # [1]
            self.gJ = cd['gJ_S12_6Li'][0]  # Breit-Rabi just for ground state valid S12 {Cesium D Line Data}, [1]
            self.deltaEhfs = BreitRabi.hbf_splitting(cd['A_S12_6Li'][0], self.I)  # [J]
        else:
            raise ParCaError("Species not implemented")

        if state_kind == "F,mF":
            self.F = state[0]
            self.mF = state[1]
            HyperfineFcts.state_check_F(self.I, self.J, self.F, self.mF)
        elif state_kind == "mI,mJ":
            self.mI = state[0]
            self.mJ = state[1]
            self.F, self.mF = BreitRabi.state_mapping_IJ_to_F(self.I, self.mI, self.mJ)
        else:
            raise ParCaError("Invalid kind of state")

    def get_state_F(self):
        return self.F, self.mF

    def get_state_IJ(self):
        try:
            self.mI
        except:
            self.mI, self.mJ = BreitRabi.state_mapping_F_to_IJ(self.I, self.F, self.mF)
        return self.mI, self.mJ

    def get_energy(self, B):
        """
        :param B: magnetic field in Tesla [T]
        :return: energy in Joule [J]
        """
        return BreitRabi.breit_rabi(self.F, self.I, self.mF, self.deltaEhfs, self.gJ, self.gI, B)

class PolarizationFcts:
    """
    Collection of functions to calculate the polarization, connection potential, scattering rate
    """
    @ staticmethod
    def real_alpha(d1_line_strength, d2_line_strength, linewidth_d1, linewidth_d2, omega_d1, omega_d2, omega):
        """
        In case of a large detuning (in comparison with the hyperfine splitting) the real part of the atomic
        polarizability can be calculated with the following formula (fine in the LiCs experiment).
        Also it most yield that the scattering rate is much smaller then the damping rate (very low saturation). This is
        also fine for typical optical dipole traps.
        {Grimm-Optical Traps}
        :param d1_line_strength, d2_line_strength: line strength factors, dependent of the polarization state of
                                                    light for the D1 and D2 transition respectively [1]
        :param linewidth_d1, linewidth_d2: Natural Line Width (FWHM)/decay rate of D2 transition. This can be also
                                            understood as a damping rate. Angular frequency [1/s]
        :param omega_d1, omega_d2: angular frequency of the D1 and D2 transition respectively [1/s]
        :param omega: probe angular frequency [1/s]
        :return: real part of the polarizability alpha [A*s/(V*m)]
        """
        d1_part = d1_line_strength * linewidth_d1 / omega_d1 ** 3 \
                  * (1 / (omega_d1 - omega)
                     + 1 / (omega_d1 + omega))

        d2_part = d2_line_strength * linewidth_d2 / omega_d2 ** 3 \
                  * (1 / (omega_d2 - omega)
                     + 1 / (omega_d2 + omega))

        return 3 * np.pi * cd['epsilon_0'][0] * cd['c'][0]**3 * (d1_part + d2_part)

    @staticmethod
    def imag_alpha(d1_line_strength, d2_line_strength, linewidth_d1, linewidth_d2, omega_d1, omega_d2, omega):
        """
        Calculate the imaginary part of the polarizability. More info: see description real_alpha.
        {Grimm-Optical Traps}
        :param d1_line_strength, d2_line_strength: line strength factors, dependent of the polarization state of
                                                    light for the D1 and D2 transition respectively [1]
        :param linewidth_d1, linewidth_d2: Natural Line Width (FWHM)/decay rate of D2 transition. This can be also
                                            understood as a damping rate. Angular frequency [1/s]
        :param omega_d1, omega_d2: angular frequency of the D1 and D2 transition respectively [1/s]
        :param omega: probe angular frequency [1/s]
        :return: imaginary part of the polarizability alpha [A*s/(V*m)]
        """
        d1_part = d1_line_strength * linewidth_d1**2 / omega_d1**6 * (1 / (omega_d1 - omega) + 1 / (omega_d1 + omega))**2

        d2_part = d2_line_strength * linewidth_d2**2 / omega_d2**6 * (1 / (omega_d2 - omega) + 1 / (omega_d2 + omega))**2

        return 3 * np.pi * cd['epsilon_0'][0] * cd['c'][0]**3 * omega**3 / 2 * (d1_part + d2_part)

    @staticmethod
    def dipole_potential(d1_line_strength, d2_line_strength, linewidth_d1, linewidth_d2, omega_d1, omega_d2, omega, I):
        """
        Calculate the optical dipole potential out of the polarizability of the considered atoms.
        Just valid as long the optical detuning stay large in comparison with the excited-state hyperfine splitting.
        {Grimm-Optical Traps}
        :param d1_line_strength, d2_line_strength: line strength factors, dependent of the polarization state of
                                                    light for the D1 and D2 transition respectively [1]
        :param linewidth_d1, linewidth_d2: Natural Line Width (FWHM)/decay rate of D2 transition. This can be also
                                            understood as a damping rate. Angular frequency [1/s]
        :param omega_d1, omega_d2: angular frequency of the D1 and D2 transition respectively [1/s]
        :param omega: probe angular frequency [1/s]
        :param I: light intensity [W/m^2]
        :return: dipole potential [J]
        """
        real_alpha = PolarizationFcts.real_alpha(d1_line_strength, d2_line_strength, linewidth_d1, linewidth_d2,
                                                 omega_d1, omega_d2, omega)
        return - 1 / (2 * cd['epsilon_0'][0] * cd['c'][0]) * real_alpha * I

    @staticmethod
    def scattering_rate(d1_line_strength, d2_line_strength, linewidth_d1, linewidth_d2, omega_d1, omega_d2, omega, I):
        """
        "Considering the light as a steam of photons, the absorption [imaginary part of alpha] can be interpreted in
        terms of photon scattering"
        {Grimm-Optical Traps}
        :param d1_line_strength, d2_line_strength: line strength factors, dependent of the polarization state of
                                                    light for the D1 and D2 transition respectively [1]
        :param linewidth_d1, linewidth_d2: Natural Line Width (FWHM)/decay rate of D2 transition. This can be also
                                            understood as a damping rate. Angular frequency [1/s]
        :param omega_d1, omega_d2: angular frequency of the D1 and D2 transition respectively [1/s]
        :param omega: probe angular frequency [1/s]
        :param I: light intensity [W/m^2]
        :return: scattering rate [1/s]
        """
        imag_alpha = PolarizationFcts.imag_alpha(d1_line_strength, d2_line_strength, linewidth_d1, linewidth_d2,
                                                 omega_d1, omega_d2, omega)
        return 1 / (cd['hbar'][0] * cd['epsilon_0'][0] * cd['c'][0]) * imag_alpha * I

    @staticmethod
    def line_strength_factors(P, F=None, mF=None, I=None, gJ=None, gI=None):
        """
        calculate the line strength factors for alkali atoms for the D1 and D2 line. For linear polarized light, the
        situation is easy, for circular a bit more involved and requires the specify the ground state of the transition.
        For this ground state F, mF, I, gJ, and gI is necessary. Please note that gJ is state dependent (fine structure
        state)
        (For Cs typically |F=3,mF=3>, for Li typically |F=1/2,mF=1/2> or |F=1/2,mF=-1/2>)
        {Grimm-Optical Traps}
        :param P: laser polarization:
                    0-> linear, 1 -> sigma^+, -1 -> sigma^-
        :param F: total angular momentum [1]
        :param mF: z projection of the total angular momentum [1]
        :param I: nuclear angular momentum [1]
        :param gJ: Lande-factor (g-factor) of the fine structure [1]
        :param gI: Lande-factor (g-factor) of the nuclear [1]
        :return: array: [d1_line_strength, d2_line_strength] [1]
        """
        if P == 0:
            d1_line_strength = 1/3
            d2_line_strength = 2/3
        else:
            if F is None or mF is None or I is None or gJ is None or gI is None:
                raise ParCaError("Calculation with circular polarization need the ground state to calculate the line strength factors")

            J = 1/2  # the expression is only valid for alkali atoms

            HyperfineFcts.state_check_F(I, J, F, mF)

            gF = HyperfineFcts.hf_lande(F, I, J, gJ, gI)

            if P == 1:
                d1_line_strength = 1/3 * (1 - gF * mF)
                d2_line_strength = 1/3 * (2 + gF * mF)
            elif P == -1:
                d1_line_strength = 1/3 * (1 + gF * mF)
                d2_line_strength = 1/3 * (2 - gF * mF)
            else:
                raise ParCaError("Invalid value for P")

        return [d1_line_strength, d2_line_strength]

class DipoleTrapCal:
    """
    Calculate the polarizability and its effect. Specific values for are inserted here.
    """
    def __init__(self, species, P, F=None, mF=None):
        """
        Define a certain situation to give values for the dipole trap potential
        :param species: char, 'Cs' or '6Li' is implemented
        :param P: laser polarization:
                    0-> linear, 1 -> sigma^+, -1 -> sigma^-
        :param F: total angular momentum [1]
        :param mF: z projection of the total angular momentum [1]
        :param I: nuclear angular momentum [1]
        """
        if species == 'Cs':
            self.linewidth_d1 = cd['linewidthD1_Cs'][0]
            self.omega_d1 = cd['omegaD1_Cs'][0]
            self.linewidth_d2 = cd['linewidthD2_Cs'][0]
            self.omega_d2 = cd['omegaD2_Cs'][0]
            self.gJ = cd['gJ_S12_Cs'][0]  # to calculate the line strength factors for circular polarized light the
                                            # properties of the ground state are necessary
            self.gI = cd['gI_Cs'][0]
            self.I = cd['I_Cs'][0]

        elif species == '6Li':
            self.linewidth_d1 = cd['linewidthD1_6Li'][0]
            self.omega_d1 = cd['omegaD1_6Li'][0]
            self.linewidth_d2 = cd['linewidthD2_6Li'][0]
            self.omega_d2 = cd['omegaD2_6Li'][0]
            self.gJ = cd['gJ_S12_6Li'][0]  # to calculate the line strength factors for circular polarized light the
                                            # properties of the ground state are necessary
            self.gI = cd['gI_6Li'][0]
            self.I = cd['I_6Li'][0]
        else:
            raise ParCaError("Species not implemented")

        self.d1_line_strength, self.d2_line_strength = PolarizationFcts.line_strength_factors(P, F, mF, self.I, self.gJ,
                                                                                              self.gI)

    def get_real_alpha(self, omega):
        """
        Calculate real alpha for the situation defined with init
        :param omega: probe angular frequency [1/s]
        :return: real part of the polarizability alpha [A*s/(V*m)]
        """
        re_alpha = PolarizationFcts.real_alpha(self.d1_line_strength, self.d2_line_strength, self.linewidth_d1,
                                               self.linewidth_d2, self.omega_d1, self.omega_d2, omega)
        return re_alpha

    def get_imag_alpha(self, omega):
        """
        Calculate imag alpha for the situation defined with init
        :param omega: probe angular frequency [1/s]
        :return: imaginary part of the polarizability alpha [A*s/(V*m)]
        """
        imag_alpha = PolarizationFcts.imag_alpha(self.d1_line_strength, self.d2_line_strength, self.linewidth_d1,
                                                 self.linewidth_d2, self.omega_d1, self.omega_d2, omega)
        return imag_alpha

    def get_dipole_potential(self, omega, I):
        """
        Calculate dipole potential for the situation defined with init
        :param omega: probe angular frequency [1/s]
        :param I: light intensity [W/m^2]
        :return: dipole potential [J]
        """
        pot = PolarizationFcts.dipole_potential(self.d1_line_strength, self.d2_line_strength, self.linewidth_d1,
                                                self.linewidth_d2, self.omega_d1, self.omega_d2, omega, I)
        return pot

    def get_scattering_rate(self, omega, I):
        """
        Calculate scattering rate for the situation defined with init
        :param omega: probe angular frequency [1/s]
        :param I: light intensity [W/m^2]
        :return: scattering rate [1/s]
        """
        sct = PolarizationFcts.scattering_rate(self.d1_line_strength, self.d2_line_strength, self.linewidth_d1,
                                               self.linewidth_d2, self.omega_d1, self.omega_d2, omega, I)
        return sct

