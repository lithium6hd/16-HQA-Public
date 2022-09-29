from scipy.constants import k, c, h, hbar, g, epsilon_0, physical_constants, tera, giga, mega, kilo, centi, milli, \
    micro, nano
from numpy import pi  # , sqrt, exp, log
import warnings

from Calculations.ParCa.Constants.transfer import TransferFcts
tF = TransferFcts

class ConstantsWarning():
    """
    create warning message
    These warnings should rise up if things are not as precise as one expect.
    For instance if old values of experimental data ist used or if the calculation method itself is not as precise as
    one could thing
    """

    def __init__(self, message):
        message = "Constants WARNING: %s" % str(message)
        warnings.warn(message)


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

    const_dict['gS'] = [2.0023193043622, '1', 'electron spin g-factor {Cesium D Line Data}']
    const_dict['gL'] = [0.99999587, '1', 'electron orbital g-factor {Cesium D Line Data}']

    # Cesium
    const_dict['m_Cs'] = [132.905451931 * const_dict['u0'][0], 'kg', 'mass of 133Cs {Cesium D Line Data}']
    const_dict['I_Cs'] = [7 / 2, '1', 'nuclear spin of 133Cs {Cesium D Line Data}']
    const_dict['gI_Cs'] = [-0.00039885395, '1', 'nuclear g-factor {Cesium D Line Data}']
    const_dict['A_S12_Cs'] = [h * 2.2981579425 * giga, 'J',
                              'magnetic dipole constant of the 6 2S_1/2 state of 133Cs {Cesium D Line Data}']
    const_dict['A_P12_Cs'] = [h * 291.9201 * mega, 'J',
                              'magnetic dipole constant of the 6 2P_1/2 state of 133Cs {Cesium D Line Data}']
    const_dict['A_P32_Cs'] = [h * 50.28827 * mega, 'J',
                              'magnetic dipole constant of the 6 2P_3/2 state of 133Cs {Cesium D Line Data}']
    const_dict['gJ_S12_Cs'] = [2.00254032, '1',
                               'fine structure Lande g-factor of the 6 2S_1/2 state of 133Cs {Cesium D Line Data}']
    const_dict['gJ_P12_Cs'] = [0.665900, '1',
                               'fine structure Lande g-factor of the 6 2P_1/2 state of 133Cs {Cesium D Line Data}']
    const_dict['gJ_P32_Cs'] = [1.33400, '1',
                               'fine structure Lande g-factor of the 6 2P_3/2 state of 133Cs {Cesium D Line Data}']

    # Cs D1 transition (6 2S_1/2 -> 6 2P_1/2)
    const_dict['omegaD1_Cs'] = [2 * pi * 335.116048807 * tera, '1/s',
                                'angular frequency of D1 transition of 133Cs {Cesium D Line Data}']
    const_dict['linewidthD1_Cs'] = [2 * pi * 4.575 * mega, '1/s',
                                    'Natural Line Width (FWHM)/decay rate of D1 transition in angular frequency (Cs) {Cesium D Line Data}']

    # Cs D2 transition (6 2S_1/2 -> 6 2P_3/2)
    const_dict['omegaD2_Cs'] = [2 * pi * 351.72571850 * tera, '1/s',
                                'angular frequency of D2 transition of 133Cs {Cesium D Line Data}']
    const_dict['linewidthD2_Cs'] = [2 * pi * 5.234 * mega, '1/s',
                                    'Natural Line Width (FWHM)/decay rate of D2 transition in angular frequency (Cs) {Cesium D Line Data}']

    # Lithium
    const_dict['m_6Li'] = [6.0151214 * const_dict['u0'][0], 'kg', 'mass of 6Li {Properties of 6Li}']
    const_dict['I_6Li'] = [1, '1', 'nuclear spin of 6Li {Properties of 6Li}']
    const_dict['gI_6Li'] = [-0.0004476540, '1', 'nuclear g-factor {Properties of 6Li}']
    const_dict['A_S12_6Li'] = [h * 152.1368407 * mega, 'J',
                               'magnetic dipole constant of the 2 2S_1/2 state of 6Li {Properties of 6Li}']
    const_dict['A_P12_6Li'] = [h * 17.386 * mega, '1/s',
                               'magnetic dipole constant of the 2 2P_1/2 state of 6Li {Properties of 6Li}']
    const_dict['A_P32_6Li'] = [h * (-1.155) * mega, '1/s',
                               'magnetic dipole constant of the 6 2P_3/2 state of 6Li {Properties of 6Li}']
    const_dict['gJ_S12_6Li'] = [2.0023010, '1',
                                'fine structure Lande g-factor of the 6 2S_1/2 state of 6Li {Properties of 6Li}']
    const_dict['gJ_P12_6Li'] = [0.6668, '1',
                                'fine structure Lande g-factor of the 6 2P_1/2 state of 6Li {Properties of 6Li}']
    const_dict['gJ_P32_6Li'] = [1.335, '1',
                                'fine structure Lande g-factor of the 6 2P_3/2 state of 6Li {Properties of 6Li}']

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


def experimental_parameter():
    """
    Collection of experimental parameters. To use this collection use get_experimental_parameter
    This collection should also include experimental uncertainties (the 3rd position of the dict is in some cases deltaX
    if the full expression is: x +- deltaX. It is important that one comment on the
    uncertainties, how they are calculated
    :return: dict
    """
    expara_dict = dict()

    # Reservoir Trap
    expara_dict['RT_pol_power_red'] = [0.8575, '1', "Reservoir Trap: Reduction factor due to Polarization.\n"
                                                    "Further information can be found in\n"
                                                    "{LiCs_Pro/eLabNotes/A. Cesium/Reservoir Trap Characterization/"
                                                    "Polarisation in Cesium Optical potential}\n"
                                                    "notes taken on January 18, 2021"]
    expara_dict['RT_power_los_power_red'] = [0.85, '1', "Reservoir Trap: Reduction factor due to power loss of second "
                                                        "beam.\n"
                                                        "This value was used in Manuel's mathematica script\n"
                                                        "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['RT_waist_first_dz'] = [600 * micro, 'm',
                                        "Reservoir Trap: beam waist radius in dz direction of the "
                                        "first beam, see beam-profile-coordinates.png.\n"
                                        "Value was used in Manuel's mathematica script.\n"
                                        "ATTENTION: After Rebuild of the Reservoir Trap this has to "
                                        "be verified again\n"
                                        "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['RT_waist_first_d1'] = [600 * micro, 'm',
                                        "Reservoir Trap: beam waist radius in d1 direction of the "
                                        "first beam, see beam-profile-coordinates.png."
                                        "Value was used in Manuel's mathematica script.\n"
                                        "ATTENTION: After Rebuild of the Reservoir Trap this has to"
                                        " be verified again\n"
                                        "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['RT_waist_second_dz'] = [600 * micro, 'm',  "Reservoir Trap: beam waist radius in dz direction of the "
                                                            "second beam, see beam-profile-coordinates.png.\n"
                                                            "Value was used in Manuel's mathematica script\n"
                                                            "ATTENTION: After Rebuild of the Reservoir Trap this has "
                                                            "to be verified again\n"
                                                            "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['RT_waist_second_d1'] = [600 * micro, 'm',
                                         "Reservoir Trap: beam waist radius in d1 direction of the second beam, "
                                         "see beam-profile-coordinates.png. \n"
                                         "Value was used in Manuel's mathematica script.\n"
                                         "ATTENTION: After Rebuild of the Reservoir Trap this has to be verified "
                                         "again!\n"
                                         "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['RT_waist_first_rho'] = [0, '1', "Reservoir Trap: correlation coefficient of the first beam. Value"
                                                 "was used in Manuel's mathematica script.\n"
                                                 "ATTENTION: After Rebuild of the Reservoir Trap this has to be "
                                                 "verified again!\n"
                                                 "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['RT_waist_second_rho'] = [0, '1',
                                          "Reservoir Trap: correlation coefficient of the second beam. Value "
                                          "was used in Manuel's mathematica script.\n"
                                          "ATTENTION: After Rebuild of the Reservoir Trap this has to be"
                                          " verified again!\n"
                                          "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['RT_beam_angle'] = [270 + 45, 'degree', "Reservoir Trap: propagation direction of the first incoming"
                                                       "beam in degree. The angle is measured in the mathematical "
                                                       "positive direction from the y-axis on, see xyz-coordinates.png"]
    expara_dict['RT_crossing_angle'] = [90, 'degree', "Reservoir Trap: crossing angle in the xy plane in degree. "
                                                      "The angle is measured in the mathematical positive direction,"
                                                      "starting from the first incoming beam and end with the second"
                                                      "incoming beam"]
    expara_dict['RT_frequency'] = [tF.wlen_to_omega(1070 * nano), '1/s', "Reservoir Trap: angular frequency of the "
                                                                         "used laser"]

    # Dimple Trap
    expara_dict['DT_waist_first_dz'] = [62 * micro, 'm', "Dimple Trap: beam waist radius in dz direction of the "
                                        "first beam, see beam-profile-coordinates.png.\n"
                                        "Value was used in Manuel's mathematica script.\n"
                                        "ATTENTION: This should be done again!\n"
                                        "{LiCs_Pro/eLabNotes/B. Lithium/Dimple Trap Loading/Dimple-Trap Trapping "
                                        "frequency} written on March 12, 2020"
                                        "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['DT_waist_first_d1'] = [62 * micro, 'm',  "Dimple Trap: beam waist radius in d1 direction of the "
                                        "first beam, see beam-profile-coordinates.png.\n"
                                        "Value was used in Manuel's mathematica script.\n"
                                        "ATTENTION: This should be done again!\n"
                                        "{LiCs_Pro/eLabNotes/B. Lithium/Dimple Trap Loading/Dimple-Trap Trapping "
                                        "frequency} written on March 12, 2020"
                                        "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['DT_waist_second_dz'] = [62 * micro, 'm', "Dimple Trap: beam waist radius in dz direction of the "
                                         "second beam, see beam-profile-coordinates.png.\n"
                                         "Value was used in Manuel's mathematica script.\n"
                                         "ATTENTION: This should be done again!\n"
                                         "{LiCs_Pro/eLabNotes/B. Lithium/Dimple Trap Loading/Dimple-Trap Trapping "
                                         "frequency} written on March 12, 2020"
                                         "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['DT_waist_second_d1'] = [62 * micro, 'm',  "Dimple Trap: beam waist radius in d1 direction of the "
                                         "second beam, see beam-profile-coordinates.png.\n"
                                         "Value was used in Manuel's mathematica script.\n"
                                         "ATTENTION: This should be done again!\n"
                                         "{LiCs_Pro/eLabNotes/B. Lithium/Dimple Trap Loading/Dimple-Trap Trapping "
                                         "frequency} written on March 12, 2020"
                                         "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['DT_waist_first_rho'] = [0, '1',
                                         "Dimple Trap: correlation coefficient of the first beam. Value "
                                         "was used in Manuel's mathematica script.\n"
                                         "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['DT_waist_second_rho'] = [0, '1',
                                          "Dimple Trap: correlation coefficient of the second beam. Value "
                                          "was used in Manuel's mathematica script.\n"
                                          "{MasterScript_ManuelGerken_16062021.nb}"]
    expara_dict['DT_beam_angle'] = [90 - 4.2, 'degree', "Dimple Trap: propagation direction of the first incoming"
                                                         "beam in degree. The angle is measured in the mathematical "
                                                         "positive direction from the y-axis on, see xyz-coordinates"
                                                         ".png"]
    expara_dict['DT_crossing_angle'] = [180 + 8.3, 'degree', 0.9, "Dimple Trap: crossing angle in the xy plane in "
                                                                  "degree. The angle is measured in the mathematical "
                                                                  "positive direction, starting from the first "
                                                                  "incoming beam and end with the second incoming "
                                                                  "beam. Value was updated after insert the dichroic "
                                                                  "mirror.The error is an estimation on the measured "
                                                                  "lengths and by using gaussian uncertainty "
                                                                  "propagation the uncertainty above occur. "
                                                                  "(2022-01-19)"]
    expara_dict['DT_crossing_angle-old-1'] = [180 + 8.4, 'degree',  "Dimple Trap: crossing angle in the xy plane in "
                                                                    "degree. The angle is measured in the mathematical "
                                                                    "positive direction, starting from the first "
                                                                    "incoming beam and end with the second incoming "
                                                                    "beam. Value was taken before the dichroic mirror"
                                                                    "was inserted (2022-01-19)."]
    expara_dict['DT_frequency'] = [tF.wlen_to_omega(1070 * nano), '1/s', "Dimple Trap: angular frequency of the used "
                                                                         "laser"]

    # Alternative Dimple Trap
    expara_dict['ADT_waist_first_dz'] = [65 * micro, 'm', "Alternative Dimple Trap: beam waist radius in dz direction "
                                                          "of the first beam, see beam-profile-coordinates.png.\n"
                                                          "Value was used in Eleonoras's mathematica script.\n"
                                                          "{PolaronNumbers_EleonoraLippi_18082021.nb}"]
    expara_dict['ADT_waist_first_d1'] = [65 * micro, 'm', "Alternative Dimple Trap: beam waist radius in d1 direction "
                                                          "of the first beam, see beam-profile-coordinates.png.\n"
                                                          "Value was used in Eleonoras's mathematica script.\n"
                                                          "{PolaronNumbers_EleonoraLippi_18082021.nb}"]
    expara_dict['ADT_waist_second_dz'] = [65 * micro, 'm', "Alternative Dimple Trap: beam waist radius in dz direction "
                                                           "of the second beam, see beam-profile-coordinates.png.\n"
                                                           "Value was used in Eleonoras's mathematica script.\n"
                                                           "{PolaronNumbers_EleonoraLippi_18082021.nb}"]
    expara_dict['ADT_waist_second_d1'] = [65 * micro, 'm', "Alternative Dimple Trap: beam waist radius in d1 direction "
                                                           "of the second beam, see beam-profile-coordinates.png.\n"
                                                           "Value was used in Eleonoras's mathematica script.\n"
                                                           "{PolaronNumbers_EleonoraLippi_18082021.nb}"]
    expara_dict['ADT_waist_first_rho'] = [0, '1',
                                          "Alternative Dimple Trap: correlation coefficient of the first beam."
                                          "Value was used in Eleonoras's mathematica script.\n"
                                          "{PolaronNumbers_EleonoraLippi_18082021.nb}"]
    expara_dict['ADT_waist_second_rho'] = [0, '1',
                                           "Alternative Dimple Trap: correlation coefficient of the second beam."
                                           "Value was used in Eleonoras's mathematica script.\n"
                                           "{PolaronNumbers_EleonoraLippi_18082021.nb}"]
    expara_dict['ADT_beam_angle'] = [90 + 45 + 4, 'degree', "Dimple Trap: propagation direction of the first incoming"
                                                         "beam in degree. The angle is measured in the mathematical "
                                                         "positive direction from the y-axis on, see xyz-coordinates"
                                                         ".png\n"
                                                         "Value was used in Eleonoras's mathematica script.\n"
                                                         "{PolaronNumbers_EleonoraLippi_18082021.nb}"]
    expara_dict['ADT_crossing_angle'] = [90, 'degree', "Alternative Dimple Trap: crossing angle in the xy plane in "
                                                       "degree. The angle is measured in the mathematical positive "
                                                       "direction, starting from the first incoming beam and end "
                                                       "with the second incoming beam"]
    expara_dict['ADT_frequency'] = [tF.wlen_to_omega(1064 * nano), '1/s', "Alternative Dimple Trap: angular frequency "
                                                                          "of the used laser"]
    # Micro Trap
    expara_dict['MT_waist_dz'] = [10.5 * micro, 'm', "Microtrap: beam waist in the dz direction.\n"
                                                     "ATTENTION: This value is just an estimation!\n"
                                                     "see beam-profile-coordinates.png\n"
                                                     "Value was used in Eleonoras's mathematica script.\n"
                                                     "{PolaronNumbers_EleonoraLippi_18082021.nb}"]
    expara_dict['MT_waist_d1'] = [10.5 * micro, 'm', "Microtrap: beam waist in the d1 direction.\n"
                                                     "ATTENTION: This value is just an estimation!\n"
                                                     "see beam-profile-coordinates.png\n"
                                                     "Value was used in Eleonoras's mathematica script.\n"
                                                     "{PolaronNumbers_EleonoraLippi_18082021.nb}"]
    expara_dict['MT_waist_rho'] = [0, '1', "Mixrotrap: correlation coefficient\n"
                                           "Value was used in Eleonoras's mathematica script.\n"
                                           "{PolaronNumbers_EleonoraLippi_18082021.nb}"]
    expara_dict['MT_beam_angle'] = [90, 'degree', "Microtrap: propagation direction of the beam in degree. The angle"
                                                   "is measured in the mathematical positive direction from the y-axis"
                                                   "on, see xyz-coordinates.png\n"
                                                   "Value was used in Eleonoras's mathematica scrict.\n"
                                                   "{PolaronNumbers_EleonoraLippi_18082021.nb}"]
    expara_dict['MT_frequency'] = [tF.wlen_to_omega(880.1821249976 * nano), '1/s', "Microtrap: angular frequency of "
                                                                                   "the used laser"]

    # Feshbach coils
    expara_dict['FC_alpha'] = [274, '1/m^2', 6, "Feshbach coils: the Helmholtz configuration of the feshbach coils is"
                                                "far from perfect, this lead to a certain quadratic contriubtion.\n"
                                                "The uncertainty is just due to the fit, we expect much higher"
                                                "uncertainties in the real experiment."
                                                "{OneNote/GroupNotes/07 Coding/ParCa/Old Mathematica corrections} "
                                                "and {PhD-Gerken}"]

    # curvature coils
    expara_dict['CC_x0'] = [-1.2 * milli, 'm', 0.1 * milli, "Curvature coils: The center of the curvature coils is not"
                                                            "equal to the center of the Feshbach coils, which define"
                                                            "the center of the whole system\n"
                                                            "{OneNote/GroupNotes/07 Coding/ParCa/Additional Information"
                                                            "on the Equations in the code}"]

    expara_dict['CC_y0'] = [-1.2 * milli, 'm', 0.1 * milli, "Curvature coils: The center of the curvature coils is not"
                                                            "equal to the center of the Feshbach coils, which define"
                                                            "the center of the whole system\n"
                                                            "{OneNote/GroupNotes/07 Coding/ParCa/Additional Information"
                                                            "on the Equations in the code}"]

    return expara_dict


def get_experimental_parameter(name):
    """
    Function to pass the dictionary information and raises a Warning if necessary, warnings are though for values which
    are not as precise as they should be
    :param name: name of the variable
    :return: dic entry
    """
    if name == 'DT_waist_first_dz':
        ConstantsWarning('VALUE IS OLD: %s. According to {MasterScript_ManuelGerken_16062021.nb} this value'
                         ' should be measured again' % str(name))
    elif name == 'DT_waist_first_d1':
        ConstantsWarning('VALUE IS OLD: %s. According to {MasterScript_ManuelGerken_16062021.nb} this value'
                         ' should be measured again' % str(name))
    elif name == 'DT_waist_second_dz':
        ConstantsWarning('VALUE IS OLD: %s. According to {MasterScript_ManuelGerken_16062021.nb} this value'
                         ' should be measured again' % str(name))
    elif name == 'DT_waist_second_dz':
        ConstantsWarning('VALUE IS OLD: %s. According to {MasterScript_ManuelGerken_16062021.nb} this value'
                         ' should be measured again' % str(name))
    elif name == "MT_waist_dz":
        ConstantsWarning('VALUE IS AN ESTIMATE: %s' % str(name))
    elif name == "MT_waist_d1":
        ConstantsWarning('VALUE IS AN ESTIMATE: %s' % str(name))

    expara_dict = experimental_parameter()
    return expara_dict[name]
