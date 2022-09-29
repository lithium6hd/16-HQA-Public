import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from numpy import array as arr
from numpy import pi
from numpy import sqrt
from scipy import integrate
from scipy.special import erf
import fundamentals as fundi
from constants import load_constants, get_experimental_parameter, ConstantsWarning
from Formulas import ThermalFacts

uc_colors = {'red': '#C61A27', 'blue': '#3891A6', 'orange': '#F79D65', 'yellow': '#FDE74C', 'lightblue': '#C1FFF2',
             'green': '#9BC53D', 'white': '#FFFFFF', 'black': '#000000', 'grey': '#999999', 'lightgrey': 'lightgrey'}

# ------------------------------- Iteration Dependency?  --------------------------------

# Set False for a single run with iterations specified below. Set True for testing the iteration dependency
iteration_dependency = True

# ------------------------------- Experimental parameters  --------------------------------

# Laser
wavelength = 532e-9  # [m] Laser wavelength
P = 3  # [W] Laser power
omega = 2 * np.pi * load_constants()["c"][0] / wavelength  # [Rad] probe angular frequency
print("Laser P = {} W at {:.0f} nm".format(P, wavelength * 1e9))

# Atoms
nat_linewidth = load_constants()["linewidthD1_6Li"][0]  # [Hz]
I_sat = 75.9e-3  # [W/m^2] saturation intensity
mass = 9.9883414 * 1e-27  # [kg] atomic mass

# Adiabatic approximation:
kappa = 5 / 3  # approx. via d.o.f. (monoatomic gas)

# %% Set the steps for single or multiple iterations:

if iteration_dependency:
    steps_list = np.logspace(0.4, 3.7, num=650, dtype=int, endpoint=True)
else:
    steps_list = [5000]

print("Steps:", steps_list)

# %%
final_scatteringrates = []
final_totaltimes = []
final_atomnumbers = []
final_densities = []
final_temperatures = []

# ------------------------------- Simulation  --------------------------------

for steps in steps_list:

    if iteration_dependency:
        print("--------------New Iteration Number", steps, "-----------------")

    # ------------------------------- Initial facts  --------------------------------

    r = arr([])  # [m] Inner ring radius
    V = arr([])  # [m^3] Trap volume
    A = arr([])  # [m^2] Ring area
    Intensity = arr([])  # [W/m^2] Laser intensity in ring
    T = arr([])  # [K] Temperature before evaporation
    T_adiab = arr([])  # [K] Iteration intermediate step: temperature due to compression, used for truncated boltzmann
    U = arr([])  # Trap depth in [J]
    U_inK = arr([])  # Trap depth in [K]
    v_U = arr([])  # Trap depth in [m/s]
    n = arr([])  # density
    Number_of_atoms = arr([])
    ScatteringRate = arr([])
    ScatteringTime = arr([])

    # MOT density and temperature
    n = np.append(n, [1e17])  # [m^-3] initial density -> particle per µm^3 from MOT
    T = np.append(T, [1e-3])  # [K] Temperature from MOT before compression
    T_adiab = np.append(T_adiab, [1e-3])  # [K] Temperature due to compression, used only for calculations

    # Correct number of steps for correct output:
    steps = steps - 1

    # Ring facts:
    r_initial_set = 50e-6  # SET INITIAL RADIUS HERE [m]
    r_final = 5e-6  # [m] final radius
    delta_r = (r_initial_set - r_final) / steps
    d = 5e-6  # [m] ring radial thickness

    # Output initial radius
    r_initial = r_final + steps * delta_r
    r = np.append(r, [r_initial])
    print("Initial ring radius: {:.0f} µm".format(r[0] * 1e6))
    print("Final ring radius: {:.0f} µm".format(r_final * 1e6))
    print("Ring thickness: {:.0f} µm".format(d * 1e6))
    print("Step Size: {:.2e}".format(delta_r))

    # Z-confinement thickness
    z_thickness = 20e-6  # [m] size of the vertical space confined by the z-sheets

    # Iteration variable:
    i = 0

    # ------------------------------- LOOP START  --------------------------------

    while i <= steps:
        # ------------------------------- Ring calculations  --------------------------------

        # Calculate volume of confinement
        V = np.append(V, [np.pi * r[i] ** 2 * z_thickness])  # [m^3] ->  r=r1 cylindrical ODT

        # Distribute laser power over the initial area:
        A = np.append(A, [np.pi * ((r[i] + d) ** 2 - r[i] ** 2)])  # [m^2]  -> outer radius: 215µm
        Intensity = np.append(Intensity, [P / A[i]])  # Intensity in the initial ring

        # ------------------------------- Calculate temperature after compression  --------------------------------
        if i > 0:
            T_adiab = np.append(T_adiab,
                                [T[i - 1] * ((V[i - 1] / V[i]) ** (kappa - 1))])  # [K] adiabatic approximation?

        # ------------------------------- Calculate Potentials/Trap depths --------------------------------

        # Calculate potential resulting from the achieved intensity in ring areas:
        U = np.append(U, [fundi.PolarizationFcts.dipole_potential(1, 0,
                                                                  load_constants()["linewidthD1_6Li"][0],
                                                                  load_constants()["linewidthD2_6Li"][0],
                                                                  load_constants()['omegaD1_6Li'][0],
                                                                  load_constants()['omegaD2_6Li'][0],
                                                                  omega, Intensity[i])])  # [J]

        # Conversion of trap potential to temperature:
        U_inK = np.append(U_inK, [U[i] / load_constants()["kB"][0]])  # [K]

        # Conversion of trap depth into [m/s]:
        v_U = np.append(v_U, [np.sqrt(2 * U[i] / mass)])  # -> Cut in boltzmann due to ring potential

        # ----------------------------- Calculate Boltzmann distributions and atom density/numbers---------------------
        """
        Process (approximated):
        Boltzmann_i cut at U[i] -> thermalization -> Boltzmann'_i forms and yields T[i]
        """

        # 1st evaporation: Particle density after loss due to initial ring:
        if i > 0:
            n = np.append(n,
                          [integrate.quad(lambda v: ThermalFacts.boltzmann(v, T_adiab[i], mass), 0, v_U[i])[0] *
                           n[i - 1] * V[i - 1] / V[i]])  # Integral * the density (adjusted for volume decrease)

            # n = np.append(n,
            #               [integrate.quad(lambda e: ThermalFacts.boltzmann_energy(e, T_adiab[i]), 0, U[i])[0] *
            #                n[i - 1] * V[i - 1] / V[i]])  # Integral * the density (adjusted for volume decrease)

            # Calculate new T[i] from new boltzmann thermalization
            root = sqrt(U_inK[i] / T_adiab[i])
            T = np.append(T,
                          T_adiab[i] / 3 * (
                                  (3 * sqrt(pi) * erf(root) - (6 * root + 4 * (U_inK[i] / T_adiab[i]) ** (3 / 2)) *
                                   np.exp(-U_inK[i] / T_adiab[i]))
                                  /
                                  (sqrt(pi) * erf(root) - 2 * root * np.exp(-U_inK[i] / T_adiab[i])))
                          )

        # Prevent increase of radius in the last iteration:
        if i < steps:
            r = np.append(r, [r[i] - delta_r])

        Number_of_atoms = np.append(Number_of_atoms, n[i] * V[i])

        # Integer incremation
        i += 1

    # Scattering Rate Calculation
    ScatteringRate = ThermalFacts.scattering_rate(T, mass, n)
    ScatteringTime = 2.8 / ScatteringRate

    # Output
    print("Number of steps: ", i)

    print("Final radius: {:.2e} m".format(r[-1]))
    print("Final particle density: {:.2f} µm^-3".format(n[-1] * 1e-18))
    print("Final number of atoms: {:.0f}".format(Number_of_atoms[-1]))
    print("Final scattering rate: {:.2e} Hz".format(ScatteringRate[-1]))
    print("Final temperature: {:.2e} K".format(T[-1]))
    print("Final Trap Depth: {:.2e} K".format(U_inK[-1]))
    print("Summed Thermalization Time: {:.2e} s".format(np.sum(ScatteringTime)))

    final_atomnumbers = np.append(final_atomnumbers, Number_of_atoms[-1])
    final_densities = np.append(final_densities, n[-1] * 1e-18)
    final_temperatures = np.append(final_temperatures, T[-1])
    final_scatteringrates = np.append(final_scatteringrates, ScatteringRate[-1])
    final_totaltimes = np.append(final_totaltimes, np.sum(ScatteringTime))

# ------------------------------------------ END OF LOOP -------------------------------------------


# ------------------------------------------ OUTPUT -------------------------------------------


"""
Plot ouputs.
The first (initial) datapoint for r doesn't have a value, since only radius/volume CHANGES are used.
"""


def color_left_right(color1, color2, lw):
    ax_left.spines['left'].set_color(color1)
    ax_left.spines['left'].set_lw(lw)
    ax_right.spines['right'].set_color(color2)
    ax_right.spines['right'].set_lw(lw)
    # ax_left.spines['left'].set_bounds(10, 8620)
    # ax_right.spines['right'].set_bounds(0.001, 0.505)
    ax_right.spines['left'].set_visible(False)
    ax_left.spines['right'].set_visible(False)
    # ax_left.tick_params(axis='y', colors=uc_colors["blue"])


labelsize = 16
fontsize = labelsize
labelpad = 12

if iteration_dependency:
    fig, ax_left = plt.subplots()
    ax_right = ax_left.twinx()

    ticks = ticker.FuncFormatter(lambda x, pos: '{:.0f}e3'.format(x * 1e-3))
    # ax_left.yaxis.set_major_formatter(ticks)

    ax_left.set_ylabel(r"Final Scattering Rate $\Gamma_{th, fin}$ [kHz]", fontsize=fontsize, labelpad=labelpad)
    # ax_left.set_ylim(0, 1.1*Number_of_atoms[0])
    # ax_left.set_yscale('log')
    l1, = ax_left.plot(steps_list, final_scatteringrates * 1e-3,
                       label="Final Scattering Rate", color=uc_colors["blue"], lw=2)

    ax_right.set_ylabel(r"Total Time $t_{tot}$ [s]", fontsize=fontsize, labelpad=labelpad)
    ax_right.set_yscale("log")
    # ax_right.set_ylim(-0.005, 1.1*n[-1]*1e-18)
    l2, = ax_right.plot(steps_list, final_totaltimes,
                        label="Total times", color=uc_colors["red"], lw=2, ls="--")

    ax_left.spines['left'].set_color(uc_colors["blue"])
    ax_left.spines['left'].set_lw(2)
    ax_right.spines['right'].set_color(uc_colors["red"])
    ax_right.spines['right'].set_lw(2)
    # ax_right.spines['right'].set_linestyle("--")

    # ax_left.spines['left'].set_bounds(10, 8620)
    # ax_right.spines['right'].set_bounds(0.001, 0.505)

    ax_right.spines['left'].set_visible(False)
    ax_left.spines['right'].set_visible(False)
    # ax_left.tick_params(axis='y', colors=uc_colors["blue"])

    ax_left.set_xlabel("Iterations", fontsize=fontsize, labelpad=labelpad)
    # ax_left.set_xlim(r_initial*1.1e6, 0)

    ax_left.set_xscale('log')
    # ax_left.set_ylim([10e-1, 10000])
    ax_left.tick_params(axis='both', which='major', labelsize=labelsize)
    ax_right.tick_params(axis='both', which='major', labelsize=labelsize)
    ax_left.grid(alpha=0.5)
    ax_right.grid(linestyle="--", alpha=0.5)
    # plt.title("Atoms in a compressed trap: Number of atoms & Particle density \n")
    plt.legend([l1, l2], [r"Final Scattering Rate $\Gamma_{th, fin}$ ", r"Total time $t_{tot}$"],
               loc="lower right", ncol=1, fontsize=fontsize * .8)
    fig.tight_layout()
    # plt.savefig("../Blue_detuned_ODT_estimation/IterationsScRateTotalTime.pdf", format="pdf")
    plt.show()
    plt.clf()

    # %%
    fig, ax_left = plt.subplots()
    ax_right = ax_left.twinx()

    ax_left.set_ylabel(r"Final Density $\rho_{fin}$", fontsize=fontsize, labelpad=labelpad)
    # ax_left.set_ylim(0, 1.1*Number_of_atoms[0])
    # ax_left.set_yscale('log')
    l1, = ax_left.plot(steps_list, final_densities,
                       label="Final Particle Density", color=uc_colors["blue"], lw=2)

    ax_right.set_ylabel(r"Final Temperature $T_{fin}$ [mK]", fontsize=fontsize, labelpad=labelpad)
    # ax_right.set_ylim(0.055, 1.2*max(final_temperatures)*1e3)
    l2, = ax_right.plot(steps_list, final_temperatures * 1e3,
                        label="Final Temperature", color=uc_colors["red"], lw=2, ls="--")

    ax_left.spines['left'].set_color(uc_colors["blue"])
    ax_left.spines['left'].set_lw(2)
    ax_right.spines['right'].set_color(uc_colors["red"])
    ax_right.spines['right'].set_lw(2)
    # ax_right.spines['right'].set_linestyle("--")

    # ax_left.spines['left'].set_bounds(10, 8620)
    # ax_right.spines['right'].set_bounds(0.001, 0.505)

    ax_right.spines['left'].set_visible(False)
    ax_left.spines['right'].set_visible(False)
    # ax_left.tick_params(axis='y', colors=uc_colors["blue"])

    ax_left.set_xlabel("Iterations", fontsize=fontsize, labelpad=labelpad)
    # ax_left.set_xlim(r_initial*1.1e6, 0)

    ax_left.set_xscale('log')
    # ax_left.set_ylim([10e-1, 10000])
    ax_left.tick_params(axis='both', which='major', labelsize=labelsize)
    ax_right.tick_params(axis='both', which='major', labelsize=labelsize)
    ax_left.grid(alpha=0.5)
    ax_right.grid(linestyle="--", alpha=0.5)
    # plt.title("Atoms in a compressed trap: Number of atoms & Particle density \n")
    plt.legend([l1, l2], [r"Final Density $\rho_{fin}$", r"Final Temp. $T_{fin}$ [mK]"],
               loc="right", ncol=1, fontsize=fontsize * .8)
    fig.tight_layout()
    # plt.savefig("../Blue_detuned_ODT_estimation/IterationsAtomsTemp.pdf", format="pdf")
    plt.show()
    plt.clf()

# %%
if not iteration_dependency:
    # ------------------------------- Plot of Density & Atom Number --------------------------------
    fig, ax_left = plt.subplots()
    ax_right = ax_left.twinx()

    ticks = ticker.FuncFormatter(lambda x, pos: '{:.0f}e3'.format(x * 1e-3))
    ax_left.yaxis.set_major_formatter(ticks)

    ax_left.set_ylabel("Number of Atoms N", fontsize=fontsize, labelpad=labelpad)
    # ax_left.set_ylim(0, 1.1*Number_of_atoms[0])
    # ax_left.set_yscale('log')
    l1, = ax_left.plot(r * 1e6, Number_of_atoms, color=uc_colors["blue"], label="Number of atoms", lw=2)

    ax_right.set_ylabel(r"Particle Density $\rho$ [µm$^{-3}$]", fontsize=fontsize, labelpad=labelpad)
    ax_right.set_ylim(-0.005, 1.1 * n[-1] * 1e-18)
    l2, = ax_right.plot(r * 1e6, n * 1e-18, color=uc_colors["red"], label="density", lw=2, ls="--")

    ax_left.spines['left'].set_color(uc_colors["blue"])
    ax_left.spines['left'].set_lw(2)
    ax_right.spines['right'].set_color(uc_colors["red"])
    ax_right.spines['right'].set_lw(2)
    ax_right.spines['right'].set_linestyle("--")

    # ax_left.spines['left'].set_bounds(10, 8620)
    # ax_right.spines['right'].set_bounds(0.001, 0.505)

    ax_right.spines['left'].set_visible(False)
    ax_left.spines['right'].set_visible(False)
    # ax_left.tick_params(axis='y', colors=uc_colors["blue"])

    ax_left.set_xlabel("Inner Ring Radius $r$ [µm]", fontsize=fontsize, labelpad=labelpad)
    ax_left.set_xlim(r_initial * 1.1e6, 0)

    # ax_left.set_yscale('log')
    # ax_left.set_ylim([10e-1, 10000])
    ax_left.tick_params(axis='both', which='major', labelsize=labelsize)
    ax_right.tick_params(axis='both', which='major', labelsize=labelsize)
    ax_left.grid(alpha=0.5)
    ax_right.grid(linestyle="--", alpha=0.5)
    # plt.title("Atoms in a compressed trap: Number of atoms & Particle density \n")
    plt.legend([l1, l2], ["Number of atoms N", r"Particle density $\rho$"], loc="upper center", ncol=1,
               fontsize=fontsize * .8)
    fig.tight_layout()
    # plt.savefig("../Blue_detuned_ODT_estimation/DensityAtomNumber.pdf", format="pdf")
    plt.show()
    plt.clf()

    # %%
    # ------------------------------- Plot of Trap depth & Temperature --------------------------------

    fig, ax_left = plt.subplots()
    ax_right = ax_left.twinx()
    ylim1, ylim2 = 0, max(T[-1] * 1e3, U_inK[-1] * 1e3) * 1.35

    ax_left.set_ylabel("Gas Temperature $T$ [mK]", fontsize=fontsize, labelpad=labelpad)
    ax_left.set_ylim(ylim1, ylim2)
    l1, = ax_left.plot(r * 1e6, T * 1e3, color=uc_colors["blue"], label="Temperature", lw=2)

    ax_right.set_ylabel("Trap Depth $U(I)$ [mK]", fontsize=fontsize, labelpad=labelpad)
    ax_right.set_ylim(ylim1, ylim2)
    l2, = ax_right.plot(r * 1e6, U_inK * 1e3, color=uc_colors["red"], label="Trap depth", lw=2)

    color_left_right(uc_colors["blue"], uc_colors["red"], 2)

    ax_left.set_xlabel("Inner Ring Radius $r$ [µm]", fontsize=fontsize, labelpad=labelpad)
    ax_left.set_xlim((r_initial * 1.1e6), 0)
    ax_left.grid(alpha=0.5)
    # plt.title("Atoms in a compressed trap: Temperature & Trap depth \n")
    plt.legend([l1, l2], ["Gas Temperature $T$", "Trap Depth $U(I)$"], loc="upper center", ncol=1,
               fontsize=fontsize * .8)
    ax_left.tick_params(axis='both', which='major', labelsize=labelsize)
    ax_right.tick_params(axis='both', which='major', labelsize=labelsize)
    fig.tight_layout()
    # plt.savefig("../Blue_detuned_ODT_estimation/TemperatureTrapdepth.pdf", format="pdf")
    plt.show()
    plt.clf()

    # %%
    # ------------------------------- Plot of Scattering rate --------------------------------
    labelsize = 16
    fontsize = labelsize
    labelpad = 12
    fig, ax_left = plt.subplots()
    ax_right = ax_left.twinx()

    l1, = ax_right.plot(r * 1e6, 1e-3 * ScatteringRate, linewidth=2, color=uc_colors["red"], label="Scattering Rate")
    ax_right.set_xlabel("Inner Ring Radius $r$ [µm]", fontsize=fontsize, labelpad=labelpad)
    ax_right.set_xlim((r_initial * 1.1e6), 0)
    ax_right.tick_params(labelsize=labelsize)
    ax_right.set_ylabel(r"Scattering Rate $\Gamma_{th}$  [kHz]", fontsize=fontsize, labelpad=labelpad)
    ax_right.set_ylim(0, 1.05 * 1e-3 * ScatteringRate[-1])

    l2, = ax_left.plot(r * 1e6, 1e3 * ScatteringTime, linewidth=2, color=uc_colors["blue"], ls="--", label="test")
    ax_left.tick_params(labelsize=labelsize)
    ax_left.set_ylabel(r"Thermalization Time $t$ [ms]", fontsize=fontsize, labelpad=labelpad)
    ax_left.set_ylim(0, 1e3 * 1.1 * max(ScatteringTime))
    ax_left.set_xlabel("Inner Ring Radius $r$ [µm]", fontsize=fontsize, labelpad=labelpad)

    color_left_right(uc_colors["blue"], uc_colors["red"], 2)
    ax_right.grid(alpha=0.5)
    ax_left.grid(linestyle="--", alpha=0.5)
    # plt.title("Compressed trap: An adiabatic simulation \n Scattering rate", fontsize="x-large")

    plt.legend([l1, l2], ["Scattering Rate $\Gamma_{th}$ ", "Thermalization Time $t$ "], loc="upper center",
               ncol=1, fontsize=fontsize * .8)
    plt.tight_layout()
    # plt.savefig("../Blue_detuned_ODT_estimation/Scatteringrate.pdf", format="pdf")
    plt.show()
    plt.clf()
