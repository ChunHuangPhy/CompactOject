# Python packages
import numpy as np


# Import files
import TOVsolver.solver_code as TOV_solver
import TOVsolver.EoS_import as EoS_import
import TOVsolver.speed_of_sound as speed_of_sound
from TOVsolver.unit import g_cm_3

# Global Variables
def OutputMR(input_file="", density=[], pressure=[]):
    """Outputs the mass, radius, and tidal deformability
    Args:
        file_name (string, optional): string. CSV file to be opened.
        density (array, optional): numpy 1Darray. Passed into a check function and returned if valid.
        pressure (array, optional): numpy 1Darray. Passed into a check function and returned if valid.

    Returns:
        MR (tuple): tuple with mass, radius.
    """

    #############This is something we need to change, like the input for this EOS import should
    ############# be one file contatining Whole EOS. that first column is density and second is pressure
    energy_density, pressure = EoS_import.EOS_import(input_file, density, pressure)
    ############# Lets the user only input the EOS file path, then this EOS_import should have file
    ############# as input. and the outputMR should have a file as input too?

    Radius = []
    Mass = []

    density = np.logspace(14.3, 15.6, 50) * g_cm_3
    # This following step is to make a dicision whether the EOS ingredients is always increase. We can do that outsie of this main to the
    # EOS import.
    # if   all(x<y for x, y in zip(eps_total_poly[:], eps_total_poly[[1:])) and all(x<y for x, y in zip(pres_total_poly[j][:], pres_total_poly[j][1:])):
    for i in range(len(density)):
        try:
            M, R = TOV_solver.solveTOV(density[i], energy_density, pressure)
            Mass.append(M)
            Radius.append(R)
        # This is sentense is for avoiding the outflow of the result, like when solveTOV blow up because of ill EOS, we need to stop
        except OverflowError as e:
            # print("This EOS is ill-defined to reach an infinity result, that is not phyiscal, No Mass radius will be generated.")
            break
    MR = np.vstack((Radius, Mass)).T
    # print("Mass Radius file will be generated and stored as  2-d array. The first column is Radius, second one is mass")

    return MR


def OutputMRT(input_file="", density=[], pressure=[]):
    """Outputs the mass, radius, and tidal deformability
    Args:
        file_name (string, optional): string. CSV file to be opened.
        density (array, optional): numpy 1Darray. Passed into a check function and returned if valid.
        pressure (array, optional): numpy 1Darray. Passed into a check function and returned if valid.

    Returns:
        MRT (tuple): tuple with mass, radius, and tidal deformability.
    """

    #############This is something we need to change, like the input for this EOS import should
    ############# be one file contatining Whole EOS. that first column is density and second is pressure
    energy_density, pressure = EoS_import.EOS_import(input_file, density, pressure)
    ############# Lets the user only input the EOS file path, then this EOS_import should have file
    ############# as input. and the outputMR should have a file as input too?

    Radius = []
    Mass = []
    tidal = []
    density = np.logspace(14.3, 15.6, 50) * g_cm_3
    # This following step is to make a dicision whether the EOS ingredients is always increase. We can do that outsie of this main to the
    # EOS import.
    # if   all(x<y for x, y in zip(eps_total_poly[:], eps_total_poly[[1:])) and all(x<y for x, y in zip(pres_total_poly[j][:], pres_total_poly[j][1:])):
    for i in range(len(density)):
        try:
            M, R, T = TOV_solver.solveTOV_tidal(density[i], energy_density, pressure)
            Radius.append(R)
            Mass.append(M)
            tidal.append(T)
        # This is sentense is for avoiding the outflow of the result, like when solveTOV blow up because of ill EOS, we need to stop
        except OverflowError as e:
            # print("This EOS is ill-defined to reach an infinity result, that is not phyiscal, No Mass radius will be generated.")
            break
    MRT = np.vstack((Radius, Mass, tidal)).T
    # print("Mass Radius and tidal will be generated as the 3-d array. The first column is Radius, second one is mass,last is tidal")

    return MRT


def OutputC_s(input_file="", density=[], pressure=[]):
    """Calls function to open csv (if needed) and check equation of state validity.
        Then calls function to calculate speed of sound.

    Args:
        file_name (string, optional): string. CSV file to be opened.
        density (array, optional): numpy 1Darray. Passed into a check function and returned if valid.
        pressure (array, optional): numpy 1Darray. Passed into a check function and returned if valid.

    Returns:
        C_s (array): numpy 1D array. List of speeds of sound.
    """

    energy_density, pressure = EoS_import.EOS_import(input_file, density, pressure)
    C_s = speed_of_sound.speed_of_sound_calc(energy_density, pressure)
    return C_s


def OutputMRpoint(central_density, energy_density, pressure):
    """Outputs the mass, radius, and tidal deformability (single point)
    Args:
        central_density (float): central density that we want to compute
        density (array, optional): numpy 1Darray. Density of EoS
        pressure (array, optional): numpy 1Darray. pressure of EoS

    Returns:
        MR (tuple): tuple with mass, radius.
    """
    Radius = []
    Mass = []

    # This following step is to make a dicision whether the EOS ingredients is always increase. We can do that outsie of this main to the
    # EOS import.
    # if   all(x<y for x, y in zip(eps_total_poly[:], eps_total_poly[[1:])) and all(x<y for x, y in zip(pres_total_poly[j][:], pres_total_poly[j][1:])):
    try:
        M, R = TOV_solver.solveTOV(central_density* g_cm_3, energy_density, pressure)
        Mass.append(M)
        Radius.append(R)
    # This is sentense is for avoiding the outflow of the result, like when solveTOV blow up because of ill EOS, we need to stop
    except OverflowError as e:
        print(
            "This EOS is ill-defined to reach an infinity result, that is not phyiscal, No Mass radius will be generated."
        )
    MR = np.vstack((Radius, Mass)).T

    return MR


def OutputMRTpoint(central_density, energy_density, pressure):
    """Outputs the mass, radius, and tidal deformability (single point)
    Args:
        central_density (float): central density that we want to compute
        density (array, optional): numpy 1Darray. Density of EoS
        pressure (array, optional): numpy 1Darray. pressure of EoS

    Returns:
        MRT (tuple): tuple with mass, radius and tidal.
    """
    Radius = []
    Mass = []
    tidal = []
    # This following step is to make a dicision whether the EOS ingredients is always increase. We can do that outsie of this main to the
    # EOS import.
    # if   all(x<y for x, y in zip(eps_total_poly[:], eps_total_poly[[1:])) and all(x<y for x, y in zip(pres_total_poly[j][:], pres_total_poly[j][1:])):
    try:
        M, R, T = TOV_solver.solveTOV_tidal(central_density* g_cm_3, energy_density, pressure)
        Mass.append(M)
        Radius.append(R)
        tidal.append(T)
    # This is sentense is for avoiding the outflow of the result, like when solveTOV blow up because of ill EOS, we need to stop
    except OverflowError as e:
        print(
            "This EOS is ill-defined to reach an infinity result, that is not phyiscal, No Mass radius will be generated."
        )
    MRT = np.vstack((Radius, Mass, tidal)).T

    return MRT


def OutputMRpoint_with_EoS(central_density, Pmin, eos, inveos):
    """Outputs the mass, radius(single point)
    Args:
        central_density (float): central density that we want to compute, in unit.g_cm_3
        Pmin (float):  In unit.G / unit.c**4
        eos (function): Equation of state, rho in unit.G / unit.c**2, pressure in unit.G / unit.c**4
        inveos (function): The inverse of equation of state

    Returns:
        Mass (float): mass
        Radius (float): radius
    """
    try:
        Mass, Radius = TOV_solver.solveTOV2(central_density, Pmin, eos, inveos)
        return Mass, Radius
    # This is sentense is for avoiding the outflow of the result, like when solveTOV blow up because of ill EOS, we need to stop
    except OverflowError as e:
        print(
            "This EOS is ill-defined to reach an infinity result, that is not phyiscal, No Mass radius will be generated."
        )
        return (0, 0)
