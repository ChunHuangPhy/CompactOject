import csv
import numpy as np
import sys

def EOS_import(file_name = "", density = 0, pressure = 0):

    """EOS_import

    Imports density and pressure from csv or array, checks them, and returns them.

    Args:
        file_name (string, optional): string. CSV file to be opened.
        density (array, optional): numpy 1Darray. Passed into a check function and returned if valid.
        pressure (array, optional): numpy 1Darray. Passed into a check function and returned if valid.

    Returns:
        array: checked density and pressure.
    """
    
    if not file_name:
        density_checked, pressure_checked = EOS_check(density, pressure)
        return density_checked, pressure_checked

    input_file = file_name

    density, pressure = file_read(input_file)

    density_checked, pressure_checked = EOS_check(density, pressure)

    return density, pressure

def file_read(input_file):
    """file_read

    Reads a csv file of denisty and pressure given by the user.

    Args:
        input_file (string): string. File to be opened and parsed.

    Returns:
        array: two 1Darray numpy arrays, one corresponding to density and one corresponding to pressrure. 
    """
    
    data_list = []
    density_list = []
    pressure_list = []
    with open(input_file) as csvfile:
        file_read = csv.reader(csvfile, delimiter=',')
        data_list = [row for row in file_read]
    for row in data_list:
        density_list.append(float(row[0]))
        pressure_list.append(float(row[1]))

    # Make the lists numpy arrays 
    density_array = np.array(density_list)
    pressure_array = np.array(pressure_list)

    return density_array, pressure_array

def EOS_check(density, pressure):

    """file_read

    Checks that the derivative (drho/dp) is positive.

    Args:
        density (array): numpy 1Darray. Density array to be checked.
        pressure (array): numpy 1Darray. Pressure array to be checked.

    Returns:
        array: two arrays, one corresponding to density and one corresponding to pressrure or ends the function and prints
        invalid equation of state.
    """

    dp = np.diff(pressure) # dy
    drho = np.diff(density) # dx

    for value in drho:
        if value == 0:
            print("This is not a valid equation of state")
            sys.exit()

    dpdrho = dp/drho # dydx

    for value in dpdrho:
        if value >= 0:
            print("This is a valid equation of state")
            pass
        else:
            print("This is not a valid equation of state")
            sys.exit()
            
        

    return density, pressure
