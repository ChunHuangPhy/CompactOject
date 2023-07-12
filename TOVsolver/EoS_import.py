import csv
import numpy as np
import sys

def EOS_import(file_name = "", density = 0, pressure = 0):
    
    if not file_name:
        density_checked, pressure_checked = EOS_check(density, pressure)
        return density_checked, pressure_checked

    input_file = file_name

    density, pressure = file_read(input_file)

    density_checked, pressure_checked = EOS_check(density, pressure)

    return density, pressure

def file_read(input_file):
    data_list = []
    density_list = []
    pressure_list = []
    with open(input_file) as csvfile:
        file_read = csv.reader(csvfile, delimiter=' ')
        data_list = [row for row in file_read]
    for row in data_list:
        density_list.append(float(row[0]))
        pressure_list.append(float(row[1]))

    # Make the lists numpy arrays 
    density_array = np.array(density_list)
    pressure_array = np.array(pressure_list)

    return density_array, pressure_array

def EOS_check(density, pressure):

    dydx = np.gradient(density,pressure)
    for value in dydx:
        if value >= 0:
            pass
        else:
            print("This is not a valid equation of state")
            sys.exit()
    print("This is a valid equation of state")

    return density, pressure
