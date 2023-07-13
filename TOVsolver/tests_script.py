import main
import EoS_import
import constant
import speed_of_sound
import solver_code
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def test_constant():
    assert constant.c == 3e10
    assert constant.Msun == 1.989e33
    assert constant.dyncm2_to_MeVfm3 == 1./(1.6022e33)
    assert constant.gcm3_to_MeVfm3 == 1./(1.7827e12)
    assert constant.oneoverfm_MeV == 197.33
    assert type(constant.c) == float 
    assert type(constant.Msun) == float
    assert type(constant.dyncm2_to_MeVfm3) == float
    assert type(constant.gcm3_to_MeVfm3) == float
    assert type(constant.oneoverfm_MeV) == float

def test_EOS_import_two_files():
    pres = pd.read_csv('EoS_inference/Test Case/new_pres.csv')
    dens = pd.read_csv('EoS_inference/Test Case/new_eps.csv')
    for indices in range(len(pres.columns)):
        a, b = EoS_import.EOS_import(density = list(dens[str(indices)]), pressure = list(pres[str(indices)]))
        assert type(a) == list
        assert type(b) == list

def test_EOS_import_one_file():
    a, b = EoS_import.EOS_import('/home/nwhitsett/codeastro/EoS_inference/Test Case/Test_EOS.csv')
    assert type(a) == np.ndarray
    assert type(b) == np.ndarray
    
def test_file_read():
    a, b = EoS_import.file_read('EoS_inference/Test Case/Test_EOS.csv')
    assert type(a) == np.ndarray
    assert type(b) == np.ndarray

def test_EOS_check():
    a, b = EoS_import.file_read('EoS_inference/Test Case/Test_EOS.csv')
    c, d = EoS_import.EOS_check(a, b)
    assert type(c) == np.ndarray
    assert type(d) == np.ndarray

def test_end_to_end():
    a = main.OutputMR('EoS_inference/Test Case/Test_EOS.csv')
    b = main.OutputC_s('EoS_inference/Test Case/Test_EOS.csv')
    print(type(a), type(b))
test_end_to_end()
