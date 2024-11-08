"""
MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version: 08/11/2024
Adapated from parts of the CUTER tool, see https://zenodo.org/records/10781539 

Module that reads tabulated EoS in Compose format or in LALSuite format
"""

import numpy as np
#from TOVsolver.unit import g_cm_3, dyn_cm_2
from pathlib import Path
import sys
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join

c_si = 2.99792458E+8 	 #/< Velocity of light [m/s]
mev_si =1.602176634E-13   #/< One MeV [J]
MeVfm3_to_gcm3 = mev_si/(c_si*c_si)*1e42 # convert energy density from MeV/fm^3 to g/cm^3
MeVfm3_to_dyncm2 = mev_si*1e46  # convert pressure from MeV/fm^3 to dyne/cm^2
  # to convert from old Lal units to nuclear units (MeV/fm^3)
lalpconvert = 755397844024.145 # pressure
laleconvert = 755407089027.154 # energy density

def read_compose(eosdir="./filesCompose", eos_prefix="/eos",nptsmin = 100,eosname=None): 

    """
    Routine which reads the data from the original EoS model, Compose format is assumed

    One first check is performed on the data:
    1. Minimum number of grid points > nptsmin (default is 100)

    Args:
      eosdir (character string): name of the directory where the EoS files are stored (default is filesCompose)
      nptsmin (int): required minimum number of grid points
      eos_prefix (character string): prefix of the eos file names, default is "eos"


    Returns:
       energy density (g/cm^3) : numpy array 
       pressure (dyne/cm^2): numpy array
       eosname: character string with name of the EoS stored in the table

    """

    # initialising the return variables
    ntab = 1

    nB_compose = np.zeros(ntab)
    p_compose = np.zeros(ntab)
    eps_compose = np.zeros(ntab)
    temp = 0.

    # Check if there is a temperature file, and if it is empty
    if Path(eosdir+eos_prefix+".t").is_file():
        with open(eosdir+eos_prefix+".t") as file_temp:
          lines_t = file_temp.readlines()
        ndimt = int(lines_t[1]) - int(lines_t[0])
    else : ndimt = 0

    if ndimt == 0 : # No temperature file or empty temperature file
        with open(eosdir+eos_prefix+".nb") as file_nb:
            lines_nb = file_nb.readlines()
        with open(eosdir+eos_prefix+".thermo") as file_thermo:
            lines_th = file_thermo.readlines()
    else: # if temperature file not empty, check that the minimum T =0
          temp = float(lines_t[2]) # lowest temperature entry of the EoS model
          if(temp > 0.) :
              sys.exit("\t Nonzero temperature in EoS, aborting...")
          else:
              with open(eosdir+eos_prefix+".nb") as file_nb:
                  lines_nb = file_nb.readlines()
              with open(eosdir+eos_prefix+".thermo") as file_thermo:
                  lines_th = file_thermo.readlines()
              
    ndim = int(lines_nb[1]) - int(lines_nb[0]) + 1
    if(ndim <=nptsmin):
          print("\t EoS table with",ndim,"lines is not sufficient, need at least ",nptsmin)
          sys.exit()



    nB_compose = np.zeros(ndim)
    for i in range(ndim):
          nB_compose[i] = float(lines_nb[i+2])

    p_compose = np.zeros(ndim)
    eps_compose = np.zeros(ndim)


    # Reading the thermo file
    tmp = lines_th[0].split()
    m_n = float(tmp[0])
    m_p = float(tmp[1])
    print("\t Nucleon masses used within the CompOSE tables :\t m_n = %.5f MeV and m_p = %.5f MeV"%(m_n, m_p))

    for i in range(ndim):
          tmp = lines_th[i+1].split()
          p_compose[i]=float(tmp[3])*nB_compose[i]  # pressure
          eps_compose[i]=(float(tmp[9])+ 1.0)*m_n*nB_compose[i]  # energy density

    if eosname is None:
        eosname = read_README(eosdir)
    eps_gcm3 = eps_compose*MeVfm3_to_gcm3
    p_dyncm2 = p_compose*MeVfm3_to_dyncm2
    
    return eps_gcm3,p_dyncm2,eosname



def read_Lal(filename,eosdir = "./filesLaL/"):
    """
    Read the equation of state file in the LAL format, 
    a number of tables is available at 
    https://git.ligo.org/lscsoft/lalsuite/-/tree/master/lalsimulation/lib

    Args:
      eosdir (character string): name of the directory where the EoS files are stored (default is filesCompose)
      filename: name of the LAL format equation of state file with extension
    Returns:
       energy density (g/cm^3) : numpy array 
       pressure (dyne/cm^2): numpy array
       eosname: character string with name of the EoS stored in the table

    """
    print("\n \t \t Reading the input table...")
    with open(eosdir+filename) as f:
        line = f.readline()
    ncolumns = len(line.split())
    if ncolumns == 2: # old Lal format
        ori_p,ori_eps = np.loadtxt(eosdir+filename,unpack=True)
        ori_p = ori_p*lalpconvert*MeVfm3_to_dyncm2
        ori_eps = ori_eps*laleconvert*MeVfm3_to_gcm3
    else:    # new Lal format
        ori_nb, ori_eps, ori_p, ori_mub, ori_mue, ori_yemu = np.loadtxt(eosdir+filename,usecols=(1,2,3,4,5,7), unpack =True, skiprows=8)
    # nb[1/fm^3], eps[g/cm^3], P[dyn/cm^2], mu_B [MeV] mu_e [MeV]

    print("\t \t ... done.")

    # Define the list of strings that need to be replaced to get the name
    # of the EoS from the filename
    strings = ["LALSimNeutronStar",".dat"]
    eosname = filename
    for string in strings:
        eosname = eosname.replace(string,"")

    
    return ori_eps,ori_p,eosname


def read_README(eosdir):
  """
  Routine to read the information about the EoS model from the README(.txt) file

  Args:
  eosdir: directory in which the README file is

  Returns:
  eosname: name of the eos.

  """

  if Path(eosdir+"/README").is_file():
    file_eos_name = open(eosdir+"/README","r")
    eosname = file_eos_name.readline()
  else:
    if Path(eosdir+"/README.txt").is_file():  
      file_eos_name = open(eosdir+"/README.txt","r")
      eosname = file_eos_name.readline()
    else:
      print("Warning: could not determine EoS name, README file missing")
      print("README is in general contained in download of eos.zip from Compose")
      print("Eosname is then set to default value = core_eos")
      eosname = "core_eos"

      # Define the list of strings that need to be replaced
  strings = ['\n', "downloaded from ",'http://','https://','compose.obspm.fr','Original','with','electrons','references','cold NS',':','for','EoS','table','Table','model', 'Eos', 'cold','beta','equilibrated', 'equilibrium','neutron','star','in','star','unified', 'matter', 'NS',' ',')','.',',','crust']
  for string in strings:
    eosname = eosname.replace(string,"")
  eosname = eosname.replace("(","_")
  eosname = eosname.replace("-","_")
  return eosname


#### main for testing purposes only
def test_read_eos_tables():
    e_compose,p_compose,name_compose = read_compose()
    lalpath = "./filesLaL/"

    lalfiles = [f for f in listdir(lalpath) if isfile(join(lalpath, f))]
    print("Will read the following files in LaL format:")
    nlal = len(lalfiles)
    if nlal > 0:
        ii = 0
        print(lalfiles[ii])
        e_lal,p_lal,name_lal = read_Lal(filename=lalfiles[ii])
        epslal = [e_lal] 
        plal = [p_lal] 
        namelal = [name_lal] 
        for ii in range(1,nlal):
            print(lalfiles[ii])
            e_lal,p_lal,name_lal = read_Lal(filename=lalfiles[ii])
            epslal.append(e_lal) 
            plal.append(p_lal)
            namelal.append(name_lal) 
    else:
        print("No file in LAL format found")

    plt.yscale('log')
    plt.xscale('log')

    plt.plot(e_compose,p_compose,'-',color = 'blue',label=name_compose)
    for ii in range(nlal):
        plt.plot(epslal[ii],plal[ii],label=namelal[ii])

    plt.legend()
    plt.xlabel('energy density (g/cm$^3$)')
    plt.ylabel('pressure (dyne/cm$^2$)')
    savename = "testreadeos.pdf"
    plt.savefig(savename,bbox_inches='tight')
    plt.show()

    return
