---
title: 'CompactObject: An open-source Python package for full-scope neutron star equation of state inference'
tags:
    - Python
    - astrostatistics
    - neutron stars
authors:
  - name: Chun Huang
    orcid: 0000-0001-6406-1003
    affiliation: 1
  - name: Tuhin Malik
    orcid: 0000-0003-2633-5821
    affiliation: 2
  - name: João Cartaxo
    orcid: 0009-0001-7105-8272
    affiliation: 2
  - name: Shashwat Sourav
    orcid: 0000-0002-0169-4003
    affiliation: 1
  - name: Wenli Yuan
    orcid: 0000-0003-2771-759X
    affiliation: 3
  - name: Tianzhe Zhou
    orcid: 0009-0000-8504-9134
    affiliation: 4
  - name: Xuezhi Liu
    orcid: 
    affiliation: 5
# - name: Xieyuan Dong
#   orcid: 
#   affiliation: 6
  - name: John Groger
    orcid: 0000-0002-7054-9053
    affiliation: 1
  - name: Nicole Osborn
    orcid: 
    affiliation: 1
  - name: Nathan Whitsett
    orcid: 
    affiliation: 1
  - name: Zhiheng Wang
    orcid: 0009-0000-5088-6207
    affiliation: 1
  - name: Constan\c{c}a Provid\^{e}ncia
    orcid: 0000-0001-6464-8023
    affiliation: 2
  - name: Micaela Oertel
    orcid: 0000-0002-1884-8654
    affiliation: 7
  - name: Alexander Y. Chen
    orcid: 0000-0002-4738-1168
    affiliation: 1
  - name: Laura Tolos
    orcid: 0000-0002-6449-106X
    affiliation: 8
    affiliation: 9
    affiliation: 10
  - name: Anna Watts
    orcid: 0000-0002-1009-2354
    affiliation: 11

affiliations:
  - name: Physics Department and McDonnell Center for the Space Sciences, Washington University in St. Louis, MO 63130, USA
    index: 1
  - name: CFisUC, Department of Physics, University of Coimbra, 3004-516 Coimbra, Portugal
    index: 2
  - name: School of Physics and State Key Laboratory of Nuclear Physics and Technology, Peking University, Beijing 100871, China
    index: 3
  - name: Department of Physics, Tsinghua University, Beijing 100084, China
    index: 4
  - name: Physics Department, Central China Normal University, Luoyu Road, 430030, Wuhan, China
    index: 5
#  - name: School of Physics, Nankai University, Tianjin 300071, China
#    index: 6
  - name: Laboratoire Univers et Th\^{e}ories, CNRS, Observatoire de Paris, Universit\^{e} PSL, Universit\^{e} Paris Cit\^{e}, 5 place Jules Janssen, 92195 Meudon, France
    index: 7
  - name: Institute of Space Sciences (ICE, CSIC), Campus UAB, Carrer de Can Magrans, 08193, Barcelona, Spain
    index: 8
  - name: Institut d'Estudis Espacials de Catalunya (IEEC), 08860 Castelldefels (Barcelona), Spain
    index: 9
  - name: Frankfurt Institute for Advanced Studies, Ruth-Moufang-Str. 1, 60438 Frankfurt am Main, Germany
    index: 10
  - name: Anton Pannekoek Institute for Astronomy, University of Amsterdam, Science Park 904, 1090 GE Amsterdam, the Netherlands
    index: 11

date: Oct 23 2021
bibliography: cojoss.bib
---


# Summary

The CompactObject package is an open-source software framework developed to constrain the neutron star equation of state (EOS) through Bayesian statistical inference. It integrates astrophysical observational constraints from X-ray timing, gravitational wave events, and radio measurements, as well as nuclear experimental constraints derived from perturbative Quantum Chromodynamics (pQCD) and Chiral Effective Field Theory ($\chi$EFT). The package supports a diverse range of EOS models, including meta-models like and several physics-motivated EOS models. It comprises three independent components: an EOS generator module that currently provided seven EOS choices, a Tolman–Oppenheimer–Volkoff (TOV) equation solver, enabling solve Mass Radius and Tidal deformability as observables, and a comprehensive Bayesian inference workflow module, include whole pipline of implementing EOS Bayesian inference. Each component can be independently utilized in various scientific research contexts, like nuclear physics and astrophysics. Additionally, CompactObject is designed to synergize with existing software such as CompOSE, enabling the use of the CompOSE EOS database to expand the available EOS options.

# Statement of need

Understanding the equation of state (EOS) of neutron stars is important for understanding the fundamental physics governing ultra-dense matter. Neutron stars, with it core densities exceeding several nuclear saturation density, have a crucial role to play in studying nuclear interactions under extreme conditions. However, inferring the EOS from observational and experimental data has significant challenges due to the complex interplay of astrophysical phenomena and nuclear physics. Existing tools, such as [@Raaijmakers2023], primarily focus on meta-model constraints and lack integration with physics-motivated EOS models. In contrast, CompactObject achieves high accuracy and rapid computation for this family of EOS, thereby enabling researchers to perform inferences based on physically motivated models and apply nuclear physics-related constraints derived from nuclear experiments. 

CompactObject appears as a viable solution to these challenges by providing an open-source, robust platform designed for Bayesian inference on neutron star EOS constraints. Its comprehensive workflow integrates a wide range of EOS, including not only physical and meta-models of neutron star EOS but also strange star and quark star EOS, which have been proposed to explain the nature of these compact objects, enabling a detailed exploration of dense matter physics. The package's user-friendly interface and modular architecture facilitate easy adoption and extension, allowing researchers to customize analyses and incorporate new EOS as they become available. Furthermore, thorough documentation ensures that both novice and experienced users can effectively utilize the tool, promoting widespread accessibility and collaborative advancement in the field. By addressing the need for an integrated, flexible, and well-documented framework, CompactObject enhances the capability of nuclear astrophysicists to derive precise EOS constraints. 


# The Compactobject package and science use

CompactObject is an open-source software package designed to apply astrophysical and nuclear physics constraints to equation of state (EOS) parameters. Currently, the available EOS options include polytropic EOS and speed of sound model EOS, both of which are meta-models. Additionally, the package supports physics-motivated models such as the Relativistic Mean Field (RMF) theory and its density-dependent variant. Beyond neutron star EOS models, CompactObject also includes a strange star EOS based on the strangeon model and the widely used MIT bag model for quark stars.

The package integrates various likelihood constraints, including routines for simulating mass-radius measurements from X-ray timing observations and analyzing mass-radius likelihoods from actual observational data. It also incorporates constraints from radio timing observations, gravitational wave observations related to tidal deformability, and nuclear physics constraints derived from saturation properties, perturbative Quantum Chromodynamics (pQCD), and Chiral Effective Field Theory ($\chi$EFT).

Furthermore, CompactObject includes routines for EOS analysis that output additional properties of neutron stars, such as proton fraction and the number densities of different particles within the star. Other than these, CompactObject synergizes with the CompOSE database, allowing users to derive observational evidence directly into existing EOS models. The nested sampling pipeline implemented in CompactObject is based on `UltraNest`, providing a computational framework to extract Bayesian evidence for each integrated EOS model. For the inference pipeline, the package offers two sampling algorithm options: `UltraNest` (nested sampling) and `emcee` (Markov Chain Monte Carlo sampling).

CompactObject has been utilized to derive constraints on nucleonic RMF models [@Huang:2023grj] and hyperonic RMF models [@Huang:2024rvj]. Ongoing projects include constraining the strangeon star EOS [@Yuan2024] and exploring phase transitions and twin stars [@Huang24; @Zhou24]. Additionally, nuclear physics constraints and density-dependent RMF constraints are being developed *(NEED CITATIONS)*.

The released version of CompactObject is readily accessible through its GitHub repository [@EoS_inference] under the MIT license and is archived on the Zenodo repository [@COZenodo]. Comprehensive documentation and the complete workflow for implementing EOS inference are available in the GitHub repository. Future plans for CompactObject include expanding the range of available EOS options and conducting a detailed survey of existing EOS models to perform cross-comparisons of Bayesian evidences using current observational and experimental constraints.

Software: Python language [@10.1109/MCSE.2007.58], NumPy [@van_der_Walt_2011], MPI for Python [@DALCIN2008655], Numba [@numba], NumbaMinpack [@wogan_NumbaMinpack], Matplotlib [@Hunter:2007], Jupyter [@2016ppap.book...87K], UltraNest [@2021JOSS....6.3001B], emcee [@Foreman_Mackey_2013], SciPy [@2020SciPy-NMeth], Seaborn [@Waskom2021], corner.py [@corner]

## Acknowledgements

H.C. acknowledges support from the Art & Science Fellowship of Washington University in St Louis.

# References