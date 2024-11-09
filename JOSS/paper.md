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
  - name: John Groger
    orcid: 0000-0002-7054-9053
    affiliation: 1
  - name: Nicole Osborn
    orcid: 
    affiliation: 1
  - name: Nathan Whitsett
    orcid: 
    affiliation: 1
  - name: Xieyuan Dong
   orcid: 
   affiliation: 6
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
  - name: School of Physics, Nankai University, Tianjin 300071, China
    index: 6
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

The CompactObject package is an open-source software framework developed to constrain the neutron star equation of state (EOS) through Bayesian statistical inference. It integrates astrophysical observational constraints from X-ray timing, gravitational wave events, and radio measurements, as well as nuclear experimental constraints derived from perturbative Quantum Chromodynamics (pQCD) and Chiral Effective Field Theory ($\chi$EFT). The package supports a diverse range of EOS models, including meta-model like and several physics-motivated EOS models. It comprises three independent components: an EOS generator module that currently provided seven EOS choices, a Tolman–Oppenheimer–Volkoff (TOV) equation solver, enabling solve Mass Radius and Tidal deformability as observables, and a comprehensive Bayesian inference workflow module, including a whole pipeline of implementing EOS Bayesian inference. Each component can be independently utilized in various scientific research contexts, like nuclear physics and astrophysics. Additionally, CompactObject is designed to synergize with existing software such as [CompOSE](https://compose.obspm.fr), enabling the use of the CompOSE EOS database to expand the available EOS options.

# Statement of need

Understanding the equation of state (EOS) of neutron stars is important for understanding the fundamental physics governing ultra-dense matter. Neutron stars, with its core densities exceeding several time nuclear saturation density, have a crucial role to play in studying nuclear interactions under extreme conditions. However, inferring the EOS from observational and experimental data has significant challenges due to the complex interplay of astrophysical phenomena and nuclear physics. Many of these studies such as [@Raaijmakers2023] focus on EOS meta-models (which may be parameterized or non-parameterized) that attempt to span all reasonable mass-radius parameter space, rather than being driven by microphysics. In contrast, CompactObject achieves high accuracy and rapid computation for a family of physics-motivated EOSs, thereby enabling researchers to perform inferences based on physically motivated models and apply nuclear physics-related constraints derived from nuclear experiments.

CompactObject appears as a viable solution to these challenges by providing an open-source, robust platform designed for Bayesian inference on neutron star EOS constraints. Its comprehensive workflow integrates a wide range of EOSs, including not only physical and meta-models of neutron star EOSs but also strange star and quark star EOSs, which have been proposed to explain the nature of these compact objects, enabling a detailed exploration of dense matter physics. The package's user-friendly interface and modular architecture facilitate easy adoption and extension, allowing researchers to customize analyses and incorporate new EOSs as they become available. Furthermore, thorough documentation ensures that both novice and experienced users can effectively utilize the tool, promoting widespread accessibility and collaborative advancement in the field. By addressing the need for an integrated, flexible, and well-documented framework, CompactObject enhances the capability of nuclear astrophysicists to derive precise EOS constraints. 


# The Compactobject package and science use

CompactObject is an open-source software package designed to apply astrophysical and nuclear physics constraints to EOS parameters. Currently, the available EOS options include polytropic EOSs and speed of sound model EOSs, both of which are meta-models. Additionally, the package supports physics-motivated models such as the Relativistic Mean Field (RMF) theory [@Tolos_2016,@Tolos_Centelles_Ramos_2017] and its density-dependent variant [@Hempel_2010,@Char_2014]. We integrated features for users to define the density dependent variant form by themselves, and which span most of the possibility of this family of models. Beyond neutron star EOS models, CompactObject also includes a strange star EOS based on the strangeon model [@2003ApJ...596L..59X] and the widely used MIT bag model [@PhysRevD.9.3471] for quark stars. 

The package integrates various likelihood constraints, including routines for simulating mass-radius measurements from X-ray timing observations and analyzing mass-radius likelihoods from actual observational data. It also incorporates constraints from radio timing observations, gravitational wave observations related to tidal deformability, and nuclear physics constraints derived from saturation properties, pQCD [@Gorda:2022jvk] [@Providencia:2023rxc]
and $\chi$EFT [@Hebeler:2013nza] [@Huth:2021bsp]

Furthermore, CompactObject includes routines for EOS analysis that output additional properties of neutron stars, such as proton fraction and the number densities of different particles within the star. Other than these, CompactObject synergizes with the CompOSE database, allowing users to derive observational evidence directly into existing EOS models. The nested sampling pipeline implemented in CompactObject is based on UltraNest, providing a computational framework to extract Bayesian evidence for each integrated EOS model. For the inference pipeline, the package offers two sampling algorithm options: UltraNest (nested sampling) and emcee (Markov Chain Monte Carlo sampling). 

CompactObject has been utilized to derive constraints on nucleonic RMF models [@Huang:2023grj] and hyperonic RMF models \citep[@Huang:2024rvj]. Ongoing projects include constraining the strangeon star EOS and exploring phase transitions and twin stars. Additionally, various nuclear and astrophysical constraints on RMF model with density-dependent couplings are being developed [@Malik:2022zol]

The released version of CompactObject is readily accessible through its GitHub repository [@EoS_inference] under the MIT license and is archived on Zenodo repository [@COZenodo]. Comprehensive documentation and the complete workflow for implementing EOS inference are available in the GitHub repository. Future plans for CompactObject include expanding the range of available EOS options and conducting a detailed survey of existing EOS models to perform cross-comparisons of Bayesian evidence using current observational and experimental constraints.

Software: Python language [@10.1109/MCSE.2007.58], NumPy [@van_der_Walt_2011], MPI for Python [@DALCIN2008655], Numba [@numba], NumbaMinpack [@wogan_NumbaMinpack], Matplotlib [@Hunter:2007], Jupyter [@2016ppap.book...87K], UltraNest [@2021JOSS....6.3001B], emcee [@Foreman_Mackey_2013], SciPy [@2020SciPy-NMeth], Seaborn [@Waskom2021], corner.py [@corner]

## Acknowledgements

H.C., S.S., O.N., W.N., W.Z. and J.G. acknowledge support from the Arts \& Sciences Fellowship of Washington University in St Louis. H.C. also acknoledge support from NASA grant 80NSSC24K1095. 
C.P. received support from Fundação para a Ciência e a Tecnologia (FCT), I.P., Portugal, under the  projects UIDB/04564/2020 (doi:10.54499/UIDB/04564/2020), UIDP/04564/2020 (doi:10.54499/UIDP/04564/2020), and 2022.06460.PTDC (doi:10.54499/2022.06460.PTDC).
L.T. acknowledges support from CEX2020-001058-M (Unidad de Excelencia ``Mar\'{\i}a de Maeztu") and PID2022-139427NB-I00 financed by the Spanish MCIN/AEI/10.13039/501100011033/FEDER,UE as well as from the Generalitat de Catalunya under contract 2021 SGR 171,  from the Generalitat Valenciana under contract CIPROM/2023/59 and by the CRC-TR 211 'Strong-interaction matter under extreme conditions'- project Nr. 315477589 - TRR 211. A.L.W. acknowledges support from ERC Consolidator Grant No.~865768 AEONS. A.C. acknowlege the support from NSF grants DMS-2235457 and AST-2308111.

# References