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
    - name: Wenli Yuan
      orcid: 0000-0003-2771-759X
      affiliation: 3
    - name: Tianzhe Zhou
      orcid: 
      affiliation: 4
    - name: Xuezhi Liu
      orcid: 
      affiliation: 5
    - name: Xieyuan Dong
      orcid: 
      affiliation: 6
    - name: John Groger
      orcid: 0000-0002-7054-9053
      affiliation: 1
    - name: Constan\c{c}a Provid\^{e}ncia
      orcid: 0000-0001-6464-8023
      affiliation: 2
    - name: Micaela Oertel
      orcid: 0000-0002-1884-8654
      affiliation: 7
    - name:  Laura Tolos
      orcid: 0000-0002-6449-106X
      affiliation: 8
      affiliation: 9
      affiliation: 10
    - name: Alexander Y. Chen
      orcid: 0000-0002-4738-1168
      affiliation: 1
    - name: Anna L. Watts
      orcid: 0000-0002-1009-2354
      affiliation: 11
affiliations:
   - name: Physics Department and McDonnell Center for the Space Sciences, Washington University in St. Louis; MO, 63130, USA
     index: 1
   - name: CFisUC, Department of Physics, University of Coimbra, 3004-516 Coimbra, Portugal
     index: 2
   - name: School of Physics and State Key Laboratory of Nuclear Physics and Technology, Peking University, Beijing 100871, China;
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
   - name: Frankfurt Institute for Advanced Studies, Ruth-Moufang-Str. 1, 60438, Frankfurt am Main, Germany
     index: 10
   - name: Anton Pannekoek Institute for Astronomy, University of Amsterdam, Science Park 904, 1090 GE Amsterdam, the Netherlands
     index: 11
date: Oct 23 2021
bibliography: cojoss.bib
---


# Summary

The CompactObject package is a software package designed to open source the pipeline on constraining
neutron star equation of state(EOS) by performing Bayesian statistical inference using Astrophysical
observation constraints from X-ray timing, Gravitational wave events and Radio measurements, and Nuclear
experiment constraint from pQCD and Chiral Effective field Theory ($\chi$EFT). Facilitated EOS models span
from meta-models to several physics-motivated models. Three independent parts are included: EOS generator, 
Tolman–Oppenheimer–Volkoff (TOV) equation solver, and the full workflow Bayesian inference module. Each module 
could be independently implemented in various scientific reaserch. Synergy with exsisting software like CompOSE,
CompactObject can use the CompOSE EOS database to enrich the EOS choices.

# Statement of need

Neutron star interior
Pulsed X-ray signals from neutron stars
can be modeled to statistically estimate parameters such as stellar mass and
radius, and properties of the surface radiation field such as a map of
temperature. The mass and radius of a neutron star are a function of the
equation of state of internal matter, especially the dense matter in the core,
and the formation history of the star, which determines the central energy
density and the spin frequency. The state of the surface radiation field is
the product of a potentially long and complex stellar evolutionary history,
especially that of the stellar magnetosphere. Such parameter estimation
requires relativistic tracing of radiation as it propagates from surface to a
distant telescope. Pulse-profile modelling to infer neutron star parameters
is a major science goal for both current X-ray telescopes such as the Neutron
Star Interior Composition ExploreR [NICER, @Gendreau2016] and proposed future telescopes such as eXTP and STROBE-X
[@Watts2019;@Ray2019].

While there are some open-source libraries for simulating the X-ray
signals from rapidly spinning neutron stars and more generally from the
vicinity of general relativistic compact objects including black holes [@Nattila:2016;@Pihajoki:2018] the scope
of these projects does not include statistical modeling, which
necessitates tractable parametrised models and a modular framework for
constructing those models.  X-PSI addresses this need, coupling code for likelihood
functionality (simulation) with existing open-source software for posterior sampling (inference).

# The X-PSI package and science use

X-PSI is an open-source Python package for Bayesian modeling of time- and
energy-resolved X-ray pulsations. X-PSI provides a framework for the
implementation of custom models, including likelihood and prior functions, and for
feeding those models to open-source statistical sampling software for use on
high-performance computing systems. X-PSI supplies modules for post-processing
posterior sample sets, and supplies tools for visualisation. For example, one
can generate time- and energy-resolved images and animations of model X-ray pulsars
by tracing radiation from the stellar surface to a distant observer, such as a
space telescope, as it propagates through spacetime; a snapshot of such an
animation may be found in \autoref{fig:animation_snapshot}. Posterior summaries
can be plotted of the time- and energy-domain signals that a model pulsar is
inferred to generate, conditioned on observational data.

![A snapshot from a time- and energy-resolved animation of a toy neutron star
that generates X-ray pulsations. The top three panels are sky maps of photon
specific intensity at three representative photon energies, and the bottom
panels show instantaneous integrals over the sky maps, together with sky map
integrals at past times normalized to the maximum at each photon energy
(bottom-left panel) and additional photon energies (bottom-right
panel).\label{fig:animation_snapshot}](fig1.png){width=100%}

X-PSI is being used by the NICER collaboration for pulse-profile modeling of X-ray emission from rotation-powered
millisecond pulsars [@Riley:2019;@Riley:2021].  Many more papers have used the accompanying
open-source analysis pipeline and products published on Zenodo [@Riley:2019:Zenodo;
@Riley:2021:Zenodo]. The first were @Raaijmakers:2019 and
@Bilous:2019, respectively on the topics of dense matter inference and
multipolar magnetic fields.

The numerical likelihood routines native to X-PSI are written in Cython
[@cython2011], and are dependent on the GNU Scientific Library [GSL,
@Gough:2009]. High-level object-oriented model construction is performed by a
user in the Python language, as is the interfacing with sampling software.
Low-level customisation is encouraged in the extensions, either directly in
Cython or via calls to external C libraries.  X-PSI is Unix source code
compatible, and release versions are freely available on GitHub under the GNU General Public License.  Extensive documentation, step-by-step tutorials, and reproduction
code for existing data analyses, are available
via the GitHub repository, along with a growing suite of unit tests.  Future plans
include migration to Python 3, further improvements to post-processing software,
 and the implementation of an expanded suite of atmosphere models.



*Software:* Python/C language [@Python2007], GNU Scientific Library [GSL,
@Gough:2009], NumPy [@Numpy2011], Cython [@cython2011], OpenMP [@openmp], MPI
for Python [@mpi4py], Matplotlib [@Hunter:2007; @matplotlibv2], IPython
[@IPython2007], Jupyter [@Kluyver:2016aa], MultiNest [@MultiNest_2009],
PyMultiNest [@PyMultiNest], GetDist [@Lewis19], nestcheck
[@higson2018nestcheck;@higson2018sampling;@higson2019diagnostic], fgivenx
[@fgivenx], emcee [@emcee]. 

# Acknowledgements

All University of Amsterdam co-authors acknowledge
support from ERC Consolidator grant No. 865768 AEONS (PI: ALW).  DH is supported by the 
Women In Science Excel (WISE) programme of the Netherlands Organisation for 
Scientific Research (NWO). SG acknowledges the support of the CNES. More detailed acknowledgements are written in the project
documentation [hosted on GitHub](https://xpsi-group.github.io/xpsi/acknowledgements.html).

# References