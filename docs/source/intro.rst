.. _readme:

***********************************
CompactObject package Tutorials
***********************************

**An open-source package for neutron star**
**whole workflow Bayesian nference constraining**
**Neutron star EOS package**

CompactObject is designed to open-source a full-scope neutron star equation inference.
Currently, this equation of state framework is mainly based on relativistic mean field
(RMF) theory. The package integrates three independent modules:

1. A user-friendly Tolman–Oppenheimer–Volkoff (TOV) equation solver for determining
neutron star structure based on a given equation of state (EOS).

2. A neutron star full density range equation of state generator. At present, it 
only contains the relativistic mean field theory equation of state solver. In 
the future, it will easily accommodate polytropes and more equation of states.

3. A complete Bayesian inference workflow package for constraining the equation 
of state of neutron stars. This includes defining the likelihood of inference from
observations and nuclear experiments, as well as from simulated astrophysical 
observations. It also involves defining priors for the given parameters and running 
a nested sampling of the posterior space.


We should mention that these three parts are independent, which means users can embed
these modules as part of their work. The usage of this package extends beyond conducting
inference studies. Additionally, we welcome contributions of new features for our package.


The papers generated from this package are listed below:
[1]
`(Huang, C., Raaijmakers, G., Watts, A. L., Tolos, L., and Providência, C.,
“Constraining fundamental nuclear physics parameters using neutron star mass-radius
measurements I: Nucleonic models”, 2023. doi:10.48550/arXiv.2303.17518) <https://arxiv.org/abs/2303.17518>`_

If you are using our software, please consider citing us using the following standard citation:

.. code-block:: "The inference conducted here relies on the framework in \textit{CompactObject} \cite{CompactObject} package\footnote{\url{https://chunhuangphy.github.io/CompactObject/}}. This is an open-source, full-scope package designed to implement Bayesian constraints on the neutron star EOS. Other work based on this package is ...."

Concept
*******

Bayesian inference studies of the neutron star equation of state have become a trending
field nowadays, particularly due to significant advancements such as the Neutron Star 
Interior Composition Explorer (NICER) measuring the mass and radius of neutron stars 
through X-ray timing, as well as the detection of neutron star merger events through
gravitational wave observations by the LIGO detector.

Here below is the overall pipline of this field:

.. image:: workflow.png


As depicted in this plot, fundamental physics can provide the equation of state for 
neutron stars. By inputting the equation of state into the Tolman–Oppenheimer–Volkoff
(TOV) equation, we can obtain parameters related to the neutron star structure, such
as mass, radius, and tidal properties. Instruments like NICER and LIGO can be utilized
to measure these properties. With the information obtained from observations, it
becomes possible to perform Bayesian inference constraining, such as determining the
region in the Mass-Radius space where neutron stars can exist. Subsequently, we can
derive constraints on the equation of state space and, ideally, obtain insights into
the fundamental composition of the neutron star's interior.

**Equation of State**

The neutron star equation of state (EOS) plays a crucial role in reflecting the 
composition of a neutron star. It is closely connected to the microphysical 
properties of neutron stars. In this context, we are utilizing an equation of 
state derived from a model known as the Relativistic Mean Field theory (RMF). 
The Lagrangian of this model is represented as follows:

.. math::

   \mathcal{L}=\sum_N \mathcal{L}_N+\mathcal{L}_{\mathcal{M}}+\sum_l \mathcal{L}_l


Where the :math:`\mathcal{L}_N` is the neucleonic Lagrangian, :math:`\mathcal{L}_M`
is the meson part Lagrangian, :math:`\mathcal{L}_l` is the lepton Lagrangian.
details of the Lagrangian are

.. image:: lagrangian.png

where :math:`\Psi_{N}` and :math:`\psi_{l}` are the nucleon and lepton spinors,
and :math:`\bar{I}_{N}` is the nucleon isospin operator. The strong interaction
coupling of a meson to a nucleon is denoted by :math:`g`, while the masses of 
the nucleons, mesons, and leptons are denoted by :math:`m`. The parameters :math:`\kappa`,
:math:`\lambda_0`, :math:`\zeta` and :math:`\Lambda_{\omega}` plus the meson-nucleon
coupling constants are coupling constants to be determined by the inference method.

These free parameters represent the degrees of freedom in the RMF model and can be
determined through nuclear experiments. However, in addition to nuclear experiments, 
we can also explore the possibility of constraining these parameters through 
astrophysical observations. The complete list of parameters includes:

.. image:: free_para.png

These are the parameters that you should input to generate the equation of state 
from our EOSgenerators module, different equation of state parameter will have different
effect on mass radius like we showed here.More details about this physics can check 
`(Glendenning, 1996) <https://ui.adsabs.harvard.edu/abs/1996cost.book.....G/abstract>`_

**Tolman–Oppenheimer–Volkoff(TOV) equation**

TOV equation is a general relativity equation that constrain the structure of 
a spherical symmetrical body by gravity. This is the original equation:

.. math::

    \frac{d P}{d r}=-\frac{G m}{r^2} \rho\left(1+\frac{P}{\rho c^2}\right)\left(1+\frac{4 \pi r^3 P}{m c^2}\right)\left(1-\frac{2 G m}{r c^2}\right)^{-1}



To solve this problem, the essential ingredient is the equation of state (EOS).
Once you have the EOS, the basic strategy for solving the equation is as follows: 
at a given central density, you input it into the neutron star EOS to obtain the 
pressure. Then, you integrate the density from the center to the boundary, repeating 
this process across the entire possible density range.

In our code, we provide two different functions that you can call. The default option
allows you to solve the Tolman-Oppenheimer-Volkoff (TOV) equation within a predefined 
density range, which is log(14.3, 15.6) on a scale of 10. Alternatively, you can choose 
to solve the equation point by point, allowing you to select any central density range
you prefer. More information about TOV you could check 
`(wiki page) <https://en.wikipedia.org/wiki/Tolman–Oppenheimer–Volkoff_equation>`_

The following image integrates (and thus averages) over waveband (a
range of photon energies). We also decrease the mode frequency relative to the
stellar spin frequency, such that the mode is not as equatorially trapped.

**Bayesian Inference**

Using Bayesian inference tools to explore the constraint of neutron star equation of state
is common nowadays, the basic equation of it is Bayes theorem:


.. math::
    P(A \mid B)=\frac{P(B \mid A) P(A)}{P(B)}

That is, Posterior probablity is propotional to the prior probablity times likelihood.
Posterior is --- after correction of the new observations/experiment, the probablity of
something is true.
Prior is --- before the new observations/experiment come in, my initial thought about the
probablity of something is true.
Likelihood is --- the correction that we get from the new observation/experiments.

Here, the likelihood will be mostly come from three different families:

1. Mass Raius measurements from x-ray timing (like NICER).
2. Tidal measurements from gravitational wave detection
3. Mass measurements from radio timing.
4. Nuclear physics constraint comes from the nuclear experiments.

NICER Mass radius measurements are remarkable achievement of this centry of neutron star
physics, same as the gravitational wave detection. many references out there for this topic.
About the Nuclear physics connection between our equation of state and the nuclear quantities,
please check `(Chen & Piekarewicz 2014a) <https://journals.aps.org/prc/abstract/10.1103/PhysRevC.90.044305>`_


Here the nuclear physics quantities we cared are K, J and L, that is the decompressibility of
nuclear matter K, symmerty energy at saturation density J, and the slope of symmetry energy
at saturation density L. These all can be computed out by posterior samples (will add the nuclear
properties computation code in near future). Also they could be independent group of constraint on 
our equation of state of neutron star.

When you do a Real astrophysical sampling, the important thing is you should also sampling the 
neutron star central density of that measurement you are using, which means if you want to investigate
what the constraining effect for neutron star EOS by two mass radius measurements, then you need 
define another two free parameters ---  the central densities of these measurements, other-wise, 
this could be proved to be a not full-scope equation of state inference, that is why our likelihood
functions once you want to constraint from observation, always need a parameter d1, that is the 
density parameter of this observation. 
