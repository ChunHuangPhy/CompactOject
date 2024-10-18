<!-- ABOUT THE PROJECT -->
## About The Project

1. Dealing with complex Relativistic Mean field (RMF) theory to generate Equation of State (EOS) of neutron star. ([EOSgenerators](https://github.com/ChunHuangPhy/EoS_inference/blob/main/EOSgenerators) Package)
2. Solves the Tolman-Oppenheimer-Volkoff equation for a spherically symmetric compact object out of given equation of state of neutron star. ([TOVsolver](https://github.com/ChunHuangPhy/EoS_inference/blob/main/TOVsolver) Package)
3. Implementing Neutron state EOS inference by Nested Sampling, draw constraints from Nuclear experiments, Neutron star mass (and/or) radius observations (from X-ray timing and/or radio timing) (and/or) Tidal measurement from Gravitational wave detection. That all workflow is inside this folder. ([InferenceWorkflow](https://github.com/ChunHuangPhy/EoS_inference/blob/main/InferenceWorkflow) Package) 


Project papers list based these package: (Please consider cite them, if you are using this package)

[1]. [Huang, C., Raaijmakers, G., Watts, A. L., Tolos, L., and Providência, C., “Constraining fundamental nuclear physics parameters using neutron star mass-radius measurements I: Nucleonic models”,Monthly Notices of the Royal Astronomical Society,2024, 10.1093/mnras/stae844,529, https://academic.oup.com/mnras/article/529/4/4650/7634362

## Sample citation line:
```diff
"The inference conducted here relies on the framework in \textit{CompactObject} \cite{CompactObject} package\footnote{https://chunhuangphy.github.io/CompactOject/}. This is an open source full-scope package  designed to implement Bayesian constraints on the neutron star EOS. The other work based on this package is ...."

1.https://chunhuangphy.github.io/CompactOject/
```

[CompactObject-TOV package website](https://chunhuangphy.github.io/CompactOject/)

### Inlcudes
1. Routine to check a valid equation of state input
2. Return the mass, radius, and tidal deformability, and compute the corresponding speed of sound.
3. [Sample TOV solver Notebook](https://github.com/ChunHuangPhy/EoS_inference/blob/main/Test_Case/test_TOVsolver.ipynb), [Sample RMF Equation of state solver Notebook](https://github.com/ChunHuangPhy/EoS_inference/blob/main/Test_Case/test_EOSgenerators.ipynb) and [Sample Analysis Notebook on Equation of state Inference and tutorial](https://github.com/ChunHuangPhy/EoS_inference/blob/main/Test_Case/test_Inference.ipynb) on the github to show off what we can do currently and how to use our code. (**please read them before you start to work on your own project, to familiar with the coding routine.**)
4. Test cases and documentation
### v.1.3 new features:
5. Added computation function of generating Relativistic mean field theory(RMF) model EOS functionality. Defined two files fastRMF_EOS and RMF_EOS, which the fastRMF_EOS is speed up by numba, which need gcc compiler, could be hard to implement in windows, so we leave the options for users.
### v.1.5 new features:
6. Added Whole workflow of Bayesian inference of neutron star equation of state. Include defining prior by InferenceWorkflow.prior, which included two types: flat distribution and gaussian type. Include defining liklihood generated from nuclear and astrophysical constraint.


## Installation

_Below are commands to install and update the package as well as a link to pypi._


##### [PyPi](https://pypi.org/project/CompactObject-TOV/)



1. Install package
   ```sh
   pip install CompactObject-TOV
   ```
2. Update package
   ```sh
   pip install CompactObject-TOV --upgrade
   ```

When you call the package, if you need to do EoS computation just
   ```sh
   import EOSgenerators
   ```
if you need TOV solver, just
   ```sh
   import TOVsolver
   ```
if you need to do Bayesian inference, just
   ```sh
   import InferenceWorkflow
   ```

## Physics notations
1. CGS units is using here, for input quantity (equation of state): Pressure (P) and Energy density (rho).
P is in $MeV/fm^{-3}$, same for rho. However, to omit a lot of the repeat of c,G. We set P as rescaled:
(value in $MeV/fm^{-3}$)*G/c^4, for rho we have (value in $MeV/fm^{-3}$)*G/c^2
2. Out put M in Mass of sun, radius in km, unit-less for spped of sound and tidal deformability.
<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- CONTACT -->
## Contact

* Chun Huang - chun.h@wustl.edu
* Nicole Osborn - n.osborn@wustl.edu
* Nathan Whitsett - whitsett.n@wustl.edu

Project Link: [[https://github.com/ChunHuangPhy/EoS_inference](https://github.com/ChunHuangPhy/EoS_inference)]



[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.8145167.svg)](http://dx.doi.org/10.5281/zenodo..8145167)


<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

Use this space to list resources you find helpful and would like to give credit to. Here included a few of my favorites to kick things off! We would like to acknowledge the support of Code/Astro workshop to make this project happen, we all learned a lot from that 5 day intensive workshop.

Chun want to thank Professor Anna Watts, Dr. Geert Raaijmakers and Jeannie Kuijper for asistance on coding and help on providing me basic strategy of how to solve this problem.

* [Code Astro](https://github.com/semaphoreP/codeastro)
* [Choose an Open Source License](https://choosealicense.com)
* [Read Me Template](https://github.com/othneildrew/Best-README-Template)

<p align="right">(<a href="#readme-top">back to top</a>)</p>
