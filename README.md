<!-- ABOUT THE PROJECT -->
## About The Project

Solves the Tolman-Oppenheimer-Volkoff equation for a spherically symmetric compact object.

### Inlcudes
1. Routine to check a valid equation of state input
2. Return the mass, radius, and tidal deformability, and compute the corresponding speed of sound.
3. [Illustration Notebook](https://github.com/ChunHuangPhy/EoS_inference/blob/main/Test_Case/tests_script.ipynb) on the github to show off what we can do currently and how to use our code.
5. Test cases and documentation


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
