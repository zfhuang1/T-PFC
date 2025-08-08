# T-PFC

This is a Fortran 2003 source code for numerically solving the 2D T-PFC model equation (the two-dimensional phase field crystal model incorporating transverse interactions), which corresponds to Eq. 3 in [Anomalous grain dynamics and grain locomotion of odd crystals, arXiv2505.03957 (2025)](https://doi.org/10.48550/arXiv.2505.03957).

## **Usage**
The main source code is in T-PFC2D.F03, while ran2.f90 contains a function of random number generator given in the book "Numerical Recipes". To use and compile the codes, the FFTW3 library is needed, as well as a file fftw3.f03 for the related Fortran 2003 interfaces (as provided by FFTW; see the FFTW manual at [fftw.org](https://fftw.org/)).

## **Contributing**
This code was developed by [Zhi-Feng Huang](https://s.wayne.edu/huang/).

## **Citation**
If using the code in research, please cite:\
Z.-F. Huang, M. te Vrugt, R. Wittkowski, and H. Löwen, "Anomalous grain dynamics and grain locomotion of odd crystals", [arXiv2505.03957](https://doi.org/10.48550/arXiv.2505.03957) (2025).
