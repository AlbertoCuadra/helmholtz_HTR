# Helmholtz HTR
Routine to compute the dissipation, dissipation rate, and the solenoidal and compressive parts of a three-dimensional velocity field of a DNS obtained using the Hypersonic Task-based Research (HTR) solver [1] (see the GitHub repository [2]).
    
The calculation of the Helmholtz-Hodge decomposition using fast Fourier transform [3] is based on Ref. [4] and has been rewritten in MATLAB.

**Some remarks [4]:**
 * Only for uniform grid with $dx = dy = dz$.
 * Helmholtz-Hodge decomposition performed with the spectral method should only apply to relatively smooth fields, i.e., with little power on small scales.
 * For even $N_x$, $N_y$, and $N_z$, decomposed fields can be complex, with the imaginary part coming from the real part of the kmode at Nyquist frequency. In principle, the Nyquist frequency kmode should be dropped when doing the first derivatives to maintain symmetry. See footnote on page 4 of [2]. However, when the field is smooth enough, the imaginary part caused by the Nyquist frequency kmode should be negligible.

**References**
1. Di Renzo, M., Fu, L., and Urzay, J., HTR solver: An open-source exascale-oriented task-based multi-GPU high-order code for hypersonics aerothermodynamics, Comput. Phys. Commun, Vol. 255, 2020, p. 107262
2. Di Renzo, M., Fu, L., and Urzay, J., Hypersonic Task-based Research (HTR) solver for the Navier-Stokes equations at hypersonic Mach numbers including finite-rate chemistry for dissociating air and multicomponent transport, Github, 2022, Available: https://github.com/stanfordhpccenter/HTR-solver
3. Johnson, S. G., Notes on FFT-based differentiation. MIT Applied Mathematics, Tech. Rep., 2011, Available: http://math.mit.edu/~stevenj/fft-deriv.pdf
4. Xun, S., Helmholtz-Hodge decomposition using fft (Python), Github, 2018, Available: https://github.com/shixun22/helmholtz

## Acknowledgments
The Helmholtz-Hodge decomposition is based on the following repository: Xun Shi (2018). helmholtz (https://github.com/shixun22/helmholtz), Github. Retrieved March 15, 2023.
