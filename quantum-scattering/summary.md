# Quantum Scattering Project

This project, adapted from (Thijssen, 2013, p. 14), simulates molecule scattering experiments using the Lennard-Jones potential (Fischer, 2023, p. 2) to model interactions. These simulations are used to produce the differential cross section and total cross section.

## Theory

### Basic idea

The scattering is described by the time-independent Schrödinger equation. This is:

$$-\frac{\hbar^2}{2m}\nabla^2\psi(\vec{r}) + V_{LJ}(r)\psi(\vec{r}) = E\psi(\vec{r})\tag{1}$$

where $\vec{r}$ is the relative position of the molecules, $V_{LJ}(r)$ is the Lennard-Jones potential, $m$ is the reduced mass of the molecules, $\psi(\vec{r})$ is the wavefunction, $E$ is the total energy of the molecules[^1] (can be any value) and $\hbar$ is Plank's constant.

Due to spherical symmetry, $\psi(\vec{r})$ can be expanded in terms of spherical harmonics as:

 $$\psi(\vec{r}) = \frac{1}{r}\sum\limits_{l = 0}^{\infty} \sum\limits_{m=-l}^{l}A_{l,m}u_l(r)Y_{l,m}(\theta,\phi)\tag{2} $$

Plugging this in to $(1)$ gives the radial time-independent Schrödinger equation for each $u_l(\vec{r})$:

$$\left[\frac{\hbar^2}{2m}\frac{d^2}{dr^2} + E - V(r) - \frac{\hbar^2l(l+1)}{2mr^2}\right]u_l(r) = 0\tag{3}$$

Equation $(3)$ cannot be easily solved analytically for cases where $V_{LJ}(r)\neq 0$. However, as $V_{LJ}(r)$ is only approximately non-zero when $r$ is small, the system can be solved numerically in this region. Then, using the boundary conditions generated from the numerical part, solutions in the $V_{LJ}(r)= 0$ region can be deduced. This solution can then be used to find simulation parameters. This is the basic idea behind this project.

### Numerical simulation

The numerical simulation in the $V_{LJ}(r)\neq 0$ region is performed by the Numerov method (Thijssen, 2013, p. 573). Using an appropriate set of initial conditions, it finds a set of $u_l$ values for a given set of $r$ values. Theoretically, all $l$ values would have to be solved, however only the low $l$ parts of the solution have a non-negligible effect[^2]. This means only $l<9$ solutions to (3) need to be considered.

### Finding scattering parameters

To equate the numerical results, $V_{LJ}(r)=0$ solutions need to be found and appropriate boundary conditions applied. Solutions for $V_{LJ}(r)=0$ are:

$$u_l(r) \propto kr(\cos(\delta_l)j_l(kr) - \sin(\delta_l)n_l(kr))\tag{4}$$

where $k = \frac{\sqrt{2mE}}{\hbar}$, $j_l,n_l$ are spherical Bessel functions and $\delta_l$ is a real number (called phase shift).

From this, given $u_l(r_1)$ and $u_l(r_2)$ where $V_{LJ}(r_1)\approx V_{LJ}(r_1)\approx 0$, it can be shown that:

$$\tan(\delta_l) = \frac{K j_l(r_1) - j_l(r_2)}{K n_l(r_1) - n_l(r_2)}\tag{5}$$

where $K = \frac{r_1 u_l(r_2)}{r_1 u_l(r_2)}$.

Equation (5) is a convenient way to deduce the $\delta_l$ s from an appropriate part numerical simulation data. Conveniently, it turns out the $\delta_l$ contains all the information required to calculate the scattering parameters. This can be shown by the form of the scattering parameters for the system (for a derivation see (Thijssen, 2013, p. 16)), which are given by[^3]:

$$\frac{d\sigma}{d\omega} = \frac{1}{k^2}\bigg|\sum^{\infty}_{l = 0}(2l+1)e^{i\delta_l}sin^2(\delta_l)P_l(\cos(\theta))\bigg|^2\tag{6}$$

$$\sigma_{total}= \frac{4\pi}{k^2}\sum^{\infty}_{l = 0}(2l+1)sin^2(\delta_l)\tag{7}$$

where $P_l(x)$ is the $l^{th}$ Legendre polynomial, $\frac{d\sigma}{d\omega}$ is the differential cross section and $\sigma_{total}$ is the total cross section. As with the simulations, the terms of the sums become negligible for $l>9$, meaning these terms do not need to be evaluated.

### Summary of method

In summary, simulation parameters are deduced by numerically solving (3) in the $V_{LJ}(r) \neq 0$ region for $l<9$, finding the corresponding $\delta_l 
$ s via equation (5) and then plugging into equations (6) and (7) to deduce simulation parameters.

## Implementation

This project was implemented as a class called ScatteringSystem in spherically_symetric_quantum_scattering.py. This class takes in the initial conditions near the origin  (r_0,h,u_0,u_1)[^4], the energy of the system (E), the radius of the $V(r)\neq 0$ region (r_b1) and the max value of $l$ to calculate to (lCutoff[^5]). The class then runs the numerical simulation, deduces the $\delta_l$ s and finds the total cross section (stored as class.totalCrossSection). The differential cross section can then be found using the method differentialCrossSection().

## Example

In hydrogen-scattering-by-krypton-atoms.py the project is used to deduce the total cross section at different energies for scattering of hydrogen with krypton atoms. The distribution found is:

![alt text](https://raw.githubusercontent.com/Williame33445/physics-projects/main/quantum-scattering/total-cross-section-vs-energy-for-H-Kr-scattering.png)

This corresponds to experimental results.

## References

Thijseen J., 2013. Computational Physics. Cambridge: Cambridge University Press

Fischer, J.; Wendland, M., 2023, On the history of key empirical intermolecular potentials, Fluid Phase Equilibria. 573: 113876




[^1]: It is conventional to consider the frame where the heavier molecule is stationary at the origin and the lighter particle is moving; in this case, E is the lighter particles kinetic energy at infinity. This convention is adopted in this project.

[^2]: For a theoretical justification see (Thijssen, 2013, p. 18)

[^3]: There are two common types of scattering parameter: differential cross section and total cross section. Differential cross section is a function of theta and phi. It is proportional to the probability of a particle scattering in that direction. Total cross section is a number. It is proportional to the probability of a particle being scattered.

[^4]: Instead of r_0 and r_1, r_0 and h = r_1 - r_0 is given.

[^5]: In this situation only l<9 is needed.
