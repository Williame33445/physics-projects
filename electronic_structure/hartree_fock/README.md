# Hartree-Fock project

This project uses the Hartree-Fock method to describe simple atoms and molecules. This project used (Thijssen, 2013, p. 43) as a guide, though it goes further than the problems specified in the book. 

## Theory

### Basic approach

The full Hamiltonian of a molecular system is given by $(1)$.

![alt text](https://raw.githubusercontent.com/Williame33445/physics-projects/main/electronic_structure/hartree_fock/electronic-hamiltonian.png)

One of the main problems in molecular physics is to find the eigenstates of $(1)$. It turns out that in most cases analytical solutions are not possible, meaning numerical methods need to be applied.

The Hartree-Fock method is a possible numerical method. It uses the Born-Oppenheimer approximation to decouple the electron and nuclei parts (Griffiths, 2018, p. 428),  guesses that the eigenstate is a Slater determinant (Thijssen, 2013, p. 53) and applies variational methods to find an upper bound of the eigenstates. This process reduces $(1)$ to $N$ equations in the form:

$$\hat{F}\psi_k(\vec{r},m_s) = \epsilon_k\psi_k(\vec{r},m_s)\tag{2}$$

where $\psi_k(\vec{r},m_s)$ is a component of the Slater determinant, $\hat{F}$ is called the Fock operator (derived in  (Thijssen, 2013, p. 56)) and $\epsilon_k$ is related to the energy of the total eigenstate.

$\hat{F}$ contains $\psi_k(\vec{r},m_s)$ terms in its definition, meaning $(2)$ forms a self-consistency equation. This means $(2)$ does not take the form of an eigen-equation, making it much harder to solve.


To solve (2), the variational method is applied again. An upper bound of the total electron energy can then be deduced from $(3)$.

$$E = \frac{1}{2}\sum\limits^N_{k=1}\epsilon_k + \bra{\psi_k}-\frac{\nabla^2}{2} + \sum_n\frac{Z_n}{|\vec{r}-\vec{R}_n|}\ket{\psi_k}\tag{3}$$

This entire process is the Hatree-Fock method.

### Variational method for the Fock Equation

For the second application of the variational method, the general form of $\psi_k(\vec{r},m_s)$ is guessed. For atomic systems, the most natural guess is a linear combination of Slater-type orbitals (STO). These take the form of $(4)$ (Slater, 1930, p.57).

$$
\psi_k(\vec{r},m_s) = \sum\limits_n A_n     r^{n-1}e^{-a_nr}Y_{lm}(\theta,\phi)\chi_{m_s}\tag{4}
$$

The coefficients $a_n,A_n$ are then varied to find the stationary $\epsilon_k$; this corresponds to an upper bound of the eigenstate.

However, due to integrals involving STO being hard to calculate, its usually easier to use Gaussian-type orbitals (GTO) (Goings, 2017). These take the form of $(5)$.

$$
\psi_k(x,y,z,m_s) = \sum\limits_{ijk} A_{ijk} x^iy^jz^ke^{-\alpha_{ijk}(x^2 + y^2 + z^2)}\chi_{m_s}\tag{5}
$$  

Again, the coefficients $a_{ijk},A_{ijk}$ are then varied to find the stationary $\epsilon_k$. 

Whatever basis is chosen, the effect is to transform the problem into a non-linear variational problem. While this can be solved directly, if the $\alpha_{ijk}$ terms are deduced by some other method, the problem becomes resolvable by (easier) linear methods. This is typically done by either fitting a set of GTOs to an appropriate STO or solving a simpler but appropriate non-linear problem (Thijssen, 2013, p. 67). This project doesn't cover these methods, it only considers the linear method. The $\alpha_{ijk}$ coefficients used are taken from (Thijssen, 2013, p. 35, p. 50) and (MolSSI, 2020).

Given the $\alpha_{ijk}$ terms, the variational problems can then be written in matrix form as:

$$F^+C_k^+ = \epsilon_k^+SC_k^+,F^-C_k^- = \epsilon_k^-SC_k^-\tag{5}$$

where $F_{nm}^{\pm}$ are the matrix elements for the plus/minus parts of the GTO, $S$ is the GTO overlap matrix and $C_k^{\pm}$ is the representation of $\psi_k(\vec{r},\pm\frac{\hbar}{2})$ in the GTO basis (Thijssen, 2013, p. 64).

This equation can then be solved iteratively. $F$ is deduced for an initial guess of the $\psi_k(\vec{r},m_s)$ 's and the generalised eigenvalue problem is solved. This process is then repeated, using the $\psi_k(\vec{r},m_s)$ terms found as the initial guess, until the solution converges. These solutions can then be used to deduce an upper bound of the eigenstate.

## Implementation

As any basis can be chosen for variational solutions to $(4)$, the matrix elements required to form equation $(5)$ are described in an abstract class called Representation. A subclass of Representation then corresponds to a certain basis. Two subclasses are defined: Rep1sGTO only uses the 1s Gaussian basis and RepGTO considers the general Gaussian basis. In principle, other bases could be added (eg. STO), however this is not implemented here. The matrix elements required for Rep1sGTO are defined in GTO1s_matrix_elements.py; for a derivation of these elements, see (Thijssen, 2013, p. 64). The matrix elements required for RepGTO are found with code taken from (Goings, 2017). These matrix elements are defined in MolecularIntegrals.py.

In terms of the actual algorithm, this is implemented in Hartee_Fock.py as a recursive function called iterateHF. The eigenstate that is being calculated is determined by the getTargetEigStates function that is a parameter of iterateHF.

## Results

Some simple examples of the simulation are given in hartree-fock-calculations.ipynb. Some of the ground state energies found are given in the table below.

|  Atom/Molecule   | HTF Energy (hartree) |Experimental Energy (hartree)|
| -------- | ------- | -----|
| $He$  | -2.855    | -2.862 |
| $Li$ | -7.410 | -7.479
| $O$    | -74.582    | -75.113

Molecular parameters (like bond length) can also be deduced by finding the parameter than minimises the total energy (this is effective as at this level of approximation molecular energy states seem continuous). In hartree-fock-calculations.ipynb, $H_2$ and $H_2O$ are considered. For $H_2$, the bond length found was $1.388$ au. The experimental result is 1.401 au. For $H_2O$, the bond length found was 1.817 au and the angle between bonds found was 107.227 degrees. The experimental values are 1.809 au and 104.5 degrees.
## References

Thijseen J., 2013. Computational Physics. Cambridge: Cambridge University Press

Griffiths D., 2018. Introduction To Quantum Mechanics.  Cambridge: Cambridge University Press

Slater J., 1930. Atomic Shielding Constants. Physical Review: Volume 36

Goings, J., 2017. joshuagoings. [Online] 
Available at: https://joshuagoings.com/2017/04/28/integrals/
[Accessed 17 June 2024].

The Molecular Sciences Software Institute (MolSSI), V. T., 2020. Basis Set Exchange. [Online] 
Available at: https://www.basissetexchange.org/
[Accessed 17 June 2024].

