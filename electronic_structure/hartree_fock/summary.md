# Hartree-Fock project

This project uses the Hartree-Fock method to describe simple atoms and molecules. This project used (Thijssen, 2013, p. 43) as a guide, though note that the book does not give guidence for results beyond Helium. 

## Theory

### Basic approach

The full Hamiltonian of a molecular system is given by $(1)$.

![alt text](https://raw.githubusercontent.com/Williame33445/physics-projects/main/electronic_structure/hartree_fock/electronic-hamiltonian.png)

One of the main problems in molecular physics is find the eigenstates of $(1)$. It turns out that analytical solutions are not possible in most cases, meaning numerical methods need to be applied. 

The Hartree-Fock method uses the Born-Oppenheimer approximation to decouple the electron and nuclei parts (Griffiths, 2018, p. 428). It then guesses that the eigenstate is a Slater determinant (Thijssen, 2013, p. 53) and uses the variational princple to find an upper bound of the eigenstate. The eigenequation then reduces to $N$ equations in the form of $(2)$ where $\psi_k(\vec{r},m_s)$ is a component of the Slater determinant ($m_s$ is the spin component in the z-direction), $\hat{F}$ is called the Fock operator (derived in  (Thijssen, 2013, p. 56)) and $\epsilon_k$ is related to the energy of the total eigenstate.

$$\hat{F}\psi_k(\vec{r},m_s) = \epsilon_k\psi_k(\vec{r},m_s)\tag{2}$$

As $\hat{F}$ contians $\psi_k(\vec{r},m_s)$, $(2)$ forms a self-consistency equation. Meaning $(2)$ doesn't take the form of an eigenequation, making it much harder to solve.


To solve (2), the variational method is applied again. An upper bound of the total electron energy can then be deduced from $(3)$.

$$E = \frac{1}{2}\sum\limits^N_{k=1}\epsilon_k + \bra{\psi_k}-\frac{\nabla^2}{2} + \sum_n\frac{Z_n}{|\vec{r}-\vec{R}_n|}\ket{\psi_k}\tag{3}$$

### Variational method

The variational method guesses the general form of $\psi_k(\vec{r},m_s)$. For atomic systems, the most natural guess is a linear combination of Slater-type orbitals (STO). These take the form of $(4)$ (Slater, 1930, p.57).

$$
\psi_k(\vec{r},m_s) = \sum\limits_n A_n     r^{n-1}e^{-a_nr}Y_{lm}(\theta,\phi)\chi_{m_s}\tag{4}
$$

The coefficients $a_n,A_n$ are then varied to find the stationary $\epsilon_k$; this corresponds to an upper bound of the eigenstate.

However, due to integrals involving STO being hard to calcualte, its easier to use Guassian-type orbitals (GTO) (Goings, 2017). These take the form of $(5)$.

$$
\psi_k(x,y,z,m_s) = \sum\limits_{ijk} A_{ijk} x^iy^jz^ke^{-\alpha_{ijk}(x^2 + y^2 + z^2)}\chi_{m_s}\tag{4}
$$  

Again, the coefficients $a_n,A_n$ are then varied to find the stationary $\epsilon_k$. This is a non-linear problem.


While the non-linear problem can be solved, if the $\alpha_{ijk}$ s are deduced, the problem becomes solvable by linear methods. This is typically done by either fitting a set of GTO to approriate STO or solving a simpler but appropriate non-linear problem (Thijssen, 2013, p. 67). These methods are not covered here and the $\alpha_n$ coefficents used are taken from (Thijssen, 2013, p. 35, p. 50) and (MolSSI, 2020).

Given the $\alpha_{ijk}$ s, the variational problems can then be written in matrix form as $(5)$, where $F_{nm}^{\pm}$ are the matrix elements for the plus/minus parts of the GTO, $S$ is the GTO overlap matrix and $C_k^{\pm}$ is the representation of $\psi_k(\vec{r},\pm\frac{\hbar}{2})$ in the GTO basis (Thijssen, 2013, p. 64).

$$F^+C_k^+ = \epsilon_k^+SC_k^+,F^-C_k^- = \epsilon_k^-SC_k^-\tag{5}$$

This equation can then be solved iteratively. $F$ is deduced for an intial guess of the $\psi_k(\vec{r},m_s)$ 's and the generalised eigenvalue problem is solved. This process is then repreated using the $\psi_k(\vec{r},m_s)$ 's found as the intial guess until the solution converges. These solutions can then be used to deduce an upper bound of the eigenstate.

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

