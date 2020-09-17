# Object Oriented Trajtou Program

## Acknowledgement and citation
The current code is based on the works of F. Busnengo et al. [1] and subsequent uses of the Corrugation Reducing Procedure (CRP) method in gas-surface dynamics available in the bibliography.
This implementation adds a logistic function that controls the "amount" of CRP in the interpolation as a function of the distance of the projectile to the surface. If this code is used in any publishable study, it should include the reference: "A. S. Muzas, M. del Cueto, F. Gatti, M. F. Somers, G. J. Kroes, F. Martín and C. Diaz, Phys. Rev. B 96 205432 (2017)". That was the first paper in which our particular modified CRP method was described in detail.

## Brief functionality
1. Program for **PES** interpolation using the **CRP method**.
2. It can carry out **classical hamiltonian dynamics simulations**, adapted to scattering processes on **periodic surfaces**.
3. As the CRP method is only applicable to **monoatomic and diatomic projectiles**, the dynamics integrator has been implemented only for these two cases.

## Installation
1. Execute install.sh script and follow the instructions. You may edit that script to add new compiler compatibilities.

## References
[1] H. F. Busnengo, A. Salin and W. Dong, J. Chem. Phys. 112.17 (May 2000)
