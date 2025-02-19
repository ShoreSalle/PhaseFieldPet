# PhaseFieldPet
Open-Source Phase Field Simulation Software for Heterogeneous Architectures. It is built on top of [PETSc](https://petsc.org/release/), and hence wide variety of numerical solvers are available at run time. PhaseFieldPet models both Multiphase-field models ( phase field represents volume fraction and hence conserved, `pfe_mpfl` or `pfe_mpf`) and  Multi-order parameter models (also called continuum field models, where)order parameters are not conserved, `pfe_mop`).
In its current version (v1.0.0), PhaseFieldPet comes with three gradient energy terms (`grad_dot`, `grad_weighted`, `grad_interpolated`) and five potential energy terms (`pot_toth`, `pot_moelans`, `pot_garacke`,`pot,nestler`, `pot_steinbach`). The first three are well potentials, while the last two are obstacle potentials.  Not all combinations give a phsicaly sound simulation, so the user has to experiment with combinations and reason out why.

PhaseFieldPet runs on wide variety of hardware and is based on:
   - MPI Based (distributed computing capability by default).
   - GPUs ( NVIDIA, AMD): CUDA, HIP, Kokkos, openCL
   - openMP (using third party package integration with petsc)

# How to use

## Installation
Steps to install the project.
- See petsc.org
- Older version available from distro:
  ```bash
  sudo apt install petsc-dev
  ```

## Compile your code mpf.c
  ```bash
     make mpf
  ```

## Run The simulation
## On CPU with 24 process
  ```bash
mpiexec -n 24 mpf
  ```
This solves with default Time Step solver set in the code - Adaptive Runge Kutta ( Second order).
  ```bash
  mpiexec -n 24 mpf -ts_type bdf
  ```

 (This uses Implicit adaptive backward Differentiation Formula)
## On a GPU (cuda based)
  ```bash
mpiexec -n 1 mpf -ts_type bdf -dm_vec_type cuda -dm_mat_type aijcusparse
  ```
## Not today
  ```bash
mpiexec -n 24 mpf -ts_type bdf -sim mop
  ```

