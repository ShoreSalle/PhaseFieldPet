# PhaseFieldPet
Open-Source Phase Field Simulation Software for Heterogeneous Architectures. It is built on top of [PETSc](https://petsc.org/release/), and hence wide variety of numerical solvers are available at run time. PhaseFieldPet models both Multiphase-field models ( phase field represents volume fraction and hence conserved, `pfe_mpfl` or `pfe_mpf`) and  Multi-order parameter models (also called continuum field models, where)order parameters are not conserved, `pfe_mop`).
In its current version (v1.0.0), PhaseFieldPet comes with three gradient energy terms (`grad_dot`, `grad_weighted`, `grad_interpolated`) and five potential energy terms (`pot_toth`, `pot_moelans`, `pot_garacke`,`pot,nestler`, `pot_steinbach`). The first three are well potentials, while the last two are obstacle potentials.  Not all combinations give a phsicaly sound simulation, so the user has to experiment with combinations and reason out why. For explanations of these terms and the phase field equation as presented in PhaseFieldPet we encourage readers to see a paper [Daubner et al., (2023)](https://doi.org/10.1016/j.commatsci.2022.111995) and [Chota et al., 2025 (TBD)](https://joss.theoj.org/papers/TBD).

PhaseFieldPet runs on wide variety of hardware and it supports:
   - MPI Based (distributed computing capability by default).
   - GPUs ( NVIDIA, AMD): CUDA, HIP, Kokkos, openCL
   - openMP (using third party package integration with petsc).
   - Pthread.

# How to use
  One can take  source code PhaseFieldPet.c, and compile and run it, visualize the results (default in vtk format) using visualizations softwares such as [ParaView](https://www.paraview.org/). The steps to download and install varies, but we give a general directions to do so here.

## Installation
There are many ways to install PETSc. See [PETSc Installation](https://petsc.org/release/install/).
#### From Linux Distribution Repository ( older version )
On debian based such as Ubuntu Linux for instance
  ```bash
  sudo apt install petsc-dev
  ```
#### From Source ( most updated versions )
Steps to install vary based on what softwares are available in the machine you are connected to. For instance on the system where gcc, g++, gfortran  are available , but not mpi implmentation or lapack, PETSc will install it for you following
```bash
  - git clone -b release https://gitlab.com/petsc/petsc.git petsc
  - cd petsc
  - ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-mpich --download-fblaslapack
  - make all check
```
If in addition you have NVIDIA GPU with compute capability, add `--with-cuda`  in `./configure` above. Do similarly to kokkos , OpenMP or Pthread.\
If your system has mpi installed ( or available, say via `module load openmpi` in HPC machines), find the directory where it is installed ( eg `which mpiexec`), by optionally adding optimization flags, do the following to configure PETSc
```bash
- ./configure PETSC_ARCH=arch-optimized --with-debugging=0  COPTFLAGS='-O3 -march=native -mtune=native'  CXXOPTFLAGS='-O3  -mtune=native'  FOPTFLAGS='-O3 -march=native -mtune=native'  --download-fblaslapack --with-mpi-dir=/Path/to/your/MPI/Dir
```
 Note that  if you installed mpi implmentation with PETSc `mpiexec` or `mpirun` will be available in directory `$PETSC_DIR/$PETSC_ARCH/bin/mpiexec`. You can append `PETSC_DIR` and `PETSC_ARCH` in a start up files (like in `~/.bashrc`) and make an alias to `mpiexec` or put it in your `$PATH` variable for it to be available in every shell you open.


## Compile your code mpf.c
  ```bash
     make PhaseFieldPet
  ```

## Run your simulation
### Static Triple Junction Simulation
- Boundary condition in x direction is pinned.
#### On CPU with n # mpi processes
Solve with default time step solver set in the code (Adaptive Runge Kutta Implicit-Explicit) and see as time progress (`ts_monitor`).
  ```bash
- mpiexec -n 4 ./PhaseFieldPet -ts_monitor
  ```
Restrict each $\phi$'s to be in Gibbs simplex:
  ```bash
- mpiexec -n 4 ./PhaseFieldPet -simplex
  ```
Use dot gradient formulation with obstacle potentials due to Steinbach 
 ```bash
 - mpiexec -n 4 ./PhaseFieldPet -grad_dot  -pot_steinbach  -simplex
  ```
Save output results approximately every 10 seconds (default is 100 seconds).
 ```bash
 - mpiexec -n 4 ./PhaseFieldPet -grad_dot  -pot_steinbach  -simplex -twrite 10
  ```
Change the underlying non linear  (Newton) solver  to one iteration
```bash
- mpiexec  -n 4 ./PhaseFieldPet  -snes_type ksponly
 ```
Use Matrix Free Non linear solver
```bash
- mpiexec  -n 4 ./PhaseFieldPet  -snes_mf -snes_type ksponly
 ```
Use Fully Implicit adaptive backward Differentiation Formula
```bash
- mpiexec  -n 4 ./PhaseFieldPet  -ts_type bdf 
 ```
Use interpolated formulation to grad and Multiorder-parameter simulation
```bash
- mpiexec -n 8 ./PhaseFieldPet -grad_interpolated -pfe_mop
 ```

Increase grid points to 256 x 256 x3
```bash
- mpiexec -n 80 PhaseFieldPet -simplex -ts_type bdf -da_grid_x 256 -da_grid_y 256
 ```
#### On a GPU (cuda based)
  ```bash
- mpiexec -n 1 PhaseFieldPet -simplex -ts_type bdf -da_grid_x 256 -da_grid_y 256 -dm_mat_type aijcusparse -dm_vec_type cuda
  ```
### Steady-state motion of triple junction 
- Boundary condition in x is set to be homogenous Neumann. All other options you have experimented above can also be added in steady-state case.
```bash
- mpiexec -n 4 ./PhaseFieldPet -bcx_neumann -snes_type ksponly  -ts_monitor
```
### For Developers*
If you have your own energy expressions (gradient , potential)  and bulk driving term, you can add it to PhaseFieldPet easily by including the respective `case` clause in the stiff `IRHSLocal()` and or non stiff terms in `RHSLocal()` functions in PhaseFieldPet.c. You can also change the `InitialMicrostructure()` function in PhaseFieldPet.c  to suit other simulations than example triple junction application described here or read initial phase field data available from other software or experiment.
#### Exercise
- Include a constant driving force term to the triple junction example above. See [Hoffroge et al., 2025](https://doi.org/10.1088/1361-651X/ad8d6f) for help or contact us.




