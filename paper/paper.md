---
title: 'PhaseFieldPet: An Open-Source Phase Field Modeling Software for Heterogeneous Architectures'
tags:
  - Phase Field
  - High Performance Computing
  - Message Passing Interface
authors:
  - name: Shore Salle Chota
    corresponding: true
    orcid: 0009-0003-0529-3594
    affiliation: "1,4"
  - name: Martin Reder
    orcid: 0000-0002-7503-9351
    affiliation: "1, 2"
  - name: Daniel Schneider
    orcid: 0000-0002-9250-2918
    affiliation: "1, 2"
  - name: Harald Koestler
    orcid: 0000-0002-6992-2690
    affiliation: 4
  - name: Britta Nestler
    orcid: 0000-0002-3768-3277
    affiliation: "1, 2, 3"

affiliations:
 - name: Institute for Applied Materials (IAM), Karlsruhe Institute of Technology, Karlsruhe, 76131, Germany
   index: 1
 - name: Institute of Digital Materials Science (IDM), Karlsruhe University of Applied Sciences, Karlsruhe, 76133, Germany
   index: 2
 - name: Institute of Nanotechnology (INT), Karlsruhe Institute of Technology, Karlsruhe, 76131, Germany
   index: 3
 - name: Chair for System Simulation, Friedrich-Alexander-Universität Erlangen-Nürnberg, Erlangen, Germany
   index: 4

date: 17 February 2025
bibliography: paper.bib
---

# Summary

Phase field method has emerged as a powerful computational tool for simulating  complex, evolving interfaces at mesoscale  in materials science and Engineering , fluid dynamics, cell migration and other fields. However, achieving scalability and efficiency for multicomponent multiphase across diverse hardware architectures remains a challenge.  This paper presents PhaseFieldPet, an open-source, message passing interface (mpi) based software package for large-scale phase field simulations, specifically designed to leverage heterogeneous architectures such as CPUs  and GPUs. PhaseFieldPet is built upon the Portable, Extensible Toolkit for Scientific Computation (PETSc), to efficiently handle large-scale simulations on high-performance computing platforms. The software's modular design facilitates and easily integrates various phase field models such as Multiphase-Field and Multi-Order Parameter models with various choice of gradient and potential energy contributions. Performance benchmarks demonstrate the software’s capability to handle simulations involving from small to millions of degrees of freedom with excellent scalability.

# The Phase Field Equation

The phase-field approach  represents surfaces and interfaces implicitly using continuous scalar fields (order parameter) $\phi_{\alpha}(\mathbf{r},t)$, $\alpha  \in \{1, 2, \ldots, N\}$, which is constant in bulk phases and transition smoothly—but sharply—across a diffuse boundary. $\phi_{\alpha}(\mathbf{r},t)$  represents different grains or  phases such as solid, liquid, gas and their state like  crystallographic orientation, polarization, volume fraction. For example, in a solid-solid phase transition, a solid phase can be represented by N order parameters $\phi_{\alpha}(\mathbf{r},t)$ based on N  crystallographic orientations or grains.

Microstructure evolution (and hence evolution of order parameter $\phi_{\alpha}(\mathbf{r},t)$) can be obtained  from functionals of entropy, or free energy, or grand potential [@Hötzer:2018]. Following  energy functional, one can write total energy functional of the system as
$$\mathcal{F}(\phi,\nabla\phi) =\int_{V} f dV=  \int_{V} f_{grad}(\phi,\nabla\phi) + f_{pot}(\phi) + f_{bulk}(\phi,...)\, dV$$


The first two terms contributes to  interfacial energy , and exists during  phases are in contact with each other and is responsible to keep the diffuse interface at finite width with interplay of gradient energy density $f_{grad}(\phi,\nabla\phi)$  which diffuses the interfaces, while the potential term $f_{pot}(\phi)$counter acts it. The bulk  contribution $f_{bulk}(\phi,...)$ can be of various type depending on the problem at hand  such as chemical, thermal, mechanical, electrical, elastic, magnetic and etc.

There exists various formulations to $f_{grad}(\phi,\nabla\phi)$  and $f_{pot}(\phi)$  by many scholars in the Phase field community [@Daubner:2023]. @Nestler:2005 formulated these terms in their PACE3D software [@Hötzer:2018] by

$$f_{grad}(\phi,\nabla\phi)=\varepsilon \sum_{\alpha} \sum_{\beta > \alpha} \gamma_{\alpha \beta} |\phi_{\alpha} \nabla \phi_{\beta} - \phi_{\beta} \nabla  \phi_{\alpha}|^2 $$
$$f_{pot}(\phi)=\frac{16}{\varepsilon {\pi}^2}\sum_{\alpha} \sum_{\beta > \alpha} \gamma_{\alpha \beta} \phi_{\alpha} \phi_{\beta} + \frac{1}{\varepsilon}\sum_{\alpha} \sum_{\beta > \alpha} \sum_{\delta > \beta} \gamma_{\alpha \beta \delta} \phi_{\alpha} \phi_{\beta} \phi_{\gamma} $$

See @Daubner:2023 table B.4 for other formulations of these terms  and their explanations currently considered in PhaseFieldPet.

Phase field evolution equation is in general given by Allen-Cahn approach for each order parameter following total energy minimization principle of the system by

  $$\frac{\partial \phi_{\alpha}}{\partial t} = -L\frac{\delta\mathcal{F}}{\delta{\phi}_\alpha} = -L \left(\frac{\partial f}{\partial\phi_{\alpha}} - \nabla \cdot \frac{\partial f}{\partial\nabla\phi_{\alpha}}\right)$$

Where L is kinetic coefficient. We differentiate Phase field equation based on whether order parameters are conserved ( where each order parameter represents volume fraction and hence the sum of all  phase fields  at a point in space is 1), which are called  Multi Phase Field  (MPF) models and and those that do not put such a conservation requirement directly, called Multi order Parameter model (MOP). All these models and variants are treated and solved  in PhaseFieldPet, see the usage section below.

# Statement of need

For the past couple of decades, phase field software has been being developed and used with  in house codes, and Open source phase field software started to be available from 2007 [@Hong:2020]. Many  existing open source softwares are limited to one or two spatial dimensions, focus on binary systems, use  only one  type time step solver (usually explicit time stepping), works only on one CPU core (serial code) or are not capable of using heterogeneous compute resources available such as GPUs  for compute and energy efficiency. Notable large scale, distributed computing capable open source phase field softwares include: The open source Multiphysics Object Oriented Simulation Environment (MOOSE) [@schwen2023phasefield:2017] - which  is a powerful toolset for implementing phase field models using the finite element method and targets CPUs, PRISMS-PF - massively parallel finite element code for conducting phase field [@DeWitt:2020],   OpenPhase [@Tegeler:2017] uses  finite difference for spatial discretization and uses explicit time stepping algorithm and targets CPUs. Among non open source software, distributed machines capable software is a  Parallel Algorithms for Crystal Evolution in 3D (PACE3D)  [@Hötzer:2018] is a multiphase field software that uses explicit time stepping along with finite difference spatial discretizations.


PhaseFieldPet is a Finite difference based software built on top of  TS solver from PETSc [@Abhyankar:2018] and it extends our previous work [@Daubner:2023] to 3D and N phases, and fills the above mentioned gaps in existing softwares. Among others:

1. It paves the way to opensource PACE3D.

2. Allows multiphase simulation in 1D / 2D / 3D.

3. Decouple the numerical solution methods from the physical modeling such that one can choose various solution methods without restricting to one time        step solver (i.e. one can use methods like semi implicit, implicit time stepping algorithms, various underlying nonlinear solver ,  linear solvers and preconditioners, etc) based on composability features of PETSc [@petsc-web-page:2024]. This makes it easier for newcomers to the phase field community and advanced users alike. We formulate the phase field equation using the Implicit Explicit (IMEX) scheme such that the either the whole Phase field equation is treated implicitly or the stiff part, allowing longer time steps compared to the explicit methods which require small time steps for stability.

4. Works on single core, multicore to  multi node High Performance Computing cluster/supercomputer coupled with accelerators such as GPUs [@Mills:2021]. GPUs are  increasingly available for computing purposes with thousands of compute cores and small energy usage. We  can  leverage their  compute  power using  PhaseFieldPet from one code base for both CPUs and GPUs.

# Usage

To use PhaseFieldPet, all you need is  to install  PETSc [@petsc-web-page:2024], compile it using\
     &nbsp;&nbsp;&nbsp; `make PhaseFieldPet`\
run the executable generated with default solver settings  (E.g. with 4 mpi process) by\
     &nbsp;&nbsp;&nbsp; `mpiexec -n 4 PhaseFieldPet`\
By default PhaseFieldPet simulates a static triple junction, see @Daubner:2023, using dot gradient term  (`grad_dot`) and well potential of Toth (`pot_toth`) using lagrangian based phase field equations (`pfe_mpfl`). The default time stepping solver being Adaptive Runge Kutta  Implicit-Explicit (IMEX) method, where the stiff part of the equation is  treated implicitly.  `Pfe_mpfl` adds a lagrangian term to the phase field equation to constrain the sum of $\phi_{\alpha}$ to be 1 at every point in the domain[@Nestler:2005].

To simulate with non default combinations, for instance to solve phase field equation with gradient and potential energy terms above with $f_{bulk}(\phi,...)= 0$, we give the corresponding options to the run time as\
     &nbsp;&nbsp;&nbsp;`mpiexec -n 4 PhaseFieldPet -grad_weighted -pot_nestler  -simplex`\
This means that we are using weighted (generalized) gradient energy formulation (`grad_weighted`), the obstacle potential (`pot_nestler`) as outlined above and apply Gibbs simplex constraint (`simplex`) to constrain each phase field $\phi_{\alpha}>=0$, $\sum_{\alpha} \phi_{\alpha} =1$ at each point in  the simulation domain. One can also use different Phase field equations (pfe) with optionlike  `pfe_mop` (Multiorder parameter model) which does not put any restriction on order parameters, `pfe_mpf` (Non lagrangian based Multi-Phase field Equation) [@Tegeler:2017] along with other gradient and potential energy contributions.
For details of usage not mentioned here, including your own energy contributions, see the associated [Github page](https://github.com/ShoreSalle/PhaseFieldPet) to this paper.

## Example Performance Result

Here we report the strong scalabilty of PhaseFieldPet for simulation static triple junction using second order, Adaptive Backward Euler  (fully implicit) time step solver by increasing grid points along x and y direction using\
     &nbsp; `mpiexec -n # PhaseFieldPet -simplex -ts_type bdf -da_grid_x 256 -da_grid_y 256 `\
\autoref{fig:example} shows the result on Meggie cluster at  NHR@FAU obtained by running up  to 80 mpi process on 4 compute nodes,  where each nodes have two Intel Xeon E5-2630v4 “Broadwell” chips (10 cores per chip) running at 2.2 GHz with 25 MB Shared Cache per chip and 64 GB of RAM. The result indicated excellent agreement with theoretical expectations that the log-log plot is a straight decreasing line with slope of -1.


<div style="text-align:center;">
![Log-Log Plot of Execution Time vs MPI Processes.\label{fig:example}](timing.png){ width=70% height=70%}
</div>

# Conclusions

PhaseFieldPET provides users with flexible  methods to solve phase field equations, along with various energy contributions on heterogeneous  hardware architectures.  More specifically, a  user can include or choose what gradient energy  term, potential energy  term, the type of the phase field equation, type of numerical algorithm to use in order to  solve the differential equation (pfe) (along with the choice of the underlying nonlinear equation solver, Linear equation solver, preconditioners) and etc. Inline with  PACE3D software [@Hötzer:2018] and extension of it, the future version of PhaseFieldPet will  include various other modules corresponding to different applications of  phase field.

# Acknowledgement

We acknowledge discussions we had with Dr. Simon Daubner during the genesis of this project. The authors gratefully acknowledge the scientific support and HPC resources provided by the Erlangen National High Performance Computing Center (NHR@FAU) of the Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU). The hardware is funded by the German Research Foundation (DFG).

# Funding

The authors would like  to thank the [NHR-Verein e.V](www.nhr-verein.de) for supporting this work/project within the NHR Graduate School of National High Performance Computing (NHR).

# References

