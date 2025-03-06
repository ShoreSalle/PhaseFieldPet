static char help[] = "Solves Phase Field Equation in 3D by Finite Difference Method using IMEX scheme \n" ;

/*===========================================================================
 The Phase Field Equation (pfe):

     dphi_alpha/dt = Gradient term + Potential term + Driving Force term = RHS

     See the associated paper to this code (Chota et al 2025,...) for mathematical expressions of RHS

 Reference
    1. Daubner et al 2023, https://doi.org/10.1016/j.commatsci.2022.111995.
       ( For Phase Field Equation, Gradient, Potential energy terms)
    2. P W Hoffrogge et al 2025, https://doi.org/10.1088/1361-651X/ad8d6f.
       ( For Driving Force term, section 4.2)

 By default
    Initilaze with  3 phases in [0,100] x [0,100] x [0,1]
      phi_0 = 1 if y > 80;       0 otherwise
      phi_1 = 1 if y < 80, x<50; 0 otherwise
      phi_2 = 1 if y < 80, x>50; 0 otherwise
    Boundary Conditions
      Dirichilet in X (pinned)
      Neuman BC = 0 in Y
      Periodic in Z
    Phase Field Equation(pfe)
      Lagrangian based Multiphase field equation (pfe_mpfl)
    Energy terms
      Grad term      = grad_dot
      Potential term = pot_toth
      Driving Force  = 0

 Usage (NB The following two lines are equivalent)
    - mpiexec -n 4 ./PhaseFieldPet
    - mpiexec -n 4 ./PhaseFieldPet -grad_dot -pot_toth -pfe_mpfl

 You can alter these default simulation at run time
    Model related options example
    - mpiexec -n 4 ./PhaseFieldPet -simplex
    - mpiexec -n 4 ./PhaseFieldPet -pfe_mpf -simplex
    - mpiexec -n 4 ./PhaseFieldPet -grad_weighted  -pot_nestler  -simplex
    - mpiexec -n 4 ./PhaseFieldPet -grad_dot  -pot_steinbach  -simplex
    - mpiexec -n 8 ./PhaseFieldPet -grad_interpolated -pfe_mop

    Solver related options example
    - mpiexec  -n 4 ./PhaseFieldPet  -snes_type ksponly
    - mpiexec  -n 4 ./PhaseFieldPet  -snes_mf -snes_type ksponly
    - mpiexec  -n 4 ./PhaseFieldPet  -ts_type bdf
    - mpiexec  -n 1 ./PhaseFieldPet  -dm_mat_type aijcusparse -dm_vec_type cuda
=============================================================================*/

#include "petsc.h"

#define np 3

typedef struct {
   PetscReal   phi[np];
}Field;

typedef struct {
   PetscInt    nx, ny, nz;
   PetscReal   x_ref, O_ref, Lx, Ly,Lz, dx, dy, dz, eps, mab, K;

   PetscReal   k[np][np], O[np][np],O_abc[np][np][np], M0, dx_inv, dy_inv,dz_inv, dx_sq_inv, dy_sq_inv,dz_sq_inv;
   PetscInt    twrite, pot, grad, pfe, bulk;
   PetscBool   pfe_mpf, pfe_mpfl, pfe_mop, grad_dot, grad_weighted, grad_interpolated, pot_toth,
               pot_moelans, pot_steinbach,pot_nestler, pot_garacke, bcx_neumann, bcy_dirichilet, simplex,
               bulk_b0, bulk_b1, bulk_b2;
}Param;

extern void InitParameters(Param*);
extern void InitParamRHS(PetscReal, Param*);
extern PetscErrorCode InitialMicroStructure(DM, Vec, Param*);
extern PetscErrorCode RHSLocal(DMDALocalInfo*, PetscReal, Field***, Field***, Param*);
extern PetscErrorCode IRHSLocal(DMDALocalInfo*, PetscReal ,Field*** , Field*** , Field*** , Param*);
extern PetscErrorCode ApplySimplex(TS );
extern PetscErrorCode WriteOutput(TS, PetscInt, PetscReal, Vec, void*);

   /*=========================================================================================================*/
   /*                                Main Function                                                            */
   /*=========================================================================================================*/

int main(int argc,char **argv) {
   Param          user;
   TS             ts;
   Vec            x;
   DM             da;
   DMDALocalInfo  info;
   PetscReal      dx, dy, dz, eps_fact;


   InitParameters(&user);

   PetscInitialize(&argc, &argv, NULL, help);

     DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DM_BOUNDARY_PERIODIC, DMDA_STENCIL_STAR,
               user.nx, user.ny, user.nz, PETSC_DECIDE, PETSC_DECIDE, 1, np, 1, NULL,NULL,NULL, &da);
     DMSetFromOptions(da);
     DMSetUp(da);
     DMDASetFieldName(da,0,"phi_1");
     DMDASetFieldName(da,1,"phi_2");
     DMDASetFieldName(da,2,"phi_3");
     DMDASetUniformCoordinates(da, 0, user.Lx, 0, user.Ly, 0.0, user.Lz);

     DMDAGetLocalInfo(da,&info);
     dx             = user.Lx /(info.mx-1);
     dy             = user.Ly /(info.my-1);
     dz             = user.Lz / info.mz   ;
     user.dx_inv    = 1.0/dx;
     user.dy_inv    = 1.0/dy;
     user.dz_inv    = 1.0/dz;
     user.dx_sq_inv = 1.0/(dx*dx);
     user.dy_sq_inv = 1.0/(dy*dy);
     user.dz_sq_inv = 1.0/(dz*dz);
     eps_fact       = 5.0 ;   /* 4.8634 */

     PetscOptionsGetReal(NULL,NULL,"-eps_fact",    &eps_fact   , NULL);
     user.eps         = eps_fact*dx;
     user.twrite      = 100;
     PetscOptionsGetInt (NULL,NULL,"-twrite", &user.twrite, NULL);

     user.pot_toth = PETSC_TRUE;
     user.pot = 0;
     PetscOptionsGetBool(NULL,NULL,"-pot_toth"     ,  &user.pot_toth      , NULL);
     user.pot_moelans = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-pot_moelans"  ,  &user.pot_moelans   , NULL);
     if (user.pot_moelans) { user.pot = 1;}
     user.pot_garacke = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-pot_garacke"  ,  &user.pot_garacke   , NULL);
     if (user.pot_garacke) { user.pot = 2;}
     user.pot_steinbach = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-pot_steinbach",  &user.pot_steinbach , NULL);
     if (user.pot_steinbach){ user.K = 16.0/(PETSC_PI*PETSC_PI); user.pot = 3;}
     user.pot_nestler = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-pot_nestler",  &user.pot_nestler     , NULL);
     if (user.pot_nestler){ user.K = 16.0/(PETSC_PI*PETSC_PI);user.pot = 4;}

     InitParamRHS(user.eps, &user);

     user.grad_dot = PETSC_TRUE;
     user.grad = 0;
     PetscOptionsGetBool(NULL,NULL,"-grad_dot"         ,  &user.grad_dot         , NULL);
     user.grad_weighted = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-grad_weighted"    ,  &user.grad_weighted    , NULL);
     if (user.grad_weighted)  { user.grad = 1;}
     user.grad_interpolated = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-grad_interpolated" , &user.grad_interpolated, NULL);
     if (user.grad_interpolated)  { user.grad = 2;}

     user.M0 = user.mab/user.eps;
     user.pfe_mpfl=PETSC_TRUE;
     user.pfe = 0;
     PetscOptionsGetBool(NULL,NULL,"-pfe_mpfl"     ,  &user.pfe_mpfl      , NULL);
     user.pfe_mpf = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-pfe_mpf"      ,  &user.pfe_mpf       , NULL);
     if ( user.pfe_mpf) { user.pfe = 1; user.M0 = user.mab/(np*user.eps);}
     user.pfe_mop =PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-pfe_mop"      ,  &user.pfe_mop       , NULL);
     if ( user.pfe_mop) { user.pfe = 2;}

     user.bulk = 0;
     user.bulk_b0 = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-bulk_b0"     ,  &user.bulk_b0        , NULL);
     if ( user.bulk_b0) { user.bulk = 1;}
     user.bulk_b1 = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-bulk_b1"     ,  &user.bulk_b1        , NULL);
     if ( user.bulk_b1) { user.bulk = 2;}
     user.bulk_b2 = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-bulk_b2"     ,  &user.bulk_b2        , NULL);
     if ( user.bulk_b2) { user.bulk = 3;}

     user.bcx_neumann = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-bcx_neumann"   ,    &user.bcx_neumann, NULL);
     user.bcy_dirichilet = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-bcy_dirichilet", &user.bcy_dirichilet, NULL);

     TSCreate(PETSC_COMM_WORLD,&ts);
     TSSetProblemType(ts,TS_NONLINEAR);
     TSSetDM(ts,da);
     TSSetApplicationContext(ts,&user);

     DMDATSSetRHSFunctionLocal(da,INSERT_VALUES, (DMDATSRHSFunctionLocal)RHSLocal ,&user);
     DMDATSSetIFunctionLocal  (da,INSERT_VALUES, (DMDATSIFunctionLocal)IRHSLocal  ,&user);

     TSSetType(ts,TSARKIMEX);
     TSSetTime(ts,0.0);              //to
     TSSetMaxTime(ts,3500.0);        //tf
     TSSetMaxSteps(ts,1e9);          //iter
     TSSetTimeStep(ts,0.1);          //dt
     TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);

     DMCreateGlobalVector(da,&x);

     user.simplex = PETSC_FALSE;
     PetscOptionsGetBool(NULL,NULL,"-simplex", &user.simplex, NULL);
     if (user.simplex){
        TSSetPostStep(ts, ApplySimplex);
        TSRestartStep(ts);
     }

     TSMonitorSet(ts, WriteOutput, &user, NULL);

     TSSetFromOptions(ts);

     InitialMicroStructure(da, x, &user);

     TSSolve(ts,x);

     VecDestroy(&x);
     TSDestroy(&ts);
     DMDestroy(&da);
   PetscFinalize();
   return 0;
}

  /* ==============================================================================*/
  /*                              Simulation Parameters                            */
  /* ==============================================================================*/

void InitParameters(Param* user) {

    PetscReal  w, h,l, M_ab, t_ref;

    user->nx    = 128;
    user->ny    = 128;
    user->nz    = 3  ;
    user->K     = 9.0;

     /* In SI Unit */
     w     = 100e-6;
     h     = 100e-6;  // [100,400]e-6
     l     = 1e-6  ;
     M_ab  = 1e-14 ;
    /* gamma[l][m]   see below*/

   /* Non dimensionalizing factors [s],[m] and [J/m3] */
     t_ref       = 100;
     user->x_ref = 1e-6;
     user->O_ref = 1e6;

    /* Non dimensionalized quantities */
     user->Lx    = w/user->x_ref;
     user->Ly    = h/user->x_ref;
     user->Lz    = l/user->x_ref;
     user->mab   = M_ab* user->O_ref*t_ref/user->x_ref;
}

void InitParamRHS(PetscReal eps, Param* user) {

    PetscReal gamma[np][np], g[np][np];

    for (int l=0; l< np; l++){
       for (int m=0; m< np; m++){
          if (m==l){
              continue;
          }
          gamma[l][0] = 1.0;                         /* in SI , gamma_a0 */
          if (m != 0){
              gamma[l][m] = 1.0 * gamma[l][0];        /* gamma_ab = [0.1, 2.0] * gamma_a0 */
          }
          g[l][m]     = gamma[l][m]/(user->O_ref*user->x_ref);

          /* RHS parameters  */
          user->k[l][m]   =  eps*g[l][m];
          user->O[l][m]   =  user->K*g[l][m]/eps;
          for (int n = m+1; n< np; n++){
              user->O_abc[l][m][n] =  0.01*user->O[l][m];
          }
       }
    }
}

    /* =====================================================================================*/
    /*                      Initial Microstructure (Phi)                                    */
    /* =====================================================================================*/

PetscErrorCode  InitialMicroStructure(DM da, Vec Y, Param* user) {

   DMDALocalInfo    info;
   PetscInt         i,j,k;
   DMDACoor3d       ***aC;
   Field            ***aY;

   PetscFunctionBeginUser;
   VecSet(Y,0.0);
   if(user->pfe_mop){ VecSet(Y,0.000001);}

   DMDAGetLocalInfo(da,&info);
   DMDAGetCoordinateArray(da,&aC);
   DMDAVecGetArray(da,Y,&aY); //  DMDAVecGetArrayWrite(da,Y,&aY);
   for (k = info.zs; k < info.zs + info.zm; k++) {
     for (j = info.ys; j < info.ys+info.ym; j++) {
        for (i = info.xs; i < info.xs+info.xm; i++) {
           if (aC[k][j][i].y <0.8*user->Ly) {
              if (aC[k][j][i].x < 0.5*user->Lx){
                 aY[k][j][i].phi[1] = 1.0;
              }else {
                 aY[k][j][i].phi[2] = 1.0;
              }
           }else {
              aY[k][j][i].phi[0]    = 1.0;
           }
        }
     }
   }
   DMDAVecRestoreArray(da,Y,&aY);
   DMDARestoreCoordinateArray(da,&aC);

   PetscFunctionReturn(PETSC_SUCCESS);
}

    /* ========================================================================================================*/
    /*                            Gradient Potential                                                           */
    /* ========================================================================================================*/

  PetscErrorCode IRHSLocal(DMDALocalInfo *info, PetscReal t,Field ***aY, Field ***aYdot, Field ***aF, Param *user) {

  PetscInt       i, j,k, l,m, mx = info->mx, my = info->my;
  PetscReal      ul,ur,ut,ub,ufr, uba, dFgrad_dphi[np], PhibSq, sumSquares,khat,laplace[np], grad_phix[np],grad_phiy[np],
                 grad_phiz[np],dkhat_dphi, rhs1, GradSquare, sumGradSquare, grad_kx, grad_ky,grad_kz, PhiaSq;

  PetscFunctionBeginUser;
  for (k = info->zs; k < info->zs + info->zm; k++) {
     for (j = info->ys; j < info->ys + info->ym; j++) {
        for (i = info->xs  ; i <  info->xs + info->xm; i++) {

           sumGradSquare = 0.0;
           for (l = 0; l < np; l++) {
               if (user->bcx_neumann){
                  ul =(j == 0)    ? aY[k][j][i+1].phi[l]  : aY[k][j][i-1].phi[l] ;
                  ur =(j == my-1) ? aY[k][j][i-1].phi[l]  : aY[k][j][i+1].phi[l] ;
               }else {  /*bcx_dirichilet */
                  ul =(i == 0)    ? 2.0*aY[k][j][i].phi[l] - aY[k][j][i+1].phi[l] : aY[k][j][i-1].phi[l];
                  ur =(i == mx-1) ? 2.0*aY[k][j][i].phi[l] - aY[k][j][i-1].phi[l] : aY[k][j][i+1].phi[l];
               }

               if(user->bcy_dirichilet){
                  ub =(j == 0)    ? 2.0*aY[k][j][i].phi[l] - aY[k][j+1][i].phi[l] : aY[k][j-1][i].phi[l];
                  ut =(j == my-1) ? 2.0*aY[k][j][i].phi[l] - aY[k][j-1][i].phi[l] : aY[k][j+1][i].phi[l];
               }
               else{   /*user->bcy_neumann*/
                  ub =(j == 0)    ? aY[k][j+1][i].phi[l]  : aY[k][j-1][i].phi[l] ;
                  ut =(j == my-1) ? aY[k][j-1][i].phi[l]  : aY[k][j+1][i].phi[l] ;
               }

               ufr = aY[k+1][j][i].phi[l] ; /* z, periodic */
               uba = aY[k-1][j][i].phi[l] ;

               laplace[l] = (ur -2.0*aY[k][j][i].phi[l] + ul) *user->dx_sq_inv + (ut -2.0*aY[k][j][i].phi[l] + ub)*
                             user->dy_sq_inv + (ufr -2.0*aY[k][j][i].phi[l] + uba)*user->dz_sq_inv;

               if (user->grad_interpolated || user->grad_weighted){
                   grad_phix[l] = 0.5*(ur  - ul )*user->dx_inv ;
                   grad_phiy[l] = 0.5*(ut  - ub )*user->dy_inv ;
                   grad_phiz[l] = 0.5*(ufr - uba)*user->dz_inv ;

                   sumGradSquare += grad_phix[l]*grad_phix[l] + grad_phiy[l]*grad_phiy[l] + grad_phiz[l]*grad_phiz[l];
               }

           }

           rhs1       = 0.0;
           khat       = 0.0;
           sumSquares = 0.0;
           if (  user->grad_interpolated ){
             //khat =0.0; sumSquares = 0.0;
              for (int l = 0; l < np; l++) {
                  PhiaSq = aY[k][j][i].phi[l]*aY[k][j][i].phi[l];
                  for (int m = l + 1; m < np; m++) {
                      PhibSq      =  aY[k][j][i].phi[m]*aY[k][j][i].phi[m];
                      khat       +=  user->k[l][m] * PhiaSq * PhibSq;
                      sumSquares +=  PhiaSq *PhibSq ;
                  }
              }
              if (sumSquares <= 1e-15) { sumSquares = 1.0;}
              khat = khat/sumSquares;
           }

           switch (user->grad) {
              case 1:   /* Gradient (Weighted | generalized ) - Nestler  */
                   for (l = 0; l < np; l++) { /* alpha */
                        dFgrad_dphi[l] = 0.0;
                        for (m = 0; m < np; m++) { /* betta */
                            if (m==l) {
                               continue;
                            }
                            GradSquare = grad_phix[m]*grad_phix[m] + grad_phiy[m]*grad_phiy[m] + grad_phiz[m]*grad_phiz[m] ;
                            dFgrad_dphi[l] += 2.0*user->k[l][m]*(
                                              2.0*aY[k][j][i].phi[l] * ( GradSquare)
                                            - 2.0*aY[k][j][i].phi[m] * ( grad_phix[l]* grad_phix[m] + grad_phiy[l]* grad_phiy[m]+
                                                                         grad_phiz[l]* grad_phiz[m]  )
                                            + aY[k][j][i].phi[l]*aY[k][j][i].phi[m] * laplace[m]
                                            - aY[k][j][i].phi[m]*aY[k][j][i].phi[m] * laplace[l]  );
                        }
                        rhs1  += dFgrad_dphi[l];
                   }
                   break;
              case 2:  /* Gradient (Interpolated) - Moelans */
                   grad_kx = 0.0;
                   grad_ky = 0.0;
                   grad_kz = 0.0;
                   for ( l = 0; l < np; l++) { // alpha
                        dFgrad_dphi[l] = 0.0;
                        dkhat_dphi     = 0.0;
                        for (m = 0; m < np; m++) { // betta
                            if (m==l) {
                               continue;
                            }
                            PhibSq       =  aY[k][j][i].phi[m]*aY[k][j][i].phi[m];
                            dkhat_dphi  +=  2.0 * aY[k][j][i].phi[l] * (user->k[l][m] - khat)* PhibSq ;
                        }
                        dkhat_dphi     =  dkhat_dphi /sumSquares;

                        grad_kx +=  dkhat_dphi * grad_phix[l];
                        grad_ky +=  dkhat_dphi * grad_phiy[l];
                        grad_kz +=  dkhat_dphi * grad_phiz[l];

                        dFgrad_dphi[l] = 0.5*dkhat_dphi *sumGradSquare - khat*laplace [l];
                   }

                   for (l = 0; l < np; l++) { // alpha
                        dFgrad_dphi[l] += - grad_kx*grad_phix[l] - grad_ky*grad_phiy[l] - grad_kz*grad_phiz[l];
                        rhs1  += dFgrad_dphi[l];
                   }
                   break;
              default:  /* Gradient (Dot) */
                   for (l = 0; l < np; l++) { // alpha
                        dFgrad_dphi[l] = 0.0;
                        for (m = 0; m < np; m++) { // betta
                             if (m==l) {
                                continue;
                             }
                             dFgrad_dphi[l] +=  user->k[l][m] * laplace[m];
                        }
                        rhs1 += dFgrad_dphi[l];
                   }
                   break;
           }

           for (l = 0; l < np; l++) { // alpha
              switch (user->pfe){
                   case 1:
                      aF[k][j][i].phi[l]  = aYdot[k][j][i].phi[l] + user->M0 *( np*dFgrad_dphi[l] -   rhs1   ); // MPF
                      break;
                   case 2:
                      aF[k][j][i].phi[l]  = aYdot[k][j][i].phi[l] + user->M0 *(    dFgrad_dphi[l]           );  // MOP
                      break;
                   default:
                      aF[k][j][i].phi[l]  = aYdot[k][j][i].phi[l] + user->M0 *(    dFgrad_dphi[l] -(rhs1/np)) ; // MPFL
                      break;
              }
           }

        }
     }
  }
   PetscFunctionReturn(PETSC_SUCCESS);
}

   /*================================================================================================*/
   /*                           Well / Obstacle Potentials , Driving Force                           */
   /*================================================================================================*/

PetscErrorCode RHSLocal(DMDALocalInfo *info, PetscReal t, Field ***aY, Field ***aG, Param *user) {

  PetscInt       i, j,k,l,m,n;
  PetscReal      PhiaSq, last, dFpot_dphi[np], PhibSq, sumSquares, Ohat,first, dOhat_dphi,rhs2,df_bulk, bulk_force[np];

  PetscFunctionBeginUser;
  for (k = info->zs; k < info->zs + info->zm; k++) {
     for (j = info->ys; j < info->ys + info->ym; j++) {
        for (i = info->xs  ; i <  info->xs + info->xm; i++) {

           Ohat       = 0.0;
           sumSquares = 0.0;
           if ( user->pot_moelans || user->pot_toth ){
               //Ohat =0.0; khat =0.0; sumSquares = 0.0;
               for (int l = 0; l < np; l++) {
                   PhiaSq = aY[k][j][i].phi[l]*aY[k][j][i].phi[l];
                   for (int m = l + 1; m < np; m++) {
                      PhibSq      =  aY[k][j][i].phi[m]*aY[k][j][i].phi[m];
                      Ohat       +=  user->O[l][m] * PhiaSq * PhibSq;
                      sumSquares +=  PhiaSq *PhibSq ;
                   }
               }
               if (sumSquares <= 1e-15) {sumSquares = 1.0;}
               Ohat = Ohat/sumSquares;
           }

           rhs2 =0.0;
           switch (user->pot) {
              case 1: /* Well Pot - Moelans */
                   for (l = 0; l < np; l++) { // alpha
                        first      = 0.0;
                        last       = 0.0;
                        dOhat_dphi = 0.0;
                        for (m = 0; m < np; m++) { // betta
                            if (m==l) {
                               continue;
                            }
                            PhibSq       =  aY[k][j][i].phi[m]*aY[k][j][i].phi[m]; // if pot is separate , remove this.
                            first       +=  user->O[l][m] * PhibSq ;
                            dOhat_dphi  +=  2.0 * aY[k][j][i].phi[l] * (user->O[l][m] - Ohat)* PhibSq ;
                        }
                        dOhat_dphi    =  dOhat_dphi /sumSquares;
                        last         +=  aY[k][j][i].phi[l]*aY[k][j][i].phi[l]*(0.25*aY[k][j][i].phi[l]* aY[k][j][i].phi[l] - 0.5);
                        dFpot_dphi[l] =  aY[k][j][i].phi[l]*1.5*first  +  0.5*Ohat*aY[k][j][i].phi[l]*( aY[k][j][i].phi[l]*aY[k][j][i].phi[l] -1.0) +
                                         0.5*dOhat_dphi*(0.25 + last);
                        rhs2 += dFpot_dphi[l];
                   }
                   break;
              case 2: /* Well  Pot - Garacke */
                   for (l = 0; l < np; l++) { // alpha
                        dFpot_dphi[l] = 0.0;
                        for (m = 0; m < np; m++) { // betta
                            if (m==l) {
                               continue;
                            }
                            PhibSq = aY[k][j][i].phi[m]*aY[k][j][i].phi[m] ;
                            dFpot_dphi[l]  +=  2.0*aY[k][j][i].phi[l] *user->O[l][m]*PhibSq ;
                            for (n = m+1; n < np; n++) { // gamma
                                    dFpot_dphi[l]  +=  2.0*aY[k][j][i].phi[l] * user->O_abc[l][m][n]*PhibSq* aY[k][j][i].phi[n]*aY[k][j][i].phi[n];
                            }
                        }
                        rhs2 += dFpot_dphi[l];
                   }
                   break;
              case 3: /* Obstacle Pot1 , pot_steinbach  */
                   for (l = 0; l < np; l++) { // alpha
                        dFpot_dphi[l] = 0.0;
                        for (m = 0; m < np; m++) { // betta
                            if (m==l) {
                               continue;
                            }
                            dFpot_dphi[l]  +=  user->O[l][m]*aY[k][j][i].phi[m] ;
                        }
                        rhs2 += dFpot_dphi[l];
                   }
                   break;

              case 4: /* Obstacle Pot2 , pot_nestler */
                   for (l = 0; l < np; l++) { // alpha
                        dFpot_dphi[l] = 0.0;
                        for (m = 0; m < np; m++) { // betta
                            if (m==l) {
                               continue;
                            }
                            dFpot_dphi[l]  +=  user->O[l][m]*aY[k][j][i].phi[m] ;

                            for (n = m+1; n < np; n++) { // gamma
                                    dFpot_dphi[l]  +=  user->O_abc[l][m][n]*aY[k][j][i].phi[m]*aY[k][j][i].phi[n];
                            }
                        }
                        rhs2 += dFpot_dphi[l];
                   }
                   break;
              default: /* Well Pot - Toth */
                   for (l = 0; l < np; l++) { // alpha
                        first      = 0.0;
                        last       = 0.0;
                        dOhat_dphi = 0.0;
                        for (m = 0; m < np; m++) { // betta
                            if (m==l) {
                               continue;
                            }
                            PhibSq       =  aY[k][j][i].phi[m]*aY[k][j][i].phi[m]; // if pot is separate , remove this.
                            first       +=  user->O[l][m] * PhibSq ;
                            dOhat_dphi  +=  2.0 * aY[k][j][i].phi[l] * (user->O[l][m] - Ohat)* PhibSq ;
                        }
                        dOhat_dphi    =  dOhat_dphi /sumSquares;
                        last         +=  aY[k][j][i].phi[l]*aY[k][j][i].phi[l]*aY[k][j][i].phi[l]*( 0.25* aY[k][j][i].phi[l] - 0.33333333)    ;
                        dFpot_dphi[l] =  aY[k][j][i].phi[l]*first   +  Ohat*aY[k][j][i].phi[l]* aY[k][j][i].phi[l]*(aY[k][j][i].phi[l] -1.0)  +
                                         dOhat_dphi*(0.083333333 + last);
                        rhs2 += dFpot_dphi[l];
                   }
                   break;
           }

           /* Bulk Driving Force*/
           for (l = 0; l < np; l++) { // alpha
              bulk_force[l] = 0.0;
              if (user->bulk) {
                   for (m = 0; m < np; m++) { // betta
                        if (m==l){
                           continue;
                        }
                        /*sqrt(phia*phib) ~= 0.5*(1.0 + phia*phib) */
                        switch (user->bulk){
                           case 1:
                               bulk_force[l] +=  1.0 + aY[k][j][i].phi[l]*aY[k][j][i].phi[m];                                            // B0
                               break;
                           case 2:
                               bulk_force[l] +=  (2.0/np)*(1.0 + aY[k][j][i].phi[l]*aY[k][j][i].phi[m]);                                 // B1
                               break;
                           case 3:
                               bulk_force[l] += (aY[k][j][i].phi[l] + aY[k][j][i].phi[m])*(1.0 + aY[k][j][i].phi[l]*aY[k][j][i].phi[m]); // B2
                               break;
                           default:
                               bulk_force[l]  = 0.0;
                               break;
                       }
                   }
                   df_bulk = 2.0/user->eps;  // choose df_bulk: phat_pf = etta*df_bulk/gamma0 = [-2,2]
                   bulk_force[l] *= 0.5*PETSC_PI*df_bulk;
              }
           }


           for (l = 0; l < np; l++) { // alpha
              switch (user->pfe){
                   case 1:
                   //aG[k][j][i].phi[l]  = - user->M0 *( np*dFpot_dphi[l]  -    rhs2 );  // MPF
                     aG[k][j][i].phi[l]  = - user->M0 *( np*dFpot_dphi[l]  -    rhs2  - np*bulk_force[l]);  // MPF

                       break;
                   case 2:
                       aG[k][j][i].phi[l]  = - user->M0 *(    dFpot_dphi[l]            );  // MOP
                       break;
                   default:
                       aG[k][j][i].phi[l]  = - user->M0 *(    dFpot_dphi[l] - (rhs2/np));  // MPFL
                       break;
              }
           }


        }
     }
  }
   PetscFunctionReturn(PETSC_SUCCESS);
}

   /*==============================================================================================*/
   /*                     Write out results in vtk                                                 */
   /*==============================================================================================*/

PetscErrorCode WriteOutput(TS ts, PetscInt step, PetscReal crtime, Vec x, void *ctx){

    Param *user = (Param*)ctx;
    char filename[256];
    static PetscReal  next_write_time = 0.0;
    PetscViewer viewer;

    PetscFunctionBeginUser;
    if (crtime >= next_write_time) {
        sprintf(filename, "phi_%.0f.vts", crtime);
        PetscViewerVTKOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);
        VecView(x,viewer);
        PetscViewerDestroy(&viewer);
        next_write_time = crtime + user->twrite;
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

   /*===============================================================================================*/
   /*                            Limit Phase Fields to Gibbs Simplex                                */
   /*===============================================================================================*/

PetscErrorCode ApplySimplex(TS ts){

    PetscInt        i, j,k, l ;
    Vec             Y ;
    DM              da;
    DMDALocalInfo   info;
    Field           ***u;
    PetscReal       sumphi;

    PetscFunctionBeginUser;
    TSGetSolution(ts, &Y);
    TSGetDM(ts, &da);
    DMDAGetLocalInfo(da, &info);
    DMDAVecGetArray(da, Y, &u);
    for ( k = info.zs; k < info.zs + info.zm; k++) {
        for ( j = info.ys; j < info.ys + info.ym; j++) {
           for ( i = info.xs; i < info.xs + info.xm; i++) {
              sumphi = 0.0;
              for (l =0; l<np; l++){
                 if(u[k][j][i].phi[l] < 0){ u[k][j][i].phi[l] = 0.0;}
                 sumphi += u[k][j][i].phi[l];
              }
              for (l =0; l<np; l++){
                 u[k][j][i].phi[l]  /= sumphi;
              }
           }
        }
    }
    DMDAVecRestoreArray(da, Y, &u);
    PetscFunctionReturn(PETSC_SUCCESS);
}

