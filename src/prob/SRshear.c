#include "copyright.h"
/*============================================================================*/
/*! \file jet.c
 *  \brief Sets up a jet introduced through L-x1 boundary (left edge) */
/*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2 -- random quantity from kh.c
 * grav_pot2 -- gravitational potential from rt.c
 *============================================================================*/

static double ran2(long int *idum);

/* Make radius of jet and jet variables static global so they can be accessed
   by BC functions */
static Real rjet,Bxjet;
static Prim1DS Wjet;
static Cons1DS Ujet;
static Real x1_mid,x2_mid,x3_mid;
static Real grav_pot2(const Real x1, const Real x2, const Real x3);

/*----------------------------------------------------------------------------*/
/* problem */

void problem(DomainS *pDomain){
   GridS *pGrid=(pDomain->Grid);
   int i, is = pGrid->is, ie = pGrid->ie;
   int j, js = pGrid->js, je = pGrid->je;
   int k, ks = pGrid->ks, ke = pGrid->ke;
   int il,iu,jl,ju,kl,ku;
   Prim1DS W;
   Cons1DS U;
   Real x1_min, x1_max;
   Real x2_min, x2_max;
   Real x3_min, x3_max;
   Real Bx=0.0;
   Real vxy = 0.0, vx, vy, x1,x2,x3, noise, d0, vxamp, g;
   Bxjet = 0.0;
  long int iseed = -1;

/* read parameters from input file */

   d0  = par_getd("problem", "d"); // ambient density
   U.d = d0;
   W.P  = par_getd("problem", "p"); // ambient pressure
   vx = par_getd("problem", "vx"); // ambient velocity field
   vy = par_getd("problem", "vy");
   W.Vy = vy ; 
   W.Vx = vx ;
   W.Vz = par_getd("problem", "vz");
   vxamp = par_getd("problem", "vxamp"); 
   vxy = par_getd("problem", "vxy"); // velocity shear (vy in x direction)
   noise = par_getd("problem", "noise");
#ifdef MHD
   Bx   = par_getd("problem", "bx");
   W.By = par_getd("problem", "by");
   W.Bz = par_getd("problem", "bz");
#endif
   g = par_getd("problem", "g"); // "gravity" (essentially, boolean, works if g>0)

   Wjet.d  = par_getd("problem", "djet"); // jet density
   Wjet.P  = par_getd("problem", "pjet"); // pressure inside the jet
   Wjet.Vx = par_getd("problem", "vxjet"); // velocity field along the axis
   Wjet.Vy = par_getd("problem", "vyjet"); // other velocity components
   Wjet.Vz = par_getd("problem", "vzjet");
#ifdef MHD
   Bxjet   = par_getd("problem", "bxjet");
   Wjet.By = par_getd("problem", "byjet");
   Wjet.Bz = par_getd("problem", "bzjet");
#endif

   rjet = par_getd("problem", "rjet"); // jet head radius
   
   U = Prim1D_to_Cons1D(&W, &Bx);
   Ujet = Prim1D_to_Cons1D(&Wjet, &Bxjet);

   x1_min = pDomain->RootMinX[0];
   x1_max = pDomain->RootMaxX[0];
   x2_min = pDomain->RootMinX[1];
   x2_max = pDomain->RootMaxX[1];
   x3_min = pDomain->RootMinX[2];
   x3_max = pDomain->RootMaxX[2];

   x1_mid = 0.5 * (x1_max + x1_min);
   x2_mid = 0.5 * (x2_max + x2_min);
   x3_mid = 0.5 * (x3_max + x3_min);

/* initialize index bounds assuming problem in xy plane */
   iu = ie + nghost;
   il = is - nghost;
   ju = je + nghost;
   jl = js - nghost;
   if(pGrid->Nx[2] > 1){
      ku = pGrid->ke + nghost;
      kl = pGrid->ks - nghost;
   }
   else{
      ku = pGrid->ke;
      kl = pGrid->ks;
   }

/* initialize conserved variables */
   
   for(k=kl; k<=ku; k++){
      for(j=jl; j<=ju; j++){
         for(i=il; i<=iu; i++){
	   cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	   W.Vx = vx + vxy * ((x2-x2_mid)/(x2_max - x2_min))*2.; // updating the vx
	   W.Vx = W.Vx/sqrt(1.+W.Vx*W.Vx)*vxamp;
	   W.Vy = vy * (1.+(ran2(&iseed) - 0.5)*noise);
	   W.d = d0 * (1.+(ran2(&iseed) - 0.5)*noise);
	   U = Prim1D_to_Cons1D(&W, &Bx);
	     
            pGrid->U[k][j][i].d  = U.d;
            pGrid->U[k][j][i].M1 = U.Mx;
            pGrid->U[k][j][i].M2 = U.My;
            pGrid->U[k][j][i].E  = U.E;
#ifdef MHD
            pGrid->U[k][j][i].B1c = Bx;
            pGrid->U[k][j][i].B2c = U.By;
            pGrid->U[k][j][i].B3c = U.Bz;
            pGrid->B1i[k][j][i] = Bx;
            pGrid->B2i[k][j][i] = U.By;
            pGrid->B3i[k][j][i] = U.Bz;
#endif
         }
      }
   }

   // grav. potential
   if(g>0.)
     {
       StaticGravPot = grav_pot2;
     }
   
/* Set boundary value function pointers */

   // if (pDomain->Disp[0] == 0) bvals_mhd_fun(pDomain,left_x1,jet_iib);

}


/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  int nl,nd;

  Wjet.d  = par_getd("problem", "djet");
  Wjet.P  = par_getd("problem", "pjet");
  Wjet.Vx = par_getd("problem", "vxjet");
  Wjet.Vy = par_getd("problem", "vyjet");
  Wjet.Vz = par_getd("problem", "vzjet");
#ifdef MHD
  Bxjet   = par_getd("problem", "bxjet");
  Wjet.By = par_getd("problem", "byjet");
  Wjet.Bz = par_getd("problem", "bzjet");
#endif

  rjet = par_getd("problem", "rjet");
  Ujet = Prim1D_to_Cons1D(&Wjet, &Bxjet);

  for (nl=0; nl<pM->NLevels; nl++){
  for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      if (pM->Domain[nl][nd].Disp[0] == 0) 
        bvals_mhd_fun(&(pM->Domain[nl][nd]),left_x1,jet_iib);
    }
  }}

  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}


/*===================== PRIVATE FUNCTIONS ====================================*/

/******************************************************************************/
/*! \fn void jet_iib(GridS *pGrid) 
 *  \brief Sets ghost zones to either outflow BC or
 * if cell is within jet radius, to jet values */
/******************************************************************************/

void jet_iib(GridS *pGrid){
  int i, is = pGrid->is;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real rad,x1,x2,x3;

  for(k=ks; k<=ke; k++){
    for(j=js; j<=je; j++){
      for(i=1; i<=nghost; i++){
        cc_pos(pGrid,(is-i),j,k,&x1,&x2,&x3);
        rad = sqrt(SQR(x2 - x2_mid) + SQR(x3 - x3_mid));
            
        if(rad <= rjet){
          pGrid->U[k][j][is-i].d  = Ujet.d;
          pGrid->U[k][j][is-i].M1 = Ujet.Mx;
          pGrid->U[k][j][is-i].M2 = Ujet.My;
          pGrid->U[k][j][is-i].M3 = Ujet.Mz;
          pGrid->U[k][j][is-i].E  = Ujet.E;
#ifdef MHD
          pGrid->U[k][j][is-i].B1c = Bxjet;
          pGrid->U[k][j][is-i].B2c = Ujet.By;
          pGrid->U[k][j][is-i].B3c = Ujet.Bz;
          pGrid->B1i[k][j][is-i] = Bxjet;
          pGrid->B2i[k][j][is-i] = Ujet.By;
          pGrid->B3i[k][j][is-i] = Ujet.Bz;
#endif
        } else{
          pGrid->U[k][j][is-i] = pGrid->U[k][j][is+(i-1)];
          pGrid->U[k][j][is-i].M1 = -pGrid->U[k][j][is-i].M1;
#ifdef MHD
          pGrid->U[k][j][is-i].B2c = -pGrid->U[k][j][is-i].B2c;
          pGrid->U[k][j][is-i].B3c = -pGrid->U[k][j][is-i].B3c;
          pGrid->B1i[k][j][is-i] =  pGrid->B1i[k][j][is+i];
          pGrid->B2i[k][j][is-i] = -pGrid->B2i[k][j][is+(i-1)];
          pGrid->B3i[k][j][is-i] = -pGrid->B3i[k][j][is+(i-1)];
#endif
        }
      }
    }
  }
}


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/*! \fn double ran2(long int *idum)
 *  \brief  Extracted from the Numerical Recipes in C (version 2) code. Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
 *
 * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1. 

 */
double ran2(long int *idum){
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX
/*----------------------------------------------------------------------------*/
/*! \fn static Real grav_pot2(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 0.1
 */

static Real grav_pot2(const Real x1, const Real x2, const Real x3)
{
  return 0.01*x2;
}
