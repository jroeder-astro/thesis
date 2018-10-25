/* file    tov.c
 * author  Frank Loeffler, converted from fortran thorn by Ian Hawke
 * date    2002/10/21
 * desc    TOV initial data
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "constants.h"

#include "tov.h"

#define NUMVARS 6

#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p (&vel_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p (&vel_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p (&vel_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velx_p_p (&vel_p_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely_p_p (&vel_p_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz_p_p (&vel_p_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

#define sx (&scon[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sy (&scon[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sz (&scon[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sx_p (&scon_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sy_p (&scon_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sz_p (&scon_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sx_p_p (&scon_p_p[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sy_p_p (&scon_p_p[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define sz_p_p (&scon_p_p[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

CCTK_REAL * TOV_Surface=0;
CCTK_REAL * TOV_R_Surface=0;
CCTK_REAL * TOV_RProp_Surface=0;

CCTK_REAL * TOV_r_1d=0;
CCTK_REAL * TOV_rbar_1d=0;
CCTK_REAL * TOV_press_1d=0;
CCTK_REAL * TOV_phi_1d=0;
CCTK_REAL * TOV_m_1d=0;
CCTK_REAL * TOV_mbary_1d=0;
CCTK_REAL * TOV_rprop_1d=0;
double myrho[10000];
double myp[10000];
double mymu[10000];
char line[500];
int neos;
FILE *fin;
CCTK_REAL getrho(double press); 
CCTK_REAL getmu(double press); 
CCTK_REAL getpress(double rho); 

#include "utils.inc"

void TOV_C_ParamCheck(CCTK_ARGUMENTS)                                          //adjusted
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  neos=0; 
/* read eos for simplicity, units are MeV/fm**3 and 1/fm**3 */
 fin = fopen("myeos.dat","r");
 while(!feof(fin))
 {
  fgets(line,300,fin);
  if(sscanf(line,"%lf %lf %lf",&mymu[neos],&myp[neos],&myrho[neos])!=3) break;
// mymu[neos]/=7.514e5*.4627;
// myp[neos]/=7.514e5*.4627;
// myrho[neos]*=939.5/7.514e5/.4627;
 }
  neos++;
 fclose(fin); 
  printf("Parameter Check!!!i");
    if (TOV_Solve_for_TOVs != 3) 
      {
	if (TOV_Solve_for_TOVs == 2)
	  {
	    CCTK_WARN(1, "TOV_Solve_for_TOVs is depreciated. "
                   "Use TOV_Enforce_Interpolation=\"yes\" instead.\n");
	    if (CCTK_ParameterSet("TOV_Enforce_Interpolation",
                            "TOVSolver",
                            "true"))
	      CCTK_WARN(0, "Error while steering this parameter.\n");
	    else
	      CCTK_WARN(1, "Steered this parameter for now accordingly.\n");
	  }
	else
	  CCTK_WARN(1, "TOV_Solve_for_TOVs is depreciated. "
                   "Use TOV_Enforce_Interpolation instead.\n");
      }
  if (TOV_ProperPosition)
  {
    if (TOV_Num_TOVs != 2)
      CCTK_WARN(0, "TOV_ProperPosition atm only works for TOV_Num_TOVs==2");
  }
}

/* centered differencing with one-sided differencing at the boundary */
#define DIFF_X(a) (((i==0)?(a[CCTK_GFINDEX3D(cctkGH, i+1, j, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i  , j, k)]): \
                    (i==(cctk_lsh[0]-1))?                          \
                           (a[CCTK_GFINDEX3D(cctkGH, i  , j, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i-1, j, k)]): \
                       0.5*(a[CCTK_GFINDEX3D(cctkGH, i+1, j, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i-1, j, k)]))/\
                   CCTK_DELTA_SPACE(0))
#define DIFF_Y(a) (((j==0)?(a[CCTK_GFINDEX3D(cctkGH, i, j+1, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j  , k)]): \
                    (j==(cctk_lsh[1]-1))?                          \
                           (a[CCTK_GFINDEX3D(cctkGH, i, j  , k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j-1, k)]): \
                       0.5*(a[CCTK_GFINDEX3D(cctkGH, i, j+1, k)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j-1, k)]))/\
                    CCTK_DELTA_SPACE(1))
#define DIFF_Z(a) (((k==0)?(a[CCTK_GFINDEX3D(cctkGH, i, j, k+1)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j, k  )]): \
                    (k==(cctk_lsh[2]-1))?                          \
                           (a[CCTK_GFINDEX3D(cctkGH, i, j, k  )] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j, k-1)]): \
                       0.5*(a[CCTK_GFINDEX3D(cctkGH, i, j, k+1)] - \
                            a[CCTK_GFINDEX3D(cctkGH, i, j, k-1)]))/\
                    CCTK_DELTA_SPACE(2))


/*@@
   @routine    TOV_Source_RHS
   @date       Thu Oct 24 14:30:00 2002
   @author     Frank Loeffler - converted fortran routine by Ian Hawke
   @desc
      The source terms for the ODEs. These are equations (2), (3), (4)
      and (18) from the Baumgarte notes.
      That is the vector in order is (P, m, phi, rbar).
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
@@*/
void TOV_C_Source_RHS(CCTK_REAL r, CCTK_REAL K, CCTK_REAL Gamma,
                      CCTK_REAL old_data[NUMVARS], CCTK_REAL source_data[NUMVARS]);
void TOV_C_Source_RHS(CCTK_REAL r, CCTK_REAL K, CCTK_REAL Gamma,
                      CCTK_REAL old_data[NUMVARS], CCTK_REAL source_data[NUMVARS])
{
  CCTK_REAL LOCAL_TINY, PI;                                                          
  CCTK_REAL press, rho, eps, mu, m;                                               // adjusted 
  CCTK_REAL r_minus_two_m;

  LOCAL_TINY = 1.0e-35;
  PI=4.0*atan(1.0);

  press           = old_data[0];
  if (press < LOCAL_TINY)
    press = LOCAL_TINY;
  m               = old_data[1];

  printf("HUHU\n");
  rho = pow(press / K, 1.0 / Gamma);
  eps = press / (Gamma - 1.0) / rho;
  mu  = rho * (1.0 + eps);
  rho=getrho(press);
  mu=getmu(press);
  eps=mu/rho-1.0;

  r_minus_two_m = r - 2.0 * m;

  if ((r<=0.0) && (m<=0.0))
  {
    source_data[1] = 0.0;
    source_data[2] = 0.0;
    source_data[3] = 0.0;
    source_data[4] = 0.0;
    source_data[5] = 0.0;
  }
  else
  {
    source_data[2] = (m + 4*PI * r*r*r * press) / r_minus_two_m / r;
    /* source_data[0] = -(press + mu) * source_data[2]; */
    source_data[0] = -(press + mu) *
                     (m + 4*PI * r*r*r * press) / r_minus_two_m / r;
    source_data[1] = 4*PI * r*r * mu;
    source_data[3] = (sqrt(r) - sqrt(r_minus_two_m)) / r / sqrt(r_minus_two_m);
    source_data[5] = 1.0/sqrt(1.0-2.0*m/r);
    source_data[4] = source_data[5] * 4*PI * rho * r*r;
  }
}

/*@@
   @routine    TOV_Integrate_RHS
   @date       Thu Oct 24 14:30:00 2002
   @author     Frank Loeffler, converted fortran routine by Ian Hawke
   @desc
      Integrates the ODEs using RK4.
      We rescale at the end to match to a Schwarzschild exterior.
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
@@*/
void TOV_C_Integrate_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL LOCAL_TINY;

  CCTK_INT star, star_i, i, TOV_Surface_Index;
  CCTK_REAL old_data[NUMVARS], source_data[NUMVARS],
            in_data[NUMVARS], new_data[NUMVARS],
            k1[NUMVARS], k2[NUMVARS], k3[NUMVARS], k4[NUMVARS];
  CCTK_REAL Surface_Mass, factor, local_rho;

  LOCAL_TINY = 1.0e-20;
  printf("\nIntegrate RHS!!!\n");
  assert(TOV_Surface!=0);
  assert(TOV_R_Surface!=0);
  assert(TOV_RProp_Surface!=0);

  assert(TOV_r_1d!=0);
  assert(TOV_rbar_1d!=0);
  assert(TOV_press_1d!=0);
  assert(TOV_phi_1d!=0);
  assert(TOV_m_1d!=0);
  assert(TOV_mbary_1d!=0);
  assert(TOV_rprop_1d!=0);

  /* do it for all stars */
  for (star=0; star < TOV_Num_TOVs; star++)
  {
    /* remember array index */
    star_i = star * TOV_Num_Radial;
    const CCTK_REAL rho_central=TOV_Rho_Central[star];
    printf("HALLO %f",rho_central);                                            //adjusted
    /* Set conformal state like set in parameter file if we do not use
     * the old initial data. In this case we have to use what we get */
    if (!TOV_Use_Old_Initial_Data)
      if(CCTK_EQUALS(metric_type, "static conformal"))
      {
        *conformal_state=1;
        if (CCTK_EQUALS(conformal_storage,"factor+derivs"))
          *conformal_state = 2;
        else if (CCTK_EQUALS(conformal_storage,"factor+derivs+2nd derivs"))
          *conformal_state = 3;
        CCTK_VInfo(CCTK_THORNSTRING, "conformal_state set to %d",
                   (int)*conformal_state);
      }

    /* clear arrays first */
    TOV_C_fill(&(TOV_press_1d[star_i]), TOV_Num_Radial, 0.0);
    TOV_C_fill(&(TOV_m_1d    [star_i]), TOV_Num_Radial, 0.0);
    TOV_C_fill(&(TOV_phi_1d  [star_i]), TOV_Num_Radial, 0.0);
    TOV_C_fill(&(TOV_rbar_1d [star_i]), TOV_Num_Radial, 0.0);
    TOV_C_fill(&(TOV_r_1d    [star_i]), TOV_Num_Radial, 0.0);
    TOV_C_fill(&(TOV_mbary_1d[star_i]), TOV_Num_Radial, 0.0);
    TOV_C_fill(&(TOV_rprop_1d[star_i]), TOV_Num_Radial, 0.0);

    /* set start values */
    TOV_press_1d[star_i] = TOV_K *
                            pow(rho_central, TOV_Gamma);
    TOV_press_1d[star_i] = getpress(rho_central);                             //adjusted
    /* TOV_r_1d    [star_i] = LOCAL_TINY;
    TOV_rbar_1d [star_i] = LOCAL_TINY;*/

    /* build TOV_r_1d[] */
    for (i=star_i+1; i < star_i+TOV_Num_Radial; i++)
      TOV_r_1d[i] = TOV_r_1d[i-1] + TOV_dr[star];

    TOV_Surface[star] = -1.0;
    TOV_Surface_Index = -1.0;

#define RKLOOP for (int rk=0; rk<NUMVARS; rk++)
    /* loop over all radii */
    for (i=star_i; (i < star_i+TOV_Num_Radial-1) &&
                            (TOV_Surface[star] < 0.0); i++)
    {
      /* set up RK arrays */
      old_data[0] = TOV_press_1d[i];
      old_data[1] = TOV_m_1d[i];
      old_data[2] = TOV_phi_1d[i];
      if (fabs(TOV_rbar_1d[i] - TOV_r_1d[i]) < LOCAL_TINY)
        old_data[3] = 0.0;
      else
        old_data[3] = log(TOV_rbar_1d[i] / TOV_r_1d[i]);
      old_data[4] = TOV_mbary_1d[i];
      old_data[5] = TOV_rprop_1d[i];

      /* usual RK4 */
      RKLOOP in_data[rk] = old_data[rk];
      TOV_C_fill(source_data, 6, 0.0);

      TOV_C_Source_RHS(TOV_r_1d[i],
                     TOV_K, TOV_Gamma,
                     in_data, source_data);

      RKLOOP k1[rk] = TOV_dr[star] * source_data[rk];
      RKLOOP in_data[rk] = old_data[rk] + 0.5 * k1[rk];
      TOV_C_Source_RHS(TOV_r_1d[i]+ 0.5 * TOV_dr[star],
                       TOV_K, TOV_Gamma,
                       in_data, source_data);

      RKLOOP k2[rk] = TOV_dr[star] * source_data[rk];
      RKLOOP in_data[rk] = old_data[rk] + 0.5 * k2[rk];
      TOV_C_Source_RHS(TOV_r_1d[i]+ 0.5 * TOV_dr[star],
                       TOV_K, TOV_Gamma,
                       in_data, source_data);

      RKLOOP k3[rk] = TOV_dr[star] * source_data[rk];
      RKLOOP in_data[rk] = old_data[rk] + k3[rk];
      TOV_C_Source_RHS(TOV_r_1d[i]+ TOV_dr[star],
                       TOV_K, TOV_Gamma,
                       in_data, source_data);
      RKLOOP k4[rk] = TOV_dr[star] * source_data[rk];
      RKLOOP new_data[rk] = old_data[rk] + (k1[rk] + k4[rk] + 2.0 * (k2[rk] + k3[rk])) /6.0;

      TOV_press_1d[i+1] = new_data[0];
      TOV_m_1d    [i+1] = new_data[1];
      TOV_phi_1d  [i+1] = new_data[2];
      TOV_rbar_1d [i+1] = TOV_r_1d[i+1] * exp(new_data[3]);
      TOV_mbary_1d[i+1] = new_data[4];
      TOV_rprop_1d[i+1] = new_data[5];

      /* otherwise the code crashes later */
      if (TOV_press_1d[i+1] < 0.0)
          TOV_press_1d[i+1] = 0.0;

      local_rho = pow(TOV_press_1d[i+1] / TOV_K, 1.0 / TOV_Gamma);
      local_rho = getrho(TOV_press_1d[i+1]);                                       // adjusted

      /* scan for the surface */
      if ( (local_rho <= 0.0) ||
           (TOV_press_1d[i+1] <= 0.0) )
      {
        TOV_Surface[star]   = TOV_r_1d[i];
        TOV_R_Surface[star] = TOV_rbar_1d[i];
        TOV_RProp_Surface[star] = TOV_rprop_1d[i];
        TOV_Surface_Index = i;
      }
    }
    if (TOV_Surface[star] < 0.0)
      CCTK_WARN(0, "Did not integrate out to surface of the star! "
                   "Increase TOV_dr or TOV_Num_Radial and rerun");

    Surface_Mass = TOV_m_1d[TOV_Surface_Index];
    factor = 0.5 * (sqrt(TOV_Surface[star] *
                         (TOV_Surface[star] - 2.00 * Surface_Mass)) +
                    TOV_Surface[star] - Surface_Mass) /
                    TOV_rbar_1d[TOV_Surface_Index];

    TOV_R_Surface[star] *= factor;
    for (i=star_i; i < star_i+TOV_Num_Radial; i++)
    {
      TOV_rbar_1d[i] *= factor;
      TOV_phi_1d[i]  -= TOV_phi_1d[TOV_Surface_Index] -
                        0.5 * log(1.0 - 2.0 * Surface_Mass / TOV_Surface[star]);
      /* match to Schwarzschield */
      if (i > TOV_Surface_Index)
      {
        TOV_press_1d[i] = 0.0;
        TOV_rbar_1d [i] = 0.5 *
                          (sqrt(TOV_r_1d[i]*(TOV_r_1d[i] - 2.0*Surface_Mass)) +
                          TOV_r_1d[i] - Surface_Mass);
        TOV_m_1d[i]     = Surface_Mass;
        TOV_phi_1d[i]   = 0.5 * log( 1.0 - 2.0 * Surface_Mass / TOV_r_1d[i]);
        TOV_mbary_1d[i] = TOV_mbary_1d[TOV_Surface_Index];
      }
    }
  }
  CCTK_INFO("Integrated TOV equation");
  /* do some info */
  CCTK_VInfo(CCTK_THORNSTRING, "Information about the TOVs used:");
  CCTK_VInfo("", "TOV    radius    mass  bary_mass mass(g) cent.rho rho(cgi)        K   K(cgi)    Gamma");
  for (i=0; i<TOV_Num_TOVs; i++)
    if (fabs(TOV_Gamma - 2.0) < LOCAL_TINY)
      CCTK_VInfo("","  %d  %8g %8g %8g %8.3g %8g %8.3g %8g %8.3g %8g",
                 (int)i+1, TOV_R_Surface[i],
                 TOV_m_1d[(i+1)*TOV_Num_Radial-1],
                 TOV_mbary_1d[(i+1)*TOV_Num_Radial-1],
                 TOV_m_1d[(i+1)*TOV_Num_Radial-1]*CONSTANT_Msolar_cgi,
                 TOV_Rho_Central[i],
                 TOV_Rho_Central[i]/pow(CONSTANT_G_cgi,3.0)/
                                    pow(CONSTANT_Msolar_cgi,2.0)*
                                    pow(CONSTANT_c_cgi,6.0),
                 TOV_K,
                 TOV_K*pow(CONSTANT_G_cgi,3.0)*
                          pow(CONSTANT_Msolar_cgi,2.0)/
                          pow(CONSTANT_c_cgi,4.0),
                 TOV_Gamma);
    else
      CCTK_VInfo("","  %d  %8g %8g %8.3g %8g %8.3g %8g %8g",
                 (int)i+1, TOV_R_Surface[i],
                 TOV_m_1d[(i+1)*TOV_Num_Radial-1],
                 TOV_m_1d[(i+1)*TOV_Num_Radial-1]*CONSTANT_Msolar_cgi,
                 TOV_Rho_Central[i],
                 TOV_Rho_Central[i]/pow(CONSTANT_G_cgi,3.0)/
                                    pow(CONSTANT_Msolar_cgi,2.0)*
                                    pow(CONSTANT_c_cgi,6.0),
                 TOV_K, TOV_Gamma);

}

/*----------------------------------------------------------------------------*/

/* utility routine
 * recursive search-routine for arrays
 * here used to look for the last index in an ordered array with its
 * value < goal
 */
CCTK_INT TOV_C_find_index(CCTK_INT   array_size,
                          CCTK_REAL *array,
                          CCTK_REAL  goal,
                          CCTK_INT   lower_index,
                          CCTK_INT   upper_index);
CCTK_INT TOV_C_find_index(CCTK_INT   array_size,
                          CCTK_REAL *array,
                          CCTK_REAL  goal,
                          CCTK_INT   lower_index,
                          CCTK_INT   upper_index)
{
  CCTK_INT middle_index;

  if (lower_index >= (upper_index-1))
      return lower_index;

  middle_index = (lower_index + upper_index) /2;

  if (array[middle_index] < goal)
    return TOV_C_find_index(array_size, array, goal, middle_index, upper_index);
  else
    return TOV_C_find_index(array_size, array, goal, lower_index, middle_index);
}

/* utility rountine
 * interpolates from (thorn-internal) 1D-data to Cactus 3D-grid */
/* input is all but *press_point *phi_point and *r_point */
void TOV_C_interp_tov_isotropic(
                                CCTK_INT  star,
                                CCTK_REAL *TOV_press_1d_local,
                                CCTK_REAL *TOV_phi_1d_local,
                                CCTK_REAL *TOV_rbar_1d_local,
                                CCTK_REAL *TOV_r_1d_local,
                                CCTK_REAL *r,
                                CCTK_REAL surface,
                                CCTK_REAL *press_point,
                                CCTK_REAL *phi_point,
                                CCTK_REAL *r_point);
void TOV_C_interp_tov_isotropic(
                                CCTK_INT  star,
                                CCTK_REAL *TOV_press_1d_local,
                                CCTK_REAL *TOV_phi_1d_local,
                                CCTK_REAL *TOV_rbar_1d_local,
                                CCTK_REAL *TOV_r_1d_local,
                                CCTK_REAL *r,
                                CCTK_REAL surface,
                                CCTK_REAL *press_point,
                                CCTK_REAL *phi_point,
                                CCTK_REAL *r_point)
{
  DECLARE_CCTK_PARAMETERS
  CCTK_INT  left_index;
  CCTK_REAL h, M;

  if (*r < 0.0)
    CCTK_WARN(0, "Negative radius found");
  if (*r < TOV_rbar_1d_local[1])
    *r=TOV_rbar_1d_local[1];
  if (*r > TOV_rbar_1d_local[TOV_Num_Radial-2])
    {
	{
	  *press_point= 0.0;
	  M = 0.5 * TOV_r_1d_local[TOV_Num_Radial-1] *
	    (1.0 - exp(2.0*TOV_phi_1d_local[TOV_Num_Radial-1]));
	  *r_point=(2* *r+M)*(2* *r+M)*0.25/ *r;
	  *phi_point=0.5*log(1-2*M/ *r_point);
	  return;
	}
    }

  if (TOV_Fast_Interpolation)
    left_index = TOV_C_find_index(TOV_Num_Radial-1, TOV_rbar_1d_local, *r, 0,
                                  TOV_Num_Radial-1);
  else
  {
    left_index=0;
    while( (left_index < TOV_Num_Radial-2) &&
           (TOV_rbar_1d_local[left_index+1] < *r) )
      left_index++;
  }

  h = (*r - TOV_rbar_1d_local[left_index]) /
      (TOV_rbar_1d_local[left_index+1] - TOV_rbar_1d_local[left_index]);
  *r_point     = (1.0 - h) * TOV_r_1d_local[left_index] +
                        h  * TOV_r_1d_local[left_index+1];
  *phi_point   = (1.0 - h) * TOV_phi_1d_local[left_index] +
                        h  * TOV_phi_1d_local[left_index+1];
  if (*r_point < surface)
    *press_point = (1.0 - h) * TOV_press_1d_local[left_index] +
                          h  * TOV_press_1d_local[left_index+1];
  else
    *press_point = 0.0;
}

/*@@
   @routine    TOV_Exact
   @date       Thu Oct 24 14:30:00 2002
   @author     Frank Loeffler, converted fortran routine by Ian Hawke
   @desc
       Schedule routine for interpolation of 1D to 3D grid
   @enddesc
   @calls
   @calledby
   @history
   @endhistory
@@*/
void TOV_C_Exact(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL *press_point, *rho_point, *eps_point,
            *mu_point, *phi_point, *r_point;

  CCTK_INT  LSH_MAX_I;

  CCTK_INT i,j,k, i3D, star;
  CCTK_REAL *r_to_star;
  CCTK_REAL g_diag, max_g_diag, max_rho;
  CCTK_REAL my_velx, my_vely, my_velz, my_psi4;
  CCTK_REAL PI, local_tiny;

  CCTK_INT tov_lapse, tov_shift;

  tov_lapse = CCTK_EQUALS(initial_lapse, "tov");
  tov_shift = CCTK_EQUALS(initial_shift, "tov");

  PI=4.0*atan(1.0);
  local_tiny=1.0e-14;
  /* remember index of last member of array */
  LSH_MAX_I = CCTK_GFINDEX3D(cctkGH,
                             cctk_lsh[0]-1, cctk_lsh[1]-1, cctk_lsh[2]-1);

  assert(TOV_Surface!=0);
  assert(TOV_R_Surface!=0);

  assert(TOV_r_1d!=0);
  assert(TOV_rbar_1d!=0);
  assert(TOV_press_1d!=0);
  assert(TOV_phi_1d!=0);
  assert(TOV_m_1d!=0);

  /* allocate local arrays */
  r_to_star   = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));
  press_point = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));
  rho_point   = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));
  eps_point   = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));
  mu_point    = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));
  phi_point   = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));
  r_point     = (CCTK_REAL *) calloc (TOV_Num_TOVs, sizeof(CCTK_REAL));

  /* clear initial data */
  if (TOV_Clear_Initial_Data > 0 && !(TOV_Use_Old_Initial_Data))
  {
    TOV_C_fill(kxx, LSH_MAX_I+1, 0.0);
    TOV_C_fill(kxy, LSH_MAX_I+1, 0.0);
    TOV_C_fill(kxz, LSH_MAX_I+1, 0.0);
    TOV_C_fill(kyy, LSH_MAX_I+1, 0.0);
    TOV_C_fill(kyz, LSH_MAX_I+1, 0.0);
    TOV_C_fill(kzz, LSH_MAX_I+1, 0.0);
    TOV_C_fill(gxx, LSH_MAX_I+1, 0.0);
    TOV_C_fill(gyy, LSH_MAX_I+1, 0.0);
    TOV_C_fill(gzz, LSH_MAX_I+1, 0.0);
    TOV_C_fill(gxy, LSH_MAX_I+1, 0.0);
    TOV_C_fill(gxz, LSH_MAX_I+1, 0.0);
    TOV_C_fill(gyz, LSH_MAX_I+1, 0.0);
    TOV_C_fill(alp, LSH_MAX_I+1, 1.0);

    if (*shift_state != 0)
    {
      TOV_C_fill(betax, LSH_MAX_I+1, 0.0);
      TOV_C_fill(betay, LSH_MAX_I+1, 0.0);
      TOV_C_fill(betaz, LSH_MAX_I+1, 0.0);
    }

    if (*conformal_state != 0)
    {
      TOV_C_fill(psi, LSH_MAX_I+1, 1.0);
      if (*conformal_state > 1)
      {
        TOV_C_fill(psix, LSH_MAX_I+1, 0.0);
        TOV_C_fill(psiy, LSH_MAX_I+1, 0.0);
        TOV_C_fill(psiz, LSH_MAX_I+1, 0.0);
        if (*conformal_state > 2)
        {
          TOV_C_fill(psixx, LSH_MAX_I+1, 0.0);
          TOV_C_fill(psixy, LSH_MAX_I+1, 0.0);
          TOV_C_fill(psixz, LSH_MAX_I+1, 0.0);
          TOV_C_fill(psiyy, LSH_MAX_I+1, 0.0);
          TOV_C_fill(psiyz, LSH_MAX_I+1, 0.0);
          TOV_C_fill(psizz, LSH_MAX_I+1, 0.0);
        }
      }
    }
  }
  if (!TOV_Use_Old_Matter_Initial_Data)
  {
    CCTK_INFO("Not using old matter initial data");
    TOV_C_fill(rho,        LSH_MAX_I+1, 0.0);
    TOV_C_fill(eps,        LSH_MAX_I+1, 0.0);
    TOV_C_fill(press,      LSH_MAX_I+1, 0.0);
    TOV_C_fill(w_lorentz,  LSH_MAX_I+1, 0.0);
    TOV_C_fill(velx,       LSH_MAX_I+1, 0.0);
    TOV_C_fill(vely,       LSH_MAX_I+1, 0.0);
    TOV_C_fill(velz,       LSH_MAX_I+1, 0.0);
  }
  /* use the fast interpolation? only useful for testing this */
  if (TOV_Fast_Interpolation == 0)
      CCTK_INFO("Interpolating the slow way.");

  /* go over all 3D-grid points */
  for(i=0; i<cctk_lsh[0]; i++)
   for(j=0; j<cctk_lsh[1]; j++)
    for(k=0; k<cctk_lsh[2]; k++)
    {
      i3D=CCTK_GFINDEX3D(cctkGH, i, j, k);
      /* remember the old conformal factor to the power of 4 */
      if (*conformal_state != 0)
        my_psi4=pow(psi[i3D], 4.0);
      else
        my_psi4=1.0;

      for (star=0; star<TOV_Num_TOVs; star++)
      {
        r_to_star[star] =
          sqrt( (x[i3D]-TOV_Position_x[star]) *
                (x[i3D]-TOV_Position_x[star]) +
                (y[i3D]-TOV_Position_y[star]) *
                (y[i3D]-TOV_Position_y[star]) +
                (z[i3D]-TOV_Position_z[star]) *
                (z[i3D]-TOV_Position_z[star]) );
        int star_i = star * TOV_Num_Radial;

        /* do the actual interpolation */
        TOV_C_interp_tov_isotropic(star,
                               &(TOV_press_1d[star_i]), &(TOV_phi_1d[star_i]),
                               &(TOV_rbar_1d[star_i]), &(TOV_r_1d[star_i]),
                               &(r_to_star[star]), TOV_Surface[star],
                               &(press_point[star]),
                               &(phi_point[star]), &(r_point[star]));

        /* is some perturbation wanted? */
        if (Perturb[star] == 0)
        {
          rho_point[star] = pow(press_point[star]/TOV_K,
                                1.0/TOV_Gamma);
          rho_point[star] = getrho(press_point[star]);                          //adjusted
        }
        else
        {
          rho_point[star] = pow(press_point[star]/TOV_K,
                                1.0/TOV_Gamma) *
                            (1.0 +
                             Pert_Amplitude[star] *
                               cos(PI/2.0 * r[i3D] / TOV_R_Surface[star]));
          rho_point[star] = getrho(press_point[star])*                          //adjusted
                            (1.0 +
                             Pert_Amplitude[star] *
                               cos(PI/2.0 * r[i3D] / TOV_R_Surface[star]));
        }
        if (rho_point[star] > local_tiny)
        {
          eps_point[star] = press_point[star] / (TOV_Gamma - 1.0)
                                              /  rho_point[star];
          eps_point[star] = getmu(press_point[star])/rho_point[star]-1.0;        //adjusted     
        }
        else
          eps_point[star] = 0.0;
        mu_point[star]  = rho_point[star] * (1.0 + eps_point[star]);
      }
      /* find out from which star we want to have the data */
      if (CCTK_EQUALS(TOV_Combine_Method, "maximum"))
      {
        /* to do this, we use here simply the max of the gxx-value */
        star=0;
        max_g_diag = 0.0;
        max_rho = rho_point[0];
        for (int star_i=0; star_i<TOV_Num_TOVs; star_i++)
        {
          g_diag = (r_point[star_i] / (r_to_star[star_i] + 1.0e-30)) *
                   (r_point[star_i] / (r_to_star[star_i] + 1.0e-30));
          if ((g_diag - max_g_diag) > local_tiny)
          {
            max_g_diag=g_diag;
            star=star_i;
          }
          if ((rho_point[star_i] - max_rho) > local_tiny)
          {
            max_rho=rho_point[star_i];
            star=star_i;
          }
        }
        /* handle initial data */
        if (TOV_Use_Old_Initial_Data)
        {
          /* check metric */
          if ((my_psi4 * gxx[i3D] < max_g_diag) &&
              (my_psi4 * gyy[i3D] < max_g_diag) &&
              (my_psi4 * gzz[i3D] < max_g_diag))
          {
            if (TOV_Conformal_Flat_Three_Metric)
            {
              psi[i3D] = pow(max_g_diag, 0.25);
              my_psi4 = max_g_diag;
            }
            else
            {
              gxx[i3D] = max_g_diag/my_psi4;
              gyy[i3D] = max_g_diag/my_psi4;
              gzz[i3D] = max_g_diag/my_psi4;
              gxy[i3D] = gxz[i3D] = gyz[i3D] = 0.0;
            }
          }
          /* check matter */
          if (TOV_Use_Old_Matter_Initial_Data)
          {
            if (rho[i3D] > max_rho)
            {
              /* we do not need this array element anymore, since we use
               * the initial data, so lets use it */
              star=0;
              max_rho          =rho[i3D];
              eps_point[star]  =eps[i3D];
              press_point[star]=press[i3D];
              my_velx=velx[i3D];
              my_vely=vely[i3D];
              my_velz=velz[i3D];
            }
            else
            {
              if (tov_lapse)
                alp[i3D] = exp(phi_point[star]);
              if (tov_shift)
              {
                betax[i3D] = 0.0;
                betay[i3D] = 0.0;
                betaz[i3D] = 0.0;
              }
              my_velx=TOV_Velocity_x[star];
              my_vely=TOV_Velocity_y[star];
              my_velz=TOV_Velocity_z[star];
            }
          }
          else
          {
            if (tov_lapse)
              alp[i3D] = exp(phi_point[star]);
            if (tov_shift)
            {
              betax[i3D] = 0.0;
              betay[i3D] = 0.0;
              betaz[i3D] = 0.0;
            }
            my_velx=TOV_Velocity_x[star];
            my_vely=TOV_Velocity_y[star];
            my_velz=TOV_Velocity_z[star];
          }
        }
        else /* do not use old initial data */
        {
          /* no psi, since it is 1.0 here */
          /* but maybe we want to have it != 1.0 */
          if (TOV_Conformal_Flat_Three_Metric)
          {
            psi[i3D] = pow(max_g_diag, 0.25);
            my_psi4 = max_g_diag;
            gxx[i3D] = gyy[i3D] = gzz[i3D] = 1.0;
            gxy[i3D] = gxz[i3D] = gyz[i3D] = 0.0;
          }
          else
          {
            gxx[i3D] = max_g_diag;
            gyy[i3D] = max_g_diag;
            gzz[i3D] = max_g_diag;
            gxy[i3D] = gxz[i3D] = gyz[i3D] = 0.0;
          }
          if (tov_lapse)
            alp[i3D] = exp(phi_point[star]);
          if (tov_shift)
          {
            betax[i3D] = 0.0;
            betay[i3D] = 0.0;
            betaz[i3D] = 0.0;
          }
          my_velx=TOV_Velocity_x[star];
          my_vely=TOV_Velocity_y[star];
          my_velz=TOV_Velocity_z[star];
        }

        /* set to defined velocity. default is 0.0 because other velocities
         * violate Einsteins equations */
        velx[i3D] = my_velx;
        vely[i3D] = my_vely;
        velz[i3D] = my_velz;

        w_lorentz[i3D] = 1/sqrt(1.0-(
                                  gxx[i3D] * velx[i3D] * velx[i3D]+
                                  gyy[i3D] * vely[i3D] * vely[i3D]+
                                  gzz[i3D] * velz[i3D] * velz[i3D]+
                                2*gxy[i3D] * velx[i3D] * vely[i3D]+
                                2*gxz[i3D] * velx[i3D] * velz[i3D]+
                                2*gyz[i3D] * vely[i3D] * velz[i3D])*
                                my_psi4);

        rho[i3D] = max_rho;
        eps[i3D] = eps_point[star];
        press[i3D] = press_point[star];

      }
      else if (CCTK_EQUALS(TOV_Combine_Method, "average"))
      {
        /* here we 'average' all values in a more intelligent way */
        if (TOV_Use_Old_Matter_Initial_Data)
          max_rho=rho[i3D];
        else
        {
          max_rho=0.0;
          rho[i3D] = 0.0;
        }
        star=-1;
        for (int star_i=0; star_i<TOV_Num_TOVs; star_i++)
        {
          if (tov_lapse)
            alp[i3D] *= exp(phi_point[star_i]);
          if (tov_shift)
          {
            betax[i3D] = 0.0;
            betay[i3D] = 0.0;
            betaz[i3D] = 0.0;
          }
          if (TOV_Conformal_Flat_Three_Metric)
          {
            /* This is a hack, since it does not check if the input data is
             * really conformally flat. It simply assumes this by only using
             * gxx */
            my_psi4 = (r_point[star_i] * r_point[star_i] /
                       (r_to_star[star_i] * r_to_star[star_i] + 1.0e-30)) /
                      my_psi4 + pow(psi[i3D], 4.0) * gxx[i3D];
            psi[i3D] = pow(my_psi4, 0.25);
            if (!TOV_Use_Old_Initial_Data)
            {
              gxx[i3D] = gyy[i3D] = gzz[i3D] = 1.0;
              gxy[i3D] = gxz[i3D] = gyz[i3D] = 0.0;
            }
          }
          else
            gxx[i3D] += (r_point[star_i] * r_point[star_i] /
                         (r_to_star[star_i] * r_to_star[star_i] + 1.0e-30)) /
                        my_psi4;
          rho[i3D] += rho_point[star_i];
          eps[i3D] += eps_point[star_i];
          press[i3D] += press_point[star_i];
          /* we still have to know if we are inside one star - and which */
          if (rho_point[star_i] > max_rho)
          {
            max_rho=rho_point[star_i];
            star=star_i;
          }
        }

        if (TOV_Conformal_Flat_Three_Metric)
        {
          my_psi4 -= ((TOV_Num_TOVs+TOV_Use_Old_Initial_Data-1)/my_psi4);
          psi[i3D] = pow(my_psi4, 0.25);
        }
        else
        {
          gxx[i3D] -= ((TOV_Num_TOVs+TOV_Use_Old_Initial_Data-1)/my_psi4);
          gyy[i3D] = gxx[i3D];
          gzz[i3D] = gxx[i3D];
        }

        /* set to defined velocity. default is 0.0 because other velocities
         * violate the constraints */
        if (star > -1)
        {
          velx[i3D] = TOV_Velocity_x[star];
          vely[i3D] = TOV_Velocity_y[star];
          velz[i3D] = TOV_Velocity_z[star];
        }

        w_lorentz[i3D] = 1/sqrt(1.0-(
                                  gxx[i3D] * velx[i3D] * velx[i3D]+
                                  gyy[i3D] * vely[i3D] * vely[i3D]+
                                  gzz[i3D] * velz[i3D] * velz[i3D]+
                                2*gxy[i3D] * velx[i3D] * vely[i3D]+
                                2*gxz[i3D] * velx[i3D] * velz[i3D]+
                                2*gyz[i3D] * vely[i3D] * velz[i3D]) * my_psi4);
      }
    }


  /* if used, recalculate the derivatives of the conformal factor */
  if (*conformal_state > 1)
    /* go again over all 3D-grid points */
    for(i=0; i<cctk_lsh[0]; i++)
     for(j=0; j<cctk_lsh[1]; j++)
      for(k=0; k<cctk_lsh[2]; k++)
      {
        i3D=CCTK_GFINDEX3D(cctkGH, i, j, k);
        psix[i3D]=(((i==0)?
                    (psi[CCTK_GFINDEX3D(cctkGH, i+1, j, k)] -
                     psi[CCTK_GFINDEX3D(cctkGH, i  , j, k)]):
                   (i==(cctk_lsh[0]-1))?
                    (psi[CCTK_GFINDEX3D(cctkGH, i  , j, k)] -
                     psi[CCTK_GFINDEX3D(cctkGH, i-1, j, k)]):
                0.5*(psi[CCTK_GFINDEX3D(cctkGH, i+1, j, k)] -
                     psi[CCTK_GFINDEX3D(cctkGH, i-1, j, k)])));
        psix[i3D] = DIFF_X(psi);
        psiy[i3D] = DIFF_Y(psi);
        psiz[i3D] = DIFF_Z(psi);
      }
  if (*conformal_state > 2)
    /* go again over all 3D-grid points */
    for(i=0; i<cctk_lsh[0]; i++)
     for(j=0; j<cctk_lsh[1]; j++)
      for(k=0; k<cctk_lsh[2]; k++)
      {
        i3D=CCTK_GFINDEX3D(cctkGH, i, j, k);
        psixx[i3D] = DIFF_X(psix)/psi[i3D];
        psiyy[i3D] = DIFF_Y(psiy)/psi[i3D];
        psizz[i3D] = DIFF_Z(psiz)/psi[i3D];
        psixy[i3D] = DIFF_X(psiy)/psi[i3D];
        psiyz[i3D] = DIFF_Y(psiz)/psi[i3D];
        psixz[i3D] = DIFF_Z(psix)/psi[i3D];
      }
  if (*conformal_state > 1)
    /* go again over all 3D-grid points */
    for(i=0; i<cctk_lsh[0]; i++)
     for(j=0; j<cctk_lsh[1]; j++)
      for(k=0; k<cctk_lsh[2]; k++)
      {
        i3D=CCTK_GFINDEX3D(cctkGH, i, j, k);
        psix[i3D] /= psi[i3D];
        psiy[i3D] /= psi[i3D];
        psiz[i3D] /= psi[i3D];
      }
  i3D = cctk_lsh[2]*cctk_lsh[1]*cctk_lsh[0];
  switch(TOV_Populate_Timelevels)
  {
    case 3:
        TOV_Copy(i3D, gxx_p_p,  gxx);
        TOV_Copy(i3D, gyy_p_p,  gyy);
        TOV_Copy(i3D, gzz_p_p,  gzz);
        TOV_Copy(i3D, gxy_p_p,  gxy);
        TOV_Copy(i3D, gxz_p_p,  gxz);
        TOV_Copy(i3D, gyz_p_p,  gyz);
        TOV_Copy(i3D, rho_p_p,  rho);
        TOV_Copy(i3D, eps_p_p,  eps);
        TOV_Copy(i3D, velx_p_p, velx);
        TOV_Copy(i3D, vely_p_p, vely);
        TOV_Copy(i3D, velz_p_p, velz);
        TOV_Copy(i3D, w_lorentz_p_p, w_lorentz);
        // fall through
    case 2:
        TOV_Copy(i3D, gxx_p,  gxx);
        TOV_Copy(i3D, gyy_p,  gyy);
        TOV_Copy(i3D, gzz_p,  gzz);
        TOV_Copy(i3D, gxy_p,  gxy);
        TOV_Copy(i3D, gxz_p,  gxz);
        TOV_Copy(i3D, gyz_p,  gyz);
        TOV_Copy(i3D, rho_p,  rho);
        TOV_Copy(i3D, eps_p,  eps);
        TOV_Copy(i3D, velx_p, velx);
        TOV_Copy(i3D, vely_p, vely);
        TOV_Copy(i3D, velz_p, velz);
        TOV_Copy(i3D, w_lorentz_p, w_lorentz);
        // fall through
    case 1:
        break;
    default:
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Unsupported number of TOV_Populate_TimelevelsL: %d",
                   (int)TOV_Populate_Timelevels);
        break;
  }
  CCTK_INFO("Done interpolation.");

  /* free local arrays */
  free(r_to_star);
  free(press_point);
  free(rho_point);
  free(eps_point);
  free(mu_point);
  free(phi_point);
  free(r_point);
}

void TOV_Prepare_Fake_Evolution(CCTK_ARGUMENTS)
{
    if (CCTK_ParameterSet("TOV_Populate_Timelevels",
                          "TOVSolver",
                          "1"))
        CCTK_WARN(0,
                  "Could not prepare for fake evolution - steering failed\n");
}

inline static CCTK_REAL calc_coord_dist(CCTK_REAL prop_dist, CCTK_REAL M, CCTK_REAL R)
{
    CCTK_REAL xnp = prop_dist;
    CCTK_REAL xn, f, fp;
    do {
      xn = xnp;
      f  = xn - prop_dist + M * log(xn/R - 1);
      fp = 1 + M / (xn-R);
      xnp = xn - f/fp;
    } while (fabs((xnp-xn)/xn) > 10.e-14);
    return xnp;
}

/* Only works for equal-mass binary NS systems atm */
void TOV_Set_ProperPositions(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  /* The specified parameters are given in proper distance */
  CCTK_REAL prop_dist = sqrt((TOV_Position_x[0]-TOV_Position_x[1])*
                             (TOV_Position_x[0]-TOV_Position_x[1])+
                             (TOV_Position_y[0]-TOV_Position_y[1])*
                             (TOV_Position_y[0]-TOV_Position_y[1])+
                             (TOV_Position_z[0]-TOV_Position_z[1])*
                             (TOV_Position_z[0]-TOV_Position_z[1]));
  /* Now we need to calculate the coordinate distance from that */
  CCTK_REAL coord_dist = calc_coord_dist(prop_dist,
                                         TOV_m_1d[TOV_Num_Radial-1],
                                         TOV_RProp_Surface[0]);
  /* Now steer the parameters before the stars are interpolated onto
   * the Cactus grid using those very same variables */
  /* We do that by scaling every variable by the same factor */
  CCTK_REAL factor = coord_dist / prop_dist;
  char tmp[1024];
  snprintf(tmp, 1023, "%g", coord_dist/2);
  if (CCTK_ParameterSet("par_b", "TwoPunctures", tmp))
    CCTK_WARN(0, "Could not set par_b");
  snprintf(tmp, 1023, "%g", TOV_Position_x[0] * factor);
  if (CCTK_ParameterSet("TOV_Position_x[0]", "TOVSolver", tmp))
    CCTK_WARN(0, "Could not set x[0]");
  snprintf(tmp, 1023, "%g", TOV_Position_x[1] * factor);
  if (CCTK_ParameterSet("TOV_Position_x[1]", "TOVSolver", tmp))
    CCTK_WARN(0, "Could not set x[1]");

  // printf("proper distance: %g, coordinate distance: %g\n", prop_dist, coord_dist);
}

CCTK_REAL getrho(double p)                                                                 //adjusted
{
// make it better
 int j;
 if (p<7e-7) return 0.0;
 for(j=0;j<neos;j++)
 {
  if(p<myp[j]) break;
 }
 if(j==neos) return myrho[neos-1];
 if(j==0) return 0.0;
 return (myrho[j]-myrho[j-1])/(myp[j]-myp[j-1])*(p-myp[j-1])+myrho[j-1];
} 

CCTK_REAL getmu(double p)
{
// make it better
 int j;
 if (p<1e-10) return 0.0;
 for(j=0;j<neos;j++)
 {
  if(p<myp[j]) break;
 }
 if(j==neos) return mymu[neos-1];
 if(j==0) return 0.0;
 return (mymu[j]-mymu[j-1])/(myp[j]-myp[j-1])*(p-myp[j-1])+mymu[j-1];
} 

CCTK_REAL getpress(double rho)
{
// make it better
 int j;
 if (rho<5e-7) return 0.0;
 for(j=0;j<neos;j++)
 {
  if(rho<myrho[j]) break;
 }
 if(j==neos) return myp[neos-1];
 if(j==0) return 0.0;
 return (myp[j]-myp[j-1])/(myrho[j]-myrho[j-1])*(rho-myrho[j-1])+myp[j-1];
} 


#include "external.inc"

