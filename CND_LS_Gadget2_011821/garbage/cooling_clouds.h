///LAST MODIFIED : APRIL 25, 2016
//// VAZQUEZ COOLING ////
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>


#include "../allvars.h"

#define heat 2.00e-24 //// erg/s // galactic center heating
#define MH2 (2.0*PROTONMASS)

double cooling_heating(double T) /// input physical Kelvins
{ return ( 1e7 * exp( -114800.0 / (T + 1000.0)  ) + 0.014*sqrt(T)*exp(-92.0/T) );// OUTPUT IN CM^3
}

double cooling(double T) //input phys Kelvin
{ return cooling_heating(T) * heat ; // OUTPUT IN (ERGS/S) CM^3
}

struct cool_params {double n; double du_dt; } ; //input: n-cm^-3     du_dt- erg s^-1 g^-1

double result_of_cool(double T, void *p)
{ struct cool_params * params = (struct cool_params *)p;
    double n = (params->n);
    double du_dt = (params->du_dt);
    double cooling_solve = (heat + (du_dt * MH2) ) / n ;
    double f = cooling(T) - cooling_solve;
    return f;
}

double T_eq(double n, double du_dt)
{
    double T_lowest = 70.0; /// Kelvins
    const gsl_root_fsolver_type *T;  double T_lo = 0.0;  double T_hi = 1e6;
    int status;
    int iter = 0, max_iter = 200;
    gsl_root_fsolver *s;
    gsl_function F;
    double T_result;
    struct cool_params value = {n, du_dt};
    F.function = &result_of_cool;
    F.params = &value;
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, T_lo, T_hi);
    do
    { iter++;
        status = gsl_root_fsolver_iterate (s);
        T_result = gsl_root_fsolver_root (s);
        T_lo = gsl_root_fsolver_x_lower (s);
        T_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (T_lo, T_hi, 1e-3, 1e-6);
        if (status == GSL_SUCCESS)
        {       if (T_result < T_lowest) 
                    {return T_lowest;}
                else
                   { return T_result; }// OUTPUT KELVINS
        }
        
    } while (status == GSL_CONTINUE && iter < max_iter);
      gsl_root_fsolver_free (s);
      return status;
}

/*
double cool_for_root(double T, void *p)
{ struct cool_params * params = (struct cool_params *)p;
    double n = (params->n);
    double du_dt = (params->du_dt);
    return n - (1.0/cooling_heating(T)) - (du_dt*MH2 / cooling(T)  ) ;
}

double T_eq(double n, double du_dt)
{
    double T_lowest = 70.0; /// Kelvins
    const gsl_root_fsolver_type *T;  double T_lo = 0.0;  double T_hi = 1e6;
    int status;
    int iter = 0, max_iter = 200;
    gsl_root_fsolver *s;
    gsl_function F;
    double T_result;
    struct cool_params value = {n, du_dt};
    F.function = &cool_for_root;
    F.params = &value;
    T = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, T_lo, T_hi);
    do
    { iter++;
        status = gsl_root_fsolver_iterate (s);
        T_result = gsl_root_fsolver_root (s);
        T_lo = gsl_root_fsolver_x_lower (s);
        T_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (T_lo, T_hi, 1e-3, 1e-6);
        if (status == GSL_SUCCESS)
        {       if (T_result < T_lowest) 
                    {return T_lowest;}
                else
                   { return T_result; }// OUTPUT KELVINS
        }
        
    } while (status == GSL_CONTINUE && iter < max_iter);
     if (status != GSL_SUCCESS || status != GSL_CONTINUE )
        {  return T_lowest;
        }
    gsl_root_fsolver_free (s);
  //  return status;
}
*/
        //SphP[i].Entropy  SphP[i].Density   (All.UnitLength_in_cm)^2
double A_to_e(double A_c, double rho_c, double conversion_energy) //input A and rho in code units
{   double en = A_c * pow(rho_c, GAMMA_MINUS1) / GAMMA_MINUS1 ;
    return (en * conversion_energy) ;  // OUTPUT IN ergs/g
}
                 //SphP[i].DtEntropy  SphP[i].Density  (All.UnitLength_in_cm)^2  All.UnitTime_in_s
double dtA_to_dte(double dtA, double rho, double conversion_energy, double conversion_time) 
{   double du_dt =  dtA * pow(rho, GAMMA_MINUS1) / GAMMA_MINUS1 ;  ///code untis
    return (du_dt * conversion_energy/conversion_time);  //OUTPUT erg s^-1 g^-1
}

double dte_to_dtA(double dte, double rho, double conversion_energy, double conversion_time) 
{   double dA_dt = GAMMA_MINUS1 * (dte*conversion_time/conversion_energy) / pow(rho, GAMMA_MINUS1) ;  ///code untis
    return dA_dt; ///OUTPUT code units
}

            //SphP[i].Density  SphP[i].DtEntropy  (All.UnitLength_in_cm)^2  All.UnitTime_in_s
double e_eq(double rho, double dtA, double conversion_energy, double conversion_time, double conversion_density)
{   double n = rho * conversion_density / MH2 ; // in cm^-3
    double du_dt = dtA_to_dte( dtA, rho, conversion_energy, conversion_time ) ; //erg s^-1 g^-1
    double T_equi = T_eq(n, du_dt) ; // Kelvin
    double equilibrium_energy = ( ( (BOLTZMANN * T_equi) / MH2 ) * (1.0/GAMMA_MINUS1) ) ; //  erg/g
    if (equilibrium_energy == 0)
    { printf( "WARNING: EQUILIBRIUM ENERGY = 0") ; }
    return equilibrium_energy; //OUTPUT erg/g
}

              //erg/g   SphP[i].Density    (All.UnitLength_in_cm)^2
double e_to_A(double en, double rho, double conversion_energy) 
{   double A = GAMMA_MINUS1 * (en/conversion_energy) / pow(rho, GAMMA_MINUS1) ; //code units
    return A;  // OUTPUT IN CODE UNITS
}
             /// input erg/g
double e_to_T(double en) 
{   return (en * GAMMA_MINUS1 * MH2 / BOLTZMANN ); //OUTPUT IN Kelvin
}
                //input Kelvin
double T_to_e(double T)
{ double energy = ( BOLTZMANN * T / MH2) * (1.0 / GAMMA_MINUS1) ;
    return energy; // OUTPUT erg / g
}

                //SphP[i].Entropy   SphP[i].Density   SphP[i].DtEntropy  (All.UnitLength_in_cm)^2  All.UnitTime_in_s
double cool_time(double A, double rho, double dtA, double conversion_energy, double conversion_time, double conversion_density)
{   double en = A_to_e( A,  rho, conversion_energy) ; // e_now in ergs / g
    double Temp = e_to_T(en) ; // T_now in Kelvin
    double e_equi = e_eq( rho,  dtA,  conversion_energy,  conversion_time, conversion_density) ; //equilibrium energy in erg/g
    double n = rho * conversion_density / MH2 ; //in cm^-3
    double du_dt = dtA_to_dte(dtA, rho,  conversion_energy,  conversion_time) ;  ///erg s^-1 g^-1
    double time_cool = ( (en - e_equi) / ( (n*cooling(Temp)/MH2 ) - du_dt - (heat/MH2) ) )  ;
    return (fabs(time_cool) / conversion_time ); //OUTPUT code units
    
}
         //SphP[i].Entropy   SphP[i].Density   SphP[i].DtEntropy  (All.UnitLength_in_cm)^2  All.UnitTime_in_s  dt_entr
double new_A(double A, double rho, double dtA, double conversion_energy, double conversion_time,double conversion_density, double dt)
{   double n_check =  rho * conversion_density / MH2 ; double A_new = 0;
    if (n_check < 100.0 ) { A_new = A + (dtA*dt) ; }
    else 
    {
    double e_equi = e_eq( rho,  dtA,  conversion_energy,  conversion_time, conversion_density) ; //equilibrium energy in erg/g
    double A_equi = e_to_A( e_equi, rho, conversion_energy) ; // in code units
    double t_cool = cool_time(A,  rho,  dtA,  conversion_energy, conversion_time, conversion_density) ; // in code units
    A_new = A_equi + (A - A_equi)*exp( - dt / t_cool) ;
    }
    return A_new; 
}