#include "cooling_clouds.h"
#include <gsl/gsl_spline.h>

double lambda_code(double temp, double unit_time, double unit_ergs, double unit_cm)
{   double unit_volume = pow(unit_cm,3.0);
    double lambda_c = cooling(temp) * (unit_time/ (unit_ergs*unit_volume));
    return lambda_c;
}

double heating_code(double unit_time, double unit_ergs)
{
    double h_c = heat * unit_time/ unit_ergs ;
    return h_c;
}

double compute_L(double A_c, double rho_c, double unit_mass, double unit_energy, double unit_ergs, double unit_cm, double unit_time)
{
    double mu_c = MH2 / unit_mass;
    double n_c = rho_c / mu_c ;
    double energy_now_p = A_to_e( A_c,  rho_c,  unit_energy);
    double temp_now = e_to_T(energy_now_p);
    double lambda_c = lambda_code(temp_now, unit_time, unit_ergs,  unit_cm);
    double heat_c = heating_code(unit_time, unit_ergs);
    double L_c =   (n_c*heat_c) - (pow(n_c,2.0)*lambda_c); ///negative or positive?
    return L_c;
}

double cool_value(double dtA_c, double rho_c, double unit_mass, double unit_time, double unit_ergs)
{   double mu_c = MH2 / unit_mass;
    double n_c = rho_c / mu_c ;
    double cool_i_c = -( (dtA_c * pow(rho_c, GAMMA) / GAMMA_MINUS1) - (n_c*heating_code(unit_time, unit_ergs)) ) * (1.0/pow(n_c,2.0) ) ;
    //if (cool_i_c < 0)
      //  { printf("NEGATIVE VALUE OF EQUILIBRIUM COOLING\n") ;
   // printf("Cool value = %E\n", cool_i_c);
    return fabs(cool_i_c); //negative or positive?
    
}

