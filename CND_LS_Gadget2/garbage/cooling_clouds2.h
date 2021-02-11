#include "cooling_clouds1.h"

struct cool_params2
{
    double A_c;
    double dtA_c;
    double dt_c;
    double rho_c;
    double unit_mass;
    double unit_time;
    double unit_ergs;
    double unit_cm;
} ;

double result_of_cool2(double T_new, void *p)
{ struct cool_params2 * params2 = (struct cool_params2 *)p;
    double A_c = (params2->A_c);
    double dtA_c = (params2->dtA_c);
    double dt_c = (params2->dt_c) ;
    double rho_c = (params2->rho_c);
    double unit_mass = (params2->unit_mass);
    double unit_time = (params2->unit_time);
    double unit_ergs = (params2->unit_ergs);
    double unit_cm   = (params2->unit_cm)  ;
    double boltz_c = BOLTZMANN / unit_ergs ;
    double mu_c = MH2 / unit_mass ;
    double n_c = rho_c / mu_c;
    double conversion = GAMMA_MINUS1/pow(rho_c, GAMMA) ;
    double newTcoeff = (boltz_c/mu_c) *(1.0/pow(rho_c,GAMMA_MINUS1) ) ;
    
    double f = A_c + (dtA_c*dt_c) + ((conversion*n_c*heating_code(unit_time, unit_ergs))*dt_c)  - (dt_c*conversion*n_c*n_c*lambda_code(T_new, unit_time,  unit_ergs,  unit_cm) ) - newTcoeff*T_new;
    
    return f;
}

double find_T_new(double A_c, double dtA_c, double dt_c, double rho_c, double unit_mass,
                  double unit_time, double unit_ergs, double unit_cm)
{
    const gsl_root_fsolver_type *T;  double T_lo = 0.0;  double T_hi = 1e6;
    int status;
    int iter = 0, max_iter = 200;
    gsl_root_fsolver *s;
    gsl_function F;
    double T_result;
    struct cool_params2 value2 = {A_c, dtA_c, dt_c, rho_c, unit_mass, unit_time,
                                    unit_ergs, unit_cm};
    F.function = &result_of_cool2;
    F.params = &value2;
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
        {    
            return T_result;     //// OUTPUT KELVINS
        }
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);
    return status;

}

double dtA_coolheat_c(double A_c, double dtA_c, double dt_c, double rho_c, double unit_mass,
                      double unit_time, double unit_ergs, double unit_cm, double unit_energy)
{ 
    double dtA_coolheat_code;
    double dta_cooling_c;
    double final_result;
    
    double conversion = GAMMA_MINUS1/pow(rho_c, GAMMA) ;
    double mu_c = MH2 / unit_mass ;
    double n_c = rho_c / mu_c;
    double dtA_heat_c = conversion*n_c*heating_code(unit_time, unit_ergs);
    
   // double en_now = A_to_e( A_c,  rho_c,  unit_energy) ;
    //double t_now = e_to_T(en_now);
    double T_found = find_T_new(A_c,  dtA_c,  dt_c,  rho_c,  unit_mass,  unit_time, unit_ergs, unit_cm);
    if (T_found < 80.0)
    {
        T_found = 80.0;
    }
    
    dta_cooling_c = conversion*n_c*n_c*lambda_code(T_found, unit_time,  unit_ergs,  unit_cm);    
    dtA_coolheat_code = (dtA_heat_c - dta_cooling_c)*dt_c ;
    
    final_result = (dtA_c * dt_c ) + dtA_coolheat_code;
    return final_result;
}


double gas_cooling_entropy(double A_c, double dtA_c, double dt_c, double rho_c, double unit_mass,
                       double unit_time, double unit_ergs, double unit_cm, double unit_energy, double unit_density)
{
    double dtA_coolheat_code;
    double dta_cooling_c;
    double new_A_c;
    
    double conversion = GAMMA_MINUS1/pow(rho_c, GAMMA) ;
    double mu_c = MH2 / unit_mass ;
    double n_c = rho_c / mu_c;
    double dtA_heat_c = conversion*n_c*heating_code(unit_time, unit_ergs);
    
    double T_found;
    
    double target_density = 17.9 * MH2 / unit_density;
    
    if (rho_c < target_density)
    {
        T_found = find_T_new(A_c,  dtA_c,  dt_c,  rho_c,  unit_mass,  unit_time, unit_ergs, unit_cm);
        dta_cooling_c = conversion*n_c*n_c*lambda_code(T_found, unit_time,  unit_ergs,  unit_cm);
        dtA_coolheat_code = (dtA_heat_c - dta_cooling_c)*dt_c ;
        
        new_A_c = A_c + ( (dtA_c * dt_c ) + dtA_coolheat_code );

    }
    
    else
    {
        
        T_found = find_T_new(A_c,  dtA_c,  dt_c,  rho_c,  unit_mass,  unit_time, unit_ergs, unit_cm);
        if (T_found < 80.0)
        {
            double myenergy = T_to_e(80.0);
            new_A_c = e_to_A(myenergy, rho_c, unit_energy) ;
        }
        else
        {
            dta_cooling_c = conversion*n_c*n_c*lambda_code(T_found, unit_time,  unit_ergs,  unit_cm);
            dtA_coolheat_code = (dtA_heat_c - dta_cooling_c)*dt_c ;
            
            new_A_c = A_c + ( (dtA_c * dt_c ) + dtA_coolheat_code );
        }

    
    }
    
    return new_A_c;
}



