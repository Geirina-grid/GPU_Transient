#include "ode_solver.hpp"
#include <trans_cuda.hpp>
#include <transient.hpp>
#include <map>
#include <iostream>
#include <cmath>
#include <utility>


#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>

#include <boost/numeric/odeint.hpp>

#include <boost/numeric/odeint/external/thrust/thrust.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>


namespace transient_analysis {


real__t ODE_solver::
apply_limiting(real__t val, const real__t val_min, const real__t val_max) {
   //cout << "In apply_limiting .x="  << endl;
  return (val < val_min) ? val_min : ((val > val_max) ? val_max : val);
}

real__t ODE_solver::apply_dead_band(const real__t val, const real__t tol) {
   //cout << "In apply_dead_band .x="  << endl;
  return (abs(val) <= tol) ? 0 : val;
}

/**  integrate_block integrates one step of the following ODE:
 *    (a + bs) / (c + ds) x = y
*/
void ODE_solver::
integrate_block(const vector_type& x, vector_type& dxdt, int idx_src, int idx_dst,
                real__t a, real__t b, real__t c, real__t d) {
   //cout << "In integrate_block 1.x="  << endl;
  dxdt[idx_dst] = (abs(d) < EPS) ? 0 : (a * x[idx_src] + b * dxdt[idx_src] - c * x[idx_dst]) / d;
}

/**  integrate_block integrates one step of the following ODE:
 *    a / (c + ds) x = y
*/
void ODE_solver::
integrate_block(vector_type& dxdt, real__t val_src, real__t val_dst, int idx_dst,
                real__t a, real__t c, real__t d) {
   //cout << "In integrate_block 2.x="  << endl;
  dxdt[idx_dst] = (abs(d) < EPS) ? 0 : (a * val_src - c * val_dst) / d;
}

/**  integrate_block integrates one step of the following ODE:
 *    (a + bs) / (c + ds) x = y,
 * and directly takes the values of x and dx_dt as inputs.
*/

void ODE_solver::
integrate_block(const vector_type& x, vector_type& dxdt, real__t val_src, real__t div_src,
                real__t val_dst, int idx_dst, real__t a, real__t b, real__t c, real__t d) {
  //cout << "In integrate_block 3.x="  << endl;
  dxdt[idx_dst] = (abs(d) < EPS) ? 0 : (a * val_src + b * div_src - c * val_dst) / d;
}

/**  process_PID_block does two things:
 * 1. calculate the output value
 * 2. process the ODE:    s y = Ki * x
*/
real__t ODE_solver::
process_PID_block(const vector_type& x, vector_type& dxdt, real__t val_src, real__t div_src,
                  int idx_pid, PID_DATA &pid) {
  //cout << "In process_PID_block x="  << endl;
  real__t val1 = pid.Kp * val_src;
  real__t val2 = pid.Kd * div_src;
  dxdt[idx_pid] = pid.Ki * val_src;
  real__t val3 = apply_limiting(x[idx_pid], 10. * pid.I_Min, 10. * pid.I_Max);
  real__t val_total = val1 + val2 + val3;
  return apply_limiting(val_total, 10. * pid.PID_Min, 10. * pid.PID_Max);
}

/** a helper function for debugging */
void ODE_solver::print_dxdt(const vector_type& dxdt) {
  //cout << "In print_dx_dt" << endl;
  TextTable t('-', '|', '+');

  t.addRow(vector<string>{"omega", "delta", "Eqp", "Edp", "Eqpp", "Edpp", "delta_mu"});

  t.add(to_string(dxdt[omega_idx]));
  t.add(to_string(dxdt[delta_idx]));
  t.add(to_string(dxdt[Eqp_idx]));
  t.add(to_string(dxdt[Edp_idx]));
  t.add(to_string(dxdt[Eqpp_idx]));
  t.add(to_string(dxdt[Edpp_idx]));
  t.add(to_string(dxdt[delta_mu_idx]));
  t.endOfRow();

  std::cout << t;
}

void ODE_solver::update_generator_current(const d_vector_type& x, EPRI_GEN_DATA& gen) {
  //if (parameters.GEN_type == 3 || parameters.GEN_type == 6) {
  //  real__t denom = gen.Ra * gen.Ra + gen.Xdpp * gen.Xqpp;
  //  //assert(denom > EPS);
  //  Id = (+gen.Ra * (x[Edpp_idx] - Vd) + gen.Xqpp * (x[Eqpp_idx] - Vq)) / denom;
  //  Iq = (-gen.Xdpp * (x[Edpp_idx] - Vd) + gen.Ra * (x[Eqpp_idx] - Vq)) / denom;
  //} else {
    real__t denom = gen.Ra * gen.Ra + gen.Xdp * gen.Xqp;
    //assert(abs(denom) > EPS);
    Id = (+gen.Ra * (x[Edp_idx] - Vd) + gen.Xqp * (x[Eqp_idx] - Vq)) / denom;
    Iq = (-gen.Xdp * (x[Edp_idx] - Vd) + gen.Ra * (x[Eqp_idx] - Vq)) / denom;
  //}
}

//void ODE_solver::setup(const d_vector_type& x, d_vector_type& dxdt) {
void ODE_solver::setup(const d_vector_type& x, d_vector_type& dxdt, EPRI_GEN_DATA& gen) {
  
  //Vm_refomega_ref = parameters.omega_ref;
  //Vm_ref    = parameters.Vt_ref;
  //Pe_ref    = parameters.Pe_ref;
  //Pm_ref    = parameters.Pe_ref;
  ///freq_ref  = parameters.freq_ref;
  //Efd0      = parameters.Efd0;
  //mu0       = parameters.mu0;

  printf("omega_idx=%d\n", omega_idx);
  printf("delta_idx=%d\n", delta_idx);
  printf("Eqp_idx=%d\n", Eqp_idx);
  printf("gen.Ra=%d\n", gen.Ra);
  std::cout<<"x[omega_idx]" << "="<< x[omega_idx] << std::endl;
  std::cout<<"x[delta_idx]" << "="<< x[delta_idx] << std::endl;
  std::cout<<"x[Eqp_idx]" << "="<< x[Eqp_idx] << std::endl;
  std::cout<<"gen.a" << "="<< gen.a << std::endl;
  std::cout<<"gen.b" << "="<< gen.b << std::endl;
  std::cout<<"gen.n" << "="<< gen.n << std::endl;
  std::cout<<"gen.Ra" << "="<< gen.Ra << std::endl;

  printf("\n");
  
  Vm = sqrt(Vx *Vx + Vy * Vy);
  Va = atan2(Vy, Vx);
  Vd = Vm * sin(x[delta_idx] - Va);
  Vq = Vm * cos(x[delta_idx] - Va);
  
  Ifd     = gen.a * x[Eqp_idx] + gen.b * pow(x[Eqp_idx], gen.n);
  dIfd_dt = gen.a * dxdt[Eqp_idx] + gen.b * gen.n * pow(x[Eqp_idx], gen.n - 1) * dxdt[Eqp_idx];
  
  update_generator_current(x, gen);
  
  real__t Psi_q = -(gen.Ra * Id + Vd) / x[omega_idx];
  real__t Psi_d = +(gen.Ra * Iq + Vq) / x[omega_idx];

  Telec = Psi_d * Iq - Psi_q * Id;
  
  for(value_type v: dxdt)  {v = 0;}
  //printf("current thread=%d", threadIdx.x);
}
//void ODE_solver::
//operator()(const vector_type& x, vector_type& dxdt, real__t t) {
//  printf(" parameters.gen=%s",  parameters.gen);
//  value_type &dx = thrust::get< 0 >(x);
//  setup(x, dxdt, parameters.gen);

//}
//class ODE_solver{

//public:
//  void operator() ( const d_vector_type &x , d_vector_type &dxdt , const value_type dt )
//  {
//     printf("Hello world");
//  }
//}

void ODE_solver::
operator()(const d_vector_type &x, d_vector_type &dxdt, const value_type t) {
  //printf("x[omega_idx] in solver=%f\n", x[omega_idx]);
  //printf("dxdt in solver=%f\n", dxdt);
  //printf("t in solver=%f\n", t);
  setup(x, dxdt, d_parameters.gen);
  apply_perturbation(t,d_parameters.gen);
  update_generator_current(x, d_parameters.gen);
  //setup(x, dxdt, d_parameters.gen);
  //printf("parameters.GEN_type=%d\n", parameters.GEN_type);
  
#if DEBUG
  printf("\n\nrunning ODE solver with Gen type %d, GOV type %d, AVR type %d, PSS type %d...\n",
         parameters.GEN_type, parameters.GOV_type, parameters.EXC_type, parameters.PSS_type);
#endif

  /*
  switch (parameters.PSS_type) {
    case 0: VS = 0; break;
    case 1: VS = process_EPRI_PSS_TYPE_I(x, dxdt, parameters.pss_1);    break;
    case 2: VS = process_EPRI_PSS_TYPE_II(x, dxdt, parameters.pss_2);   break;
    case 4: VS = process_EPRI_PSS_TYPE_IV(x, dxdt, parameters.pss_4_6); break;
    case 5: VS = process_EPRI_PSS_TYPE_V(x, dxdt, parameters.pss_5);    break;
    case 8: VS = process_EPRI_PSS_TYPE_VIII(x, dxdt, parameters.pss_8); break;
    default: {std::cerr << "Error: unsupported PSS type...\n"; std::terminate(); break;}
  }
  
  switch (parameters.EXC_type) {
    case 0:  Efd = Efd0; break;
    case 1:  Efd = process_EPRI_EXC_TYPE_I(x, dxdt, parameters.exc_1);       break;
    case 2:  Efd = process_EPRI_EXC_TYPE_II(x, dxdt, parameters.exc_2);      break;
    case 3:  Efd = process_EPRI_EXC_TYPE_III(x, dxdt, parameters.exc_3_10);  break;
    case 4:  Efd = process_EPRI_EXC_TYPE_IV(x, dxdt, parameters.exc_3_10);   break;
    case 5:  Efd = process_EPRI_EXC_TYPE_V(x, dxdt, parameters.exc_3_10);    break;
    case 6:  Efd = process_EPRI_EXC_TYPE_VI(x, dxdt, parameters.exc_3_10);   break;
    case 7:  Efd = process_EPRI_EXC_TYPE_VII(x, dxdt, parameters.exc_3_10);  break;
    case 8:  Efd = process_EPRI_EXC_TYPE_VIII(x, dxdt, parameters.exc_3_10); break;
    case 9:  Efd = process_EPRI_EXC_TYPE_IX(x, dxdt, parameters.exc_3_10);   break;
    case 10: Efd = process_EPRI_EXC_TYPE_X(x, dxdt, parameters.exc_3_10);    break;
    case 11: Efd = process_EPRI_EXC_TYPE_XI(x, dxdt, parameters.exc_11_12);  break;
    case 12: Efd = process_EPRI_EXC_TYPE_XII(x, dxdt, parameters.exc_11_12); break;
    default: {std::cerr << "Error: unsupported excitor (AVR) type...\n"; std::terminate(); break;}
  }
  
  switch (parameters.GOV_type) {
    case 0: Pmech = Pe_ref; break;
    case 1: Pmech = process_EPRI_GOV_TYPE_I(x, dxdt, parameters.gov_1);    break;
    case 3: Pmech = process_EPRI_GOV_TYPE_III(x, dxdt, parameters.gov_3);  break;
    case 4: Pmech = process_EPRI_GOV_TYPE_IV(x, dxdt, parameters.gov_4);   break;
    case 5: Pmech = process_EPRI_GOV_TYPE_V(x, dxdt, parameters.gov_5);    break;
//    case 6: Pmech = process_EPRI_GOV_TYPE_VI(x, dxdt, parameters.gov_6); break;
    case 7: Pmech = process_EPRI_GOV_TYPE_VII(x, dxdt, parameters.gov_7);  break;
    case 8: Pmech = process_EPRI_GOV_TYPE_VIII(x, dxdt, parameters.gov_8); break;
    case 9: Pmech = process_EPRI_GOV_TYPE_IX(x, dxdt, parameters.gov_9);   break;
    default: {std::cerr << "Error: unsupported governor (GOV) type...\n"; std::terminate(); break;}
  }
  */

  switch (parameters.GEN_type) {
    //case 0:
    //case 1: process_EPRI_GEN_TYPE_I(x, dxdt, parameters.gen);   break;
    //case 3: process_EPRI_GEN_TYPE_III(x, dxdt, parameters.gen); break;
    //case 6: process_EPRI_GEN_TYPE_VI(x, dxdt, parameters.gen);  break;
    case 6: process_EPRI_GEN_TYPE_I_D(x, dxdt, parameters.gen.TJ, parameters.gen.D);  break;
    //default: {std::cerr << "Error: unsupported generator (GEN) type...\n"; std::terminate(); break;}
  }


//  print_dxdt(dxdt);
}

}  // namespace transient_analysis
