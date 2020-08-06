#ifndef ODE_SOLVER_HPP_
#define ODE_SOLVER_HPP_

#include "GraphDyn_util.hpp"
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
namespace transient_analysis {
class ODE_solver {
 public:
  ODE_solver() {}
  //void operator()(const vector_type& x, vector_type& dxdt, real__t t);
  void operator()(const d_vector_type &x, d_vector_type &dxdt, const value_type t);
  
  real__t get_Pmech() {return Pmech;};
  real__t get_Efd() {return Efd;};
  real__t get_VS() {return VS;};
  real__t get_mu() {return mu;};
  real__t get_rectified_regulator(real__t src);

  real__t Vx, Vy;
  real__t Id, Iq;
  //double d_Vx;
  //cudaMemcpy(d_Vx, Vx, 1*sizeof(double, cudaMemcpyHostToDevice));
  //cudaMalloc();
  //thrust::device_ptr<double> d_Vx  = thrust::device_malloc<double>(0);
  //thrust::device_ptr<double> d_parameters = thrust::device_malloc()
  GENERATOR_DATA_PACK parameters;
  GENERATOR_DATA_PACK_D d_parameters;
  

 private:
  /** helper functions */
  real__t apply_limiting(real__t val, const real__t val_min, const real__t val_max);
  real__t apply_dead_band(const real__t val, const real__t tol);
  void integrate_block(const vector_type& x, vector_type& dxdt, int idx_src, int idx_dst,
                       real__t a, real__t b, real__t c, real__t d);
  void integrate_block(vector_type& dxdt, real__t val_src, real__t val_dst, int idx_dst,
                       real__t a, real__t c, real__t d);
  void integrate_block(const vector_type& x, vector_type& dxdt, real__t val_src, real__t div_src,
                       real__t val_dst, int idx_dst, real__t a, real__t b, real__t c, real__t d);
  real__t process_PID_block(const vector_type& x, vector_type& dxdt, real__t val_src, real__t div_src,
                            int idx_pid, PID_DATA& pid);

  /** Generator models */
  void process_EPRI_GEN_TYPE_I_D(const d_vector_type& x, d_vector_type& dxdt, value_type TJ, value_type D);
  void process_EPRI_GEN_TYPE_ZERO(const vector_type& x, vector_type& dxdt, EPRI_GEN_DATA& gen);
  void process_EPRI_GEN_TYPE_I(const vector_type& x, vector_type& dxdt, EPRI_GEN_DATA& gen);
  void process_EPRI_GEN_TYPE_II(const vector_type& x, vector_type& dxdt, EPRI_GEN_DATA& gen);
  void process_EPRI_GEN_TYPE_III(const vector_type& x, vector_type& dxdt, EPRI_GEN_DATA& gen);
  void process_EPRI_GEN_TYPE_VI(const vector_type& x, vector_type& dxdt, EPRI_GEN_DATA& gen);
  void process_GENROU(const vector_type& x, vector_type& dxdt, GENROU_IEEE_DATA& gen);

  /** EPRI Governor models */
  real__t electro_hydraulic_servo(const vector_type& x, vector_type& dxdt, EHS_DATA& ehs,
                                  int idx_src, int idx_pid, int idx_feedback, int idx_dst);
  real__t process_steam_machine_third_order(const vector_type& x, vector_type& dxdt, const STEAM_DATA& steam,
                                            real__t val_src, const int idx_hp, const int idx_ip, const int idx_lp);
  real__t process_EPRI_GOV_TYPE_I(const vector_type& x, vector_type& dxdt, EPRI_GOV_I_DATA& gov);
  real__t process_EPRI_GOV_TYPE_II(const vector_type& x, vector_type& dxdt, EPRI_GOV_II_DATA& gov);
  real__t process_EPRI_GOV_TYPE_III(const vector_type& x, vector_type& dxdt, EPRI_GOV_III_DATA& gov);
  real__t process_EPRI_GOV_TYPE_IV(const vector_type& x, vector_type& dxdt, EPRI_GOV_IV_DATA& gov);
  real__t process_EPRI_GOV_TYPE_V(const vector_type& x, vector_type& dxdt, EPRI_GOV_V_DATA& gov);
//  real__t process_EPRI_GOV_TYPE_VI(const vector_type& x, vector_type& dxdt, EPRI_GOV_VI_DATA& gov); // to do
  real__t process_EPRI_GOV_TYPE_VII(const vector_type& x, vector_type& dxdt, EPRI_GOV_VII_DATA& gov);
  real__t process_EPRI_GOV_TYPE_VIII(const vector_type& x, vector_type& dxdt, EPRI_GOV_VIII_DATA& gov);
  real__t process_EPRI_GOV_TYPE_IX(const vector_type& x, vector_type& dxdt, EPRI_GOV_IX_DATA& gov);

  /** IEEE Governor BPAGG Model, used for PowerWorld models */
  real__t process_GOV_BPAGG(const vector_type& x, vector_type& dxdt, GOV_BPAGG_DATA& gov);

  /** EPRI Excitor models */
  real__t get_exc_input(real__t Xc);
  void print_for_debug_EXC_III_X(const vector_type& x, vector_type& dxdt, real__t V_in, real__t val_Vdiff,
                                 real__t div_Vdiff, real__t Vfe_val, real__t Efd_val);
  real__t process_EPRI_EXC_TYPE_I(const vector_type& x, vector_type& dxdt, EPRI_EXC_I_DATA& exc);
  real__t process_EPRI_EXC_TYPE_II(const vector_type& x, vector_type& dxdt, EPRI_EXC_II_DATA& exc);
  real__t process_EPRI_EXC_TYPE_III(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc);
  real__t process_EPRI_EXC_TYPE_IV(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc);
  real__t process_EPRI_EXC_TYPE_V(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc);
  real__t process_EPRI_EXC_TYPE_VI(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc);
  real__t process_EPRI_EXC_TYPE_VII(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc);
  real__t process_EPRI_EXC_TYPE_VIII(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc);
  real__t process_EPRI_EXC_TYPE_IX(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc);
  real__t process_EPRI_EXC_TYPE_X(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc);
  real__t process_EPRI_EXC_TYPE_XI(const vector_type& x, vector_type& dxdt, EPRI_EXC_XI_TO_XII_DATA& exc);
  real__t process_EPRI_EXC_TYPE_XII(const vector_type& x, vector_type& dxdt, EPRI_EXC_XI_TO_XII_DATA& exc);

  /** IEEE Excitor Model, used for PowerWorld models */
  real__t process_EXC_IEEE_I(const vector_type& x, vector_type& dxdt, EXC_IEEE_I_DATA& exc);

  /** EPRI PSS models */
  real__t process_EPRI_PSS_TYPE_I(const vector_type& x, vector_type& dxdt, EPRI_PSS_I_DATA& pss);
  real__t process_EPRI_PSS_TYPE_II(const vector_type& x, vector_type& dxdt, EPRI_PSS_II_DATA& pss);
  real__t process_EPRI_PSS_TYPE_IV(const vector_type& x, vector_type& dxdt, EPRI_PSS_IV_VI_DATA& pss);
  real__t process_EPRI_PSS_TYPE_V(const vector_type& x, vector_type& dxdt, EPRI_PSS_V_DATA& pss);
  real__t process_EPRI_PSS_TYPE_VI(const vector_type& x, vector_type& dxdt, EPRI_PSS_IV_VI_DATA& pss);
  real__t process_EPRI_PSS_TYPE_VIII(const vector_type& x, vector_type& dxdt, EPRI_PSS_VIII_DATA& pss);

  void apply_perturbation(real__t t, EPRI_GEN_DATA& gen);
  void update_generator_current(const d_vector_type& x, EPRI_GEN_DATA& gen);
  void setup(const d_vector_type& x, d_vector_type& dxdt, EPRI_GEN_DATA& gen);
  //void setup(const d_vector_type& x, d_vector_type& dxdt);
  void print_dxdt(const vector_type& dxdt);
  
  real__t omega_ref, Vm_ref, Pe_ref, Pm_ref, freq_ref;
  real__t omega_ref0, Vm_ref0, Pe_ref0, Pm_ref0, freq_ref0;

 /* these are variables for intermediate computation */
  real__t Ifd, dIfd_dt;
  real__t Vd, Vq;
  real__t Vm, Va;
  real__t Efd0, mu0;
  
  real__t Telec;
  real__t mu = 0.;
  real__t Pmech = 0.;
  real__t Efd = 0.;
  real__t VS = 0.;
};

}  // namespace transient_analysis

#endif
