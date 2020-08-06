#ifndef ODE_SOLVER_HPP_
#define ODE_SOLVER_HPP_
#define GEN_SIZE 8

#include "GraphDyn_util.hpp"
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <thrust/device_vector.h>
#include <iostream>
namespace transient_analysis {
class ODE_solver {

 public:
  ODE_solver() {}
  //void operator()(const vector_type& x, vector_type& dxdt, real__t t);
  void operator()(const d_vector_type &x, d_vector_type &dxdt, const value_type t);
  
  real__t get_Pmech() {return Pmech[0];};
  real__t get_Efd() {return Efd[0];};
  real__t get_VS() {return VS[0];};
  real__t get_mu() {return mu[0];};
  real__t get_rectified_regulator(real__t src);

  //thrust::device_vector<int> d_vec(8);

  //int Vx = d_vec[0];
  //int Vy = d_vec[1];
  //int Id = d_vec[2];
  //int Iq = d_vec[3];

  //thrust::device_vector<double> Vm(8);
  //thrust::device_vector<double> Vx(8);
  //thrust::device_vector<double> Vy(8);
  //thrust::device_vector<double> Va(8);
  //thrust::device_vector<double> Vd(8);
  //thrust::device_vector<double> Vq(8);
  //thrust::device_vector<double> Ifd(8);

  //d_vector_type Vm(8, 0);
  //d_vector_type Vx(8, 0);
  //d_vector_type Vy(8, 0);
  //d_vector_type Va(8, 0);
  //d_vector_type Vd(8, 0);
  //d_vector_type Vq(8. 0);
  //d_vector_type Ifd(8, 0);

  //std::vector<thrust::host_vector<double>> Vm(8);
  //vector<double> Vm[8];
  //real__t Vm[8];
  //real__t Vx[8];
  //real__t Vy[8];
  //real__t Va[8];
  //real__t Vd[8];
  //real__t Vq[8];
  //real__tIfd[8];
  
  real__t Vx[GEN_SIZE], Vy[GEN_SIZE];
  real__t Id[GEN_SIZE], Iq[GEN_SIZE];

  real__t Vx_0, Vy_0;
  real__t Id_0, Iq_0;

  real__t Vx_1, Vy_1;
  real__t Id_1, Iq_1;

  real__t Vx_2, Vy_2;
  real__t Id_2, Iq_2;

  real__t Vx_3, Vy_3;
  real__t Id_3, Iq_3;

  real__t Vx_4, Vy_4;
  real__t Id_4, Iq_4;

  real__t Vx_5, Vy_5;
  real__t Id_5, Iq_5;

  real__t Vx_6, Vy_6;
  real__t Id_6, Iq_6;

  real__t Vx_7, Vy_7;
  real__t Id_7, Iq_7;
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
  
  real__t omega_ref[i], Vm_ref[i], Pe_ref, Pm_ref, freq_ref;
  real__t omega_ref0, Vm_ref0, Pe_ref0, Pm_ref0, freq_ref0;
  
  real__t omega_ref_0, Vm_ref_0, Pe_ref_0, Pm_ref_0, freq_ref_0;
  real__t omega_ref0_0, Vm_ref0_0, Pe_ref0_0, Pm_ref0_0, freq_ref0_0;

  real__t omega_ref_1, Vm_ref_1, Pe_ref_1, Pm_ref_1, freq_ref_1;
  real__t omega_ref0_1, Vm_ref0_1, Pe_ref0_1, Pm_ref0_1, freq_ref0_1;

  real__t omega_ref_2, Vm_ref_2, Pe_ref_2, Pm_ref_2, freq_ref_2;
  real__t omega_ref0_2, Vm_ref0_2, Pe_ref0_2, Pm_ref0_2, freq_ref0_2;

  real__t omega_ref_4, Vm_ref_4, Pe_ref_4, Pm_ref_4, freq_ref_4;
  real__t omega_ref0_4, Vm_ref0_4, Pe_ref0_4, Pm_ref0_4, freq_ref0_4;

  real__t omega_ref_3, Vm_ref_3, Pe_ref_3, Pm_ref_3, freq_ref_3;
  real__t omega_ref0_3, Vm_ref0_3, Pe_ref0_3, Pm_ref0_3, freq_ref0_3;

  real__t omega_ref_5, Vm_ref_5, Pe_ref_5, Pm_ref_5, freq_ref_5;
  real__t omega_ref0_5, Vm_ref0_5, Pe_ref0_5, Pm_ref0_5, freq_ref0_5;

  real__t omega_ref_6, Vm_ref_6, Pe_ref_6, Pm_ref_6, freq_ref_6;
  real__t omega_ref0_6, Vm_ref0_6, Pe_ref0_6, Pm_ref0_6, freq_ref0_6;

  real__t omega_ref_7, Vm_ref_7, Pe_ref_7, Pm_ref_7, freq_ref_7;
  real__t omega_ref0_7, Vm_ref0_7, Pe_ref0_7, Pm_ref0_7, freq_ref0_7;


 /* these are variables for intermediate computation */


  real__t Vm[GEN_SIZE];
  real__t Ifd[GEN_SIZE], dIfd_dt[GEN_SIZE];
  real__t dIfd_dt[GEN_SIZE];
  real__t Vd[GEN_SIZE], Vq[GEN_SIZE];
  real__t Vm[GEN_SIZE], Va[GEN_SIZE];
  real__t Efd0, mu0;

  real__t Telec[GEN_SIZE];
  //real__t mu[GEN_SIZE] = [0.];
  //real__t Pmech[GEN_SIZE] = [0.];
  //real__t Efd[GEN_SIZE] = [0.];
  //real__t VS[GEN_SIZE] = [0.];
  vector<real__t> mu(GEN_SIZE, 0.);
  vector<real__t> Pmech(GEN_SIZE, 0.);
  vector<real__t> Efd(GEN_SIZE, 0.);
  vector<real__t> VS(GEN_SIZE, 0.);


  real__t Ifd_0, dIfd_dt_0;
  real__t Vd_0, Vq_0;
  real__t Vm_0, Va_0;
  real__t Efd0_0, mu0_0;
  
  real__t Telec_0;
  real__t mu_0 = 0.;
  real__t Pmech_0 = 0.;
  real__t Efd_0 = 0.;
  real__t VS_0 = 0.;

  real__t Ifd_1, dIfd_dt_1;
  real__t Vd_1, Vq_1;
  real__t Vm_1, Va_1;
  real__t Efd0_1, mu0_1;

  real__t Telec_1;
  real__t mu_1 = 0.;
  real__t Pmech_1 = 0.;
  real__t Efd_1 = 0.;
  real__t VS_1 = 0.;


  real__t Ifd_2, dIfd_dt_2;
  real__t Vd_2, Vq_2;
  real__t Vm_2, Va_2;
  real__t Efd0_2, mu0_2;

  real__t Telec_2;
  real__t mu_2 = 0.;
  real__t Pmech_2 = 0.;
  real__t Efd_2 = 0.;
  real__t VS_2 = 0.;


  real__t Ifd_3, dIfd_dt_3;
  real__t Vd_3, Vq_3;
  real__t Vm_3, Va_3;
  real__t Efd0_3, mu0_3;

  real__t Telec_3;
  real__t mu_3 = 0.;
  real__t Pmech_3 = 0.;
  real__t Efd_3 = 0.;
  real__t VS_3 = 0.;

  real__t Ifd_4, dIfd_dt_4;
  real__t Vd_4, Vq_4;
  real__t Vm_4, Va_4;
  real__t Efd0_4, mu0_4;

  real__t Telec_4;
  real__t mu_4 = 0.;
  real__t Pmech_4 = 0.;
  real__t Efd_4 = 0.;
  real__t VS_4 = 0.;


  real__t Ifd_5, dIfd_dt_5;
  real__t Vd_5, Vq_5;
  real__t Vm_5, Va_5;
  real__t Efd0_5, mu0_5;

  real__t Telec_5;
  real__t mu_5 = 0.;
  real__t Pmech_5 = 0.;
  real__t Efd_5 = 0.;
  real__t VS_5 = 0.;


  real__t Ifd_6, dIfd_dt_6;
  real__t Vd_6, Vq_6;
  real__t Vm_6, Va_6;
  real__t Efd0_6, mu0_6;

  real__t Telec_6;
  real__t mu_6 = 0.;
  real__t Pmech_6 = 0.;
  real__t Efd_6 = 0.;
  real__t VS_6 = 0.;


  real__t Ifd_7, dIfd_dt_7;
  real__t Vd_7, Vq_7;
  real__t Vm_7, Va_7;
  real__t Efd0_7, mu0_7;

  real__t Telec_7;
  real__t mu_7 = 0.;
  real__t Pmech_7 = 0.;
  real__t Efd_7 = 0.;
  real__t VS_7 = 0.;

};

}  // namespace transient_analysis

#endif
