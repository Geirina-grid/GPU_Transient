/******************************************************************************
 * Copyright (c): (2019-2020) GEIRINA
 * All rights reserved.
 * Project: Power system transient simulation
 *
 * - Author:      Peng Wei, peng.wei@geirina.net
 * - Created on:  Jan. 29, 2020
 * - Last Update: Jan. 29, 2020
 *
 * - This library is part of the ODE_solver class. All generator models should be
 *   put in this file.
 *
*******************************************************************************/

#include "ode_solver.cuh"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/set_operations.h>
#include <thrust/execution_policy.h>
namespace transient_analysis {
//vector<real__t>
//typedef vector<double> value_type;
//typedef double value_type;
//typedef thrust::device_vector<value_type> d_vector_type;

//struct EPRI_GEN_DATA_D {
//  string bus_name;
//  uint__t gen_id;
//  uint__t bus_id;
//  real__t Xd, Xdp, Xdpp, Xq, Xqp, Xqpp;
//  real__t X2, Ra;
//  real__t Td0p, Td0pp, Tq0p, Tq0pp, TJ;
//  real__t a, b, n;
//  real__t D;
//};
/*
struct saxpy_functor : public thrust::binary_function<float,float>
{
    const float a;

    //saxpy_functor(float _a) : a(_a) {}

    __host__ __device__
        float operator()(const float& x, const float& y) const { 
            return x - y;
        }
};
*/
/** 1型同步机为考虑Eq'电势恒定的2阶模型。 */
void ODE_solver::
//process_EPRI_GEN_TYPE_I(const state_type& x, state_type& dxdt, EPRI_GEN_DATA& gen) {
process_EPRI_GEN_TYPE_I_D(const d_vector_type& x, d_vector_type& dxdt, real__t TJ, real__t D) {
  
  real__t omega_diff = x[omega_idx] - omega_ref;
  //thrust::device_vector<value_type> d_x = x;
  //thrust::device_vector<value_type> d_dxdt = dxdt;
  //thrust::device_vector<double> d_omega_idx = omega_idx;
  //double d_omega_ref = omega_ref;
  //real__t d_omega_diff = 0;
  //real__t d_omega_diff = saxpy_functor(d_x[omega_idx], d_omega_ref);
  //thrust::transform(d_x.begin(), d_x.end(), d_omega_ref, d_omega_diff, thrust::minus<double>());

  dxdt[omega_idx] = (TJ < EPS)
                    ? 0
                    : (Pmech - Telec - D * omega_diff) / TJ;
  
  dxdt[delta_idx] = 2 * PI * freq_ref * omega_diff;
  /*
  printf("x[omega_idx] in solver=%f\n", x[omega_idx]);
  printf("x=%f\n", x);
  printf("dxdt=%f\n", dxdt);
  printf("gen=%f\n", TJ);
  printf("gen=%f\n", D);
  printf("Pmech=%f\n", Pmech);
  printf("Telec=%f\n", Telec);
  printf("omega_idx=%f\n", omega_idx);
  printf("omega_ref=%f\n", omega_ref);
  printf("EPS=%f\n", EPS);
  printf("delta_idx=%f\n", delta_idx);
  printf("PI=%f\n", PI);
  printf("freq_ref=%f\n", freq_ref);
  printf("omega_diff=%f\n", omega_diff);
  printf("Final result=%f",  dxdt[delta_idx]);
  */

//#if DEBUG
//  cout << "\n*** GEN debugging data: ***\n";
//  cout << "Edp = " << x[Edp_idx] << endl;
//  cout << "Eqp = " << x[Eqp_idx] << endl;
//  cout << "Pmech = " << Pmech << endl;
//  cout << "Telec = " << Telec << endl;
//  cout << "Xq = " << gen.Xq << ", Xqp = " << gen.Xqp << ", Xqpp = " << gen.Xqpp << endl;
//  cout << "Xd = " << gen.Xd << ", Xdp = " << gen.Xdp << ", Xdpp = " << gen.Xdpp << endl;
//  cout << "d_omega_dt = " << dxdt[omega_idx] << endl;
//  cout << "d_delta_dt = " << dxdt[delta_idx] << endl << endl;
//#endif
}

void ODE_solver::
//process_EPRI_GEN_TYPE_I(const state_type& x, state_type& dxdt, EPRI_GEN_DATA& gen) {
process_EPRI_GEN_TYPE_I(const vector_type& x, vector_type& dxdt, EPRI_GEN_DATA& gen) {

  real__t omega_diff = x[omega_idx] - omega_ref;
  //thrust::device_vector<value_type> d_x = x;
  //thrust::device_vector<value_type> d_dxdt = dxdt;
  //thrust::device_vector<double> d_omega_idx = omega_idx;
  //double d_omega_ref = omega_ref;
  //real__t d_omega_diff = 0;
  //real__t d_omega_diff = saxpy_functor(d_x[omega_idx], d_omega_ref);
  //thrust::transform(d_x.begin(), d_x.end(), d_omega_ref, d_omega_diff, thrust::minus<double>());

  dxdt[omega_idx] = (gen.TJ < EPS)
                    ? 0
                    : (Pmech - Telec - gen.D * omega_diff) / gen.TJ;

  dxdt[delta_idx] = 2 * PI * freq_ref * omega_diff;
  /*
  printf("x=%d\n", x);
  printf("dxdt=%d\n", dxdt);
  printf("gen=%d\n", TJ);
  printf("gen=%d\n", D);
  printf("Pmech=%d\n", Pmech);
  printf("Telec=%d\n", Telec);
  printf("omega_idx=%d\n", omega_idx);
  printf("omega_ref=%d\n", omega_ref);
  printf("EPS=%d\n", EPS);
  printf("delta_idx=%d\n", delta_idx);
  printf("PI=%d\n", PI);
  printf("freq_ref=%d\n", freq_ref);
  */

//#if DEBUG
//  cout << "\n*** GEN debugging data: ***\n";
//  cout << "Edp = " << x[Edp_idx] << endl;
//  cout << "Eqp = " << x[Eqp_idx] << endl;
//  cout << "Pmech = " << Pmech << endl;
//  cout << "Telec = " << Telec << endl;
//  cout << "Xq = " << gen.Xq << ", Xqp = " << gen.Xqp << ", Xqpp = " << gen.Xqpp << endl;
//  cout << "Xd = " << gen.Xd << ", Xdp = " << gen.Xdp << ", Xdpp = " << gen.Xdpp << endl;
//  cout << "d_omega_dt = " << dxdt[omega_idx] << endl;
//  cout << "d_delta_dt = " << dxdt[delta_idx] << endl << endl;
//#endif
}
/** 2型同步机为考虑Eq'电势变化的3阶模型。 */
void ODE_solver::
process_EPRI_GEN_TYPE_II(const vector_type& x, vector_type& dxdt, EPRI_GEN_DATA& gen) {
//  Telec = x[Eqp_idx] * Iq + x[Edp_idx] * Id - (gen.Xdp - gen.Xqp) * Id * Iq;

  dxdt[omega_idx] = (gen.TJ < EPS)
                    ? 0
                    : (Pmech - Telec - gen.D * (x[omega_idx] - omega_ref)) / gen.TJ;
  
  dxdt[delta_idx] = 2 * PI * freq_ref * (x[omega_idx] - omega_ref);
  
  real__t KG = 1. + gen.b / gen.a * pow(x[Eqp_idx], gen.n - 1);
  dxdt[Eqp_idx] = (gen.Td0p < EPS)
                  ? 0.
                  : (Efd - x[Eqp_idx] - (gen.Xd - gen.Xdp) * Id - (KG - 1.) * x[Eqp_idx]) / gen.Td0p;
  
#if DEBUG
  cout << "\n*** GEN debugging data: ***\n";
  cout << "Edp = " << x[Edp_idx] << endl;
  cout << "Eqp = " << x[Eqp_idx] << endl;
  cout << "Pmech = " << Pmech << endl;
  cout << "Telec = " << Telec << endl;
  cout << "Xq = " << gen.Xq << ", Xqp = " << gen.Xqp << ", Xqpp = " << gen.Xqpp << endl;
  cout << "Xd = " << gen.Xd << ", Xdp = " << gen.Xdp << ", Xdpp = " << gen.Xdpp << endl;
  cout << "d_omega_dt = " << dxdt[omega_idx] << endl;
  cout << "d_delta_dt = " << dxdt[delta_idx] << endl;
  cout << "d_Eqp_dt = " << dxdt[Eqp_idx] << endl;
#endif
}

/** 3型同步机为考虑Eq', Eq", Ed"电势变化的5阶模型。 本模型适合于凸极转子(水轮)发电机的详细模型*/
void ODE_solver::
process_EPRI_GEN_TYPE_III(const vector_type& x, vector_type& dxdt, EPRI_GEN_DATA& gen) {
//  Telec = x[Eqpp_idx] * Iq + x[Edpp_idx] * Id - (gen.Xdpp - gen.Xqpp) * Id * Iq;

  dxdt[omega_idx] = (gen.TJ < EPS)
                    ? 0
                    : (Pmech - Telec - gen.D * (x[omega_idx] - omega_ref)) / gen.TJ;
  
  dxdt[delta_idx] = 2 * PI * freq_ref * (x[omega_idx] - omega_ref);

  real__t KG = 1. + gen.b / gen.a * pow(x[Eqp_idx], gen.n - 1);
  dxdt[Eqp_idx] = (gen.Td0p < EPS)
                  ? 0.
                  : (Efd - x[Eqp_idx] - (gen.Xd - gen.Xdp) * Id - (KG - 1.) * x[Eqp_idx]) / gen.Td0p;

  dxdt[Eqpp_idx] = (gen.Td0pp < EPS)
                   ? 0.
                   : (-x[Eqpp_idx] - (gen.Xdp - gen.Xdpp) * Id + x[Eqp_idx]) / gen.Td0pp + dxdt[Eqp_idx];

  dxdt[Edpp_idx] = (gen.Tq0pp < EPS)
                   ? 0.
                   : (-x[Edpp_idx] + (gen.Xqp - gen.Xqpp) * Iq + x[Edp_idx]) / gen.Tq0pp + dxdt[Edp_idx];
  
#if DEBUG
  cout << "\n*** GEN debugging data: ***\n";
  cout << "Edp = " << x[Edp_idx] << endl;
  cout << "Eqp = " << x[Eqp_idx] << endl;
  cout << "Edpp = " << x[Edpp_idx] << endl;
  cout << "Eqpp = " << x[Eqpp_idx] << endl;
  cout << "Pmech = " << Pmech << endl;
  cout << "Telec = " << Telec << endl;
  cout << "Xq = " << gen.Xq << ", Xqp = " << gen.Xqp << ", Xqpp = " << gen.Xqpp << endl;
  cout << "Xd = " << gen.Xd << ", Xdp = " << gen.Xdp << ", Xdpp = " << gen.Xdpp << endl;
  cout << "d_omega_dt = " << dxdt[omega_idx] << endl;
  cout << "d_delta_dt = " << dxdt[delta_idx] << endl;
  cout << "d_Eqp_dt = " << dxdt[Eqp_idx] << endl;
  cout << "d_Edpp_dt = " << dxdt[Edpp_idx] << endl;
  cout << "d_Eqpp_dt = " << dxdt[Eqpp_idx] << endl;
#endif
}

/** 6型同步机为考虑Eq", Ed", Eq', Ed'电势均发生变化的6阶同步机模型，即原始的同步电机方程。适用于任何计算精度要求较高的场合。*/
void ODE_solver::
process_EPRI_GEN_TYPE_VI(const vector_type& x, vector_type& dxdt, EPRI_GEN_DATA& gen) {
//  Telec = x[Eqpp_idx] * Iq + x[Edpp_idx] * Id - (gen.Xdpp - gen.Xqpp) * Id * Iq;
  printf("x[omega_idx] in process_EPRI_GEN_TYPE_VI =%f\n", x[omega_idx]);
  real__t omega_diff = x[omega_idx] - omega_ref;
  
  dxdt[omega_idx] = (gen.TJ < EPS)
                    ? 0
                    : (Pmech - Telec - gen.D * omega_diff) / gen.TJ;
  
  dxdt[delta_idx] = 2 * PI * freq_ref * omega_diff;

  real__t KG = 1. + gen.b / gen.a * pow(x[Eqp_idx], gen.n - 1);
  dxdt[Eqp_idx] = (gen.Td0p < EPS)
                  ? 0.
                  : (Efd - x[Eqp_idx] - (gen.Xd - gen.Xdp) * Id - 0 * (KG - 1.) * x[Eqp_idx]) / gen.Td0p;

  dxdt[Edp_idx] = (gen.Tq0p < EPS)
                  ? 0.
                  : (-x[Edp_idx] + (gen.Xq - gen.Xqp) * Iq) / gen.Tq0p;

  dxdt[Eqpp_idx] = (gen.Td0pp < EPS)
                   ? 0.
                   : (-x[Eqpp_idx] - (gen.Xdp - gen.Xdpp) * Id + x[Eqp_idx]) / gen.Td0pp + dxdt[Eqp_idx];

  dxdt[Edpp_idx] = (gen.Tq0pp < EPS)
                   ? 0.
                   : (-x[Edpp_idx] + (gen.Xqp - gen.Xqpp) * Iq + x[Edp_idx]) / gen.Tq0pp + dxdt[Edp_idx];

#if DEBUG
  cout << "\n*** GEN debugging data: ***\n";
  cout << "Efd = " << Efd << endl;
  cout << "KG = " << KG << endl;
  cout << "Id = " << Id << endl;
  cout << "Iq = " << Iq << endl;
  cout << "omega = " << x[omega_idx] << endl;
  cout << "delta (deg) = " << x[delta_idx] * 180 / PI << endl;
  cout << "Edp = " << x[Edp_idx] << endl;
  cout << "Eqp = " << x[Eqp_idx] << endl;
  cout << "Edpp = " << x[Edpp_idx] << endl;
  cout << "Eqpp = " << x[Eqpp_idx] << endl;
  cout << "Pmech = " << Pmech << endl;
  cout << "Telec = " << Telec << endl;
  cout << "Xq = " << gen.Xq << ", Xqp = " << gen.Xqp << ", Xqpp = " << gen.Xqpp << endl;
  cout << "Xd = " << gen.Xd << ", Xdp = " << gen.Xdp << ", Xdpp = " << gen.Xdpp << endl;
  cout << "d_omega_dt = " << dxdt[omega_idx] << endl;
  cout << "d_delta_dt = " << dxdt[delta_idx] << endl;
  cout << "d_Edp_dt = " << dxdt[Edp_idx] << endl;
  cout << "d_Eqp_dt = " << dxdt[Eqp_idx] << endl;
  cout << "d_Edpp_dt = " << dxdt[Edpp_idx] << endl;
  cout << "d_Eqpp_dt = " << dxdt[Eqpp_idx] << endl << endl;
#endif
}

void ODE_solver::process_GENROU(const vector_type& x, vector_type& dxdt, GENROU_IEEE_DATA& gen) {
  real__t Xl   = gen.Xl;
  real__t Xd   = gen.Xd;
  real__t Xq   = gen.Xq;
  real__t Xdp  = gen.Xdp;
  real__t Xqp  = gen.Xqp;
  real__t Xdpp = gen.Xdpp;
  real__t Xqpp = gen.Xqpp;

  Telec = x[Eqpp_idx] * Iq + x[Edpp_idx] * Id - (Xdpp - Xqpp) * Id * Iq;
  
//  Pmech = x[gov_Pmech_idx];
//  Pmech = x[PT2_idx] + gov.alpha * x[PCH_idx];

//  if (gen.gen_id == 12 || gen.gen_id == 31)
//    printf("id = %2d,  Vq = %+2.6lf,  Vd = %+2.6lf,  Eqpp = %+2.6lf,  Edpp = %+2.6lf,  Iq = %+2.6lf,  Id = %+2.6lf,  Pmech = %+2.6lf,  Telec = %+2.6lf,  diff = %+2.6lf\n",
//           gen.gen_id, Vq, Vd, x[Eqpp_idx], x[Edpp_idx], Iq, Id, Pmech, Telec, Pmech - Telec);
  
//  real__t Efd = x[Efd_idx];
//  Efd = apply_limiting(x[Efd_idx], exc.Efd_Min, exc.Efd_Max);
  real__t delta_omega = (x[omega_idx] - omega_ref);
//  delta_omega = apply_dead_band(delta_omega, 0.001);


  dxdt[omega_idx] = (gen.TJ < EPS)
                    ? 0.
                    : (  (Pmech / x[omega_idx] - Telec)
                       - gen.D * delta_omega / x[omega_idx] ) / gen.TJ;

//  if (gen.gen_id == 31 || gen.gen_id == 11 || gen.gen_id == 62 || gen.gen_id == 48 || gen.gen_id == 53)
//    dxdt[omega_idx] = 0.;
  
  dxdt[delta_idx] = 2 * PI * freq_ref * delta_omega;

  dxdt[Eqp_idx] = (gen.Td0p < EPS)
                  ? 0.
                  : (+ Efd
                     + (Xd - Xdp) / (Xdp - Xl) * x[Eqpp_idx]
                     - (Xd - Xl)  / (Xdp - Xl) * x[Eqp_idx]
                     - (Xd - Xdp) * (Xdpp - Xl) / (Xdp - Xl) * Id ) / gen.Td0p;

  dxdt[Edp_idx] = (gen.Tq0p < EPS)
                  ? 0.
                  : (+ (Xq - Xqp) / (Xqp - Xl) * x[Edpp_idx]
                     - (Xq - Xl)  / (Xqp - Xl) * x[Edp_idx]
                     + (Xq - Xqp) * (Xqpp - Xl) / (Xqp - Xl) * Iq ) / gen.Tq0p;

  dxdt[Eqpp_idx] = (gen.Td0pp < EPS)
                   ? 0.
                   :  (Xdpp - Xl) / (Xdp - Xl) * dxdt[Eqp_idx]
                    + (-x[Eqpp_idx] - (Xdp - Xdpp) * Id + x[Eqp_idx]) / gen.Td0pp;

  dxdt[Edpp_idx] = (gen.Tq0pp < EPS)
                   ? 0.
                   :  (Xqpp - Xl) / (Xqp - Xl) * dxdt[Edp_idx]
                    + (-x[Edpp_idx] + (Xqp - Xqpp) * Iq + x[Edp_idx]) / gen.Tq0pp;
  
#if DEBUG
  cout << "Edp = " << x[Edp_idx] << endl;
  cout << "Eqp = " << x[Eqp_idx] << endl;
  cout << "Edpp = " << x[Edpp_idx] << endl;
  cout << "Eqpp = " << x[Eqpp_idx] << endl;
  cout << "Pmech = " << Pmech << endl;
  cout << "Telec = " << Telec << endl;
  cout << "Xq = " << Xq << ", Xqp = " << Xqp << ", Xqpp = " << Xqpp << endl;
  cout << "Xd = " << Xd << ", Xdp = " << Xdp << ", Xdpp = " << Xdpp << endl;
  cout << "d_omega_dt = " << dxdt[omega_idx] << endl;
  cout << "d_delta_dt = " << dxdt[delta_idx] << endl;
  cout << "d_Edp_dt = " << dxdt[Edp_idx] << endl;
  cout << "d_Eqp_dt = " << dxdt[Eqp_idx] << endl;
  cout << "d_Edpp_dt = " << dxdt[Edpp_idx] << endl;
  cout << "d_Eqpp_dt = " << dxdt[Eqpp_idx] << endl;
#endif
}

}  // namespace transient_analysis
