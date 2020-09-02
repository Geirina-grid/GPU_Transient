/******************************************************************************
 * Copyright (c): (2019-2020) GEIRINA
 * All rights reserved.
 * Project: Power system transient simulation
 *
 * - Author:      Peng Wei, peng.wei@geirina.net
 * - Created on:  Jan. 29, 2020
 * - Last Update: Jan. 29, 2020
 *
 * - This library is part of the ODE_solver class. All governor models should be
 *   put in this file.
 *
*******************************************************************************/

#include "ode_solver.cuh"

namespace transient_analysis {

real__t ODE_solver::
process_steam_machine_third_order(const vector_type& x, vector_type& dxdt, const STEAM_DATA& steam,
                                  real__t val_src, const int idx_hp, const int idx_ip, const int idx_lp) {

  /*
  integrate_block(dxdt, val_src,   x[idx_hp], idx_hp, 1. , 1., steam.Tch);
  integrate_block(dxdt, x[idx_hp], x[idx_ip], idx_ip, 1. , 1., steam.Trh);
  integrate_block(dxdt, x[idx_ip], x[idx_lp], idx_lp, 1. , 1., steam.Tco);
  return (steam.Flp * x[idx_lp] + steam.Fip * x[idx_ip] +
          steam.Fhp * (x[idx_hp] + steam.lambda * (x[idx_hp] - x[idx_ip])));
  */
  return 3.14;
}

/** 1型GOV是一种水、火电机组均适用的通用调速器模型 */
real__t ODE_solver::
process_EPRI_GOV_TYPE_I(const vector_type& x, vector_type& dxdt, EPRI_GOV_I_DATA& gov) {
  /*
  assert((gov.gov_type == 1 && gov.TW == 0) || (gov.gov_type == 2 && gov.alpha == 1));

  real__t KmH = parameters.Rate_MW / 100.;
  mu = mu0 + x[delta_mu_idx];
  mu = apply_limiting(mu, gov.mu_Min, gov.mu_Max);
  real__t mult_sign = (mu >= gov.mu_Max || mu <= gov.mu_Min) ? 0. : 1.;

  real__t Pf2_val = gov.Ki * x[delta_mu_idx];
  real__t sigma = gov.K_delta * (omega_ref - x[omega_idx]) - x[GPf_idx] - Pf2_val;
  sigma = apply_dead_band(sigma, 0.5 * gov.dead_band_tol);
  sigma = apply_limiting(sigma, gov.sigma_Min, gov.sigma_Max);
  integrate_block(dxdt, sigma, x[delta_mu_idx], delta_mu_idx, 1., 0., gov.TS);
  */
  /** 软负反馈 */
  //real__t dmudt = dxdt[delta_mu_idx];
  //integrate_block(x, dxdt, mu, dmudt, x[GPf_idx], GPf_idx, 0., gov.Kbeta * gov.Ti, 1., gov.Ti);

  /** the turbine equations */
  //integrate_block(x, dxdt, KmH * mu, KmH * dmudt, x[GP1_idx], GP1_idx, 1., -gov.TW, 1., gov.T0);
  //integrate_block(dxdt, x[GP1_idx], x[GP2_idx], GP2_idx, 1. - gov.alpha, 1., gov.TRH);
  
#if DEBUG
  cout << "***GOV debugging data:***\n";
  cout << "gov.type = " << gov.gov_type << endl;
  cout << "KmH = " << KmH << endl;
  cout << "sigma = " << sigma << endl;
  cout << "omega = " << x[omega_idx] << endl;
  cout << "mu_Max = " << gov.mu_Max << endl;
  cout << "GPf1 = " << x[GPf_idx] << endl;
  cout << "GPf2 = " << Pf2_val << endl;
  cout << "mu0 = " << mu0 << endl;
  cout << "mu = " << mu << endl;
  cout << "delta_mu = " << x[delta_mu_idx] << endl;
  cout << "KmH * mu = " << KmH * mu << endl;
  cout << "GP1 = " << x[GP1_idx] << endl;
  cout << "GP2 = " << x[GP2_idx] << endl;
  cout << "d_mu_dt = " << dmudt << endl;
  cout << "d_GP1_dt = " << dxdt[GP1_idx] << endl;
  cout << "d_GP2_dt = " << dxdt[GP2_idx] << endl;
  cout << "d_delta_mu_dt = " << dxdt[delta_mu_idx] << endl;
  cout << "d_Gpf_dt = " << dxdt[GPf_idx] << endl << endl;
#endif
  
  //return gov.alpha * x[GP1_idx] + x[GP2_idx];
  return 3.14;
}

/** electro_hydraulic_servo 电液伺服机构 is used for Governors 3, 4, and 5 */
real__t ODE_solver::
electro_hydraulic_servo(const vector_type& x, vector_type& dxdt, EHS_DATA& ehs,
                        int idx_src, int idx_pid, int idx_feedback, int idx_dst) {
  /** 电液转换PID模块 */
  /*
  real__t val_pid_in = x[idx_src] - x[idx_feedback];
  real__t div_pid_in = dxdt[idx_src] - dxdt[idx_feedback];
  real__t pid_output = process_PID_block(x, dxdt, val_pid_in, div_pid_in, idx_pid, ehs.pid);
  pid_output = apply_limiting(pid_output, ehs.VEL_Close, ehs.VEL_Open);
  pid_output *= ehs.is_open ? 1. / ehs.TO : 1. / ehs.TC;
  */
  /** 油动机 */
  /*
  dxdt[idx_dst] = pid_output;
  real__t val = apply_limiting(x[idx_dst], ehs.P_Min, ehs.P_Max);
  */
  /** LVDT */
  //integrate_block(dxdt, val, x[idx_feedback], idx_feedback, 1., 1., ehs.T2);
  
  //return val;
   return 3.14;
}

///** 2型调速器属于汽轮机机械液压式调速系统模型,由“液压调节系统”、“汽轮机模型”、“主汽压力变化模型”组成 */
//void ODE_solver::
//process_EPRI_GOV_TYPE_II(const vector_type& x, vector_type& dxdt, EPRI_GOV_II_DATA& gov) {
//  /** 液压调节系统: input - delta_omega, output - delta_mu */
//  real__t delta_omega = omega_ref - x[omega_idx];
//  integrate_block(dxdt, delta_omega, x[GP1_idx], GP1_idx, gov.K, 1., gov.Tr); //转速测量
//
//  delta_omega = apply_dead_band(omega_ref - x[omega_idx], 0.5 * gov.dead_band_tol);
//
//  real__t val_pid_in = -x[GP1_idx] + Pe - Pe_ref;
//  real__t delta_mu = process_PID_block(x, dxdt, val_pid_in, dxdt[GP1_idx], GP2_idx, gov.pid_load);
//
//  /** 主汽压力变化模型 */
//  /** Note: 在现阶段使用的调速器模型中，由于主汽压力变化模型的响应时间很长，
//   *        因此在一般的暂态稳定分析计算中很少考虑，此时认为主汽压力为1
//   */
//  real__t P_pressure = 1.;
//
//  /** 汽轮机模型: input - mu (gate opening, 调门开度), output - P_mech (mechanical power, 机械功率) */
//  mu = P_pressure * (mu_zero + delta_mu);
//  real__t div_mu = P_pressure * (dxdt[GP2_idx] - gov.pid_load.Kp * dxdt[GP1_idx]);
//  real__t P_mech = process_steam_machine_third_order(x, dxdt, gov.steam, mu, hp_idx, ip_idx, lp_idx);
//
//  return P_mech;
//}

/** 3型GOV属于汽轮机电气液压式调速系统模型 */
real__t ODE_solver::
process_EPRI_GOV_TYPE_III(const vector_type& x, vector_type& dxdt, EPRI_GOV_III_DATA& gov) {
  /** 电液调节系统1: input - delta_omega, output - mu */
  //real__t delta_omega = apply_dead_band(omega_ref - x[omega_idx], 0.5 * gov.dead_band_tol);
  //integrate_block(dxdt, delta_omega, x[GP1_idx], GP1_idx, gov.K1, 1., gov.T1); //转速测量
  
  //real__t val_pid_in = -x[GP1_idx] - Telec + Pe_ref;
  //mu = process_PID_block(x, dxdt, val_pid_in, dxdt[GP1_idx], GP2_idx, gov.pid_load);

  /** 电液伺服机构 -- not implemented yet */
  
  /** 主汽压力变化模型 */
  /** Note: 在现阶段使用的调速器模型中，由于主汽压力变化模型的响应时间很长，
   *        因此在一般的暂态稳定分析计算中很少考虑，此时认为主汽压力为1。具体参考电科院文档。
   */
  //real__t P_pressure = 1.;
  
  /** 汽轮机模型: input - mu (gate opening, 调门开度), output - P_mech (mechanical power, 机械功率) */
  //mu *= P_pressure;
//  real__t mu_div = P_pressure * (dxdt[GP2_idx] - gov.pid_load.Kp * dxdt[GP1_idx]);
  //real__t P_mech = process_steam_machine_third_order(x, dxdt, gov.steam, mu, hp_idx, ip_idx, lp_idx);
  
#if DEBUG
  cout << "***GOV debugging data:***\n";
  cout << "Telec = " << Telec << endl;
  cout << "Pe_ref = " << Pe_ref << endl;
  cout << "val_pid_in = " << val_pid_in << endl;
  cout << "div_pid_in = " << dxdt[GP1_idx] << endl;
  cout << "mu = " << mu << endl;
  cout << "P_mech = " << P_mech << endl;
  cout << "GP1 = " << x[GP1_idx] << endl;
  cout << "GP2 = " << x[GP2_idx] << endl;
  cout << "x[hp_idx] = " << x[hp_idx] << endl;
  cout << "x[ip_idx] = " << x[ip_idx] << endl;
  cout << "x[lp_idx] = " << x[lp_idx] << endl;
  cout << "lambda = " << gov.steam.lambda << endl;
  cout << "d_GP1_dt = " << dxdt[GP1_idx] << endl;
  cout << "d_GP2_dt = " << dxdt[GP2_idx] << endl;
  cout << "dxdt[hp_idx] = " << dxdt[hp_idx] << endl;
  cout << "dxdt[ip_idx] = " << dxdt[ip_idx] << endl;
  cout << "dxdt[lp_idx] = " << dxdt[lp_idx] << endl << endl;
#endif
  
  //return P_mech;
  return 3.14;
}

/** 4型GOV属于汽轮机电气液压式调速系统模型 */
real__t ODE_solver::
process_EPRI_GOV_TYPE_IV(const vector_type& x, vector_type& dxdt, EPRI_GOV_IV_DATA& gov) {
  /** 电液调节系统2: input - delta_omega, output - mu */
  /*
  real__t delta_omega = omega_ref - x[omega_idx];
  integrate_block(dxdt, delta_omega, x[GP1_idx], GP1_idx, 1., 1., gov.T1); //转速测量
  
  real__t val_pid_in = - gov.K * apply_dead_band(x[GP1_idx], 0.5 * gov.dead_band_tol);
  real__t div_pid_in = - gov.K * apply_dead_band(dxdt[GP1_idx], 0.5 * gov.dead_band_tol);
  mu = 0.;
  real__t mu_div = 0.;
  
  if (gov.control_option == 1) {          //0,调节器压力反馈控制
    val_pid_in += Pm_ref - Pmech;
    mu = process_PID_block(x, dxdt, val_pid_in, div_pid_in, GP2_idx, gov.pid);
    mu_div = dxdt[GP2_idx] - gov.pid.Kp * dxdt[GP1_idx];
  } else if (gov.control_option == 2) {   //1,DEH开环控制
    mu = mu0 + val_pid_in;
    mu_div = div_pid_in;
  } else if (gov.control_option == 3) {   //2,负荷反馈控制
    mu = gov.K2 * val_pid_in;
    mu_div = gov.K2 * div_pid_in;
    val_pid_in += Pe_ref - Telec;
    mu += process_PID_block(x, dxdt, val_pid_in, div_pid_in, GP2_idx, gov.pid);
    mu_div += dxdt[GP2_idx] - gov.pid.Kp * dxdt[GP1_idx];
  }
  */
  /** 电液伺服机构 -- not implemented */
  
  /** 主汽压力变化模型 */
  /** Note: 在现阶段使用的调速器模型中，由于主汽压力变化模型的响应时间很长，
   *        因此在一般的暂态稳定分析计算中很少考虑，此时认为主汽压力为1
   */
  real__t P_pressure = 1.;
  
  /** 汽轮机模型: input - mu (gate opening, 调门开度), output - P_mech (mechanical power, 机械功率) */
  //mu *= P_pressure;
  //mu_div *= P_pressure;
  //real__t P_mech = process_steam_machine_third_order(x, dxdt, gov.steam, mu, hp_idx, ip_idx, lp_idx);
  
#if DEBUG
  cout << "***GOV debugging data:***\n";
  cout << "Telec = " << Telec << endl;
  cout << "Pe_ref = " << Pe_ref << endl;
  cout << "Pmech = " << Pmech << endl;
  cout << "delta_omega = " << delta_omega << endl;
  cout << "val_pid_in = " << val_pid_in << endl;
  cout << "div_pid_in = " << dxdt[GP1_idx] << endl;
  cout << "mu = " << mu << endl;
  cout << "P_mech = " << P_mech << endl;
  cout << "GP1 = " << x[GP1_idx] << endl;
  cout << "GP2 = " << x[GP2_idx] << endl;
  cout << "x[hp_idx] = " << x[hp_idx] << endl;
  cout << "x[ip_idx] = " << x[ip_idx] << endl;
  cout << "x[lp_idx] = " << x[lp_idx] << endl;
  cout << "lambda = " << gov.steam.lambda << endl;
  cout << "d_GP1_dt = " << dxdt[GP1_idx] << endl;
  cout << "d_GP2_dt = " << dxdt[GP2_idx] << endl;
  cout << "dxdt[hp_idx] = " << dxdt[hp_idx] << endl;
  cout << "dxdt[ip_idx] = " << dxdt[ip_idx] << endl;
  cout << "dxdt[lp_idx] = " << dxdt[lp_idx] << endl << endl;
#endif
  
  //return P_mech;
  return 3.14;
}

/** 5型GOV属于汽轮机电气液压式调速系统模型 */
real__t ODE_solver::
process_EPRI_GOV_TYPE_V(const vector_type& x, vector_type& dxdt, EPRI_GOV_V_DATA& gov) {
  /** 电液调节系统3: input - delta_omega, output - mu */
   /*
  real__t delta_omega = omega_ref - x[omega_idx];
  integrate_block(dxdt, delta_omega, x[GP1_idx], GP1_idx, 1., 1., gov.T1); //转速测量
  mu = 0.;
  real__t mu_div = 0.;
  real__t val_pid_in = 0.;
  real__t div_pid_in = 0;
  
  if (gov.control_option == 0) {          //0,CCS自动控制, ****** TO DO
    mu = Pe_ref - gov.K1 * x[GP1_idx];
    mu_div = - gov.K1 * dxdt[GP1_idx];
  } else if (gov.control_option == 1) {   //1,负荷开环控制
    mu = Pe_ref - gov.K1 * x[GP1_idx];
    mu_div = - gov.K1 * dxdt[GP1_idx];
  } else if (gov.control_option == 2) {   //2,带主汽压力修正负荷控制
    val_pid_in = Pe_ref - Pmech - gov.K1 * x[GP1_idx];
    div_pid_in = - gov.K1 * dxdt[GP1_idx];
    mu = process_PID_block(x, dxdt, val_pid_in, div_pid_in, GP2_idx, gov.pid);
    mu_div = dxdt[GP2_idx] - gov.pid.Kp * dxdt[GP1_idx];
  }
  */
  /** 电液伺服机构 -- not implemented */
  
  /** 主汽压力变化模型 */
  /** Note: 在现阶段使用的调速器模型中，由于主汽压力变化模型的响应时间很长，
   *        因此在一般的暂态稳定分析计算中很少考虑，此时认为主汽压力为1
   */
  real__t P_pressure = 1.;
  
  /** 汽轮机模型: input - mu (gate opening, 调门开度), output - P_mech (mechanical power, 机械功率) */
   //mu *= P_pressure;
   //mu_div *= P_pressure;
   //real__t P_mech = process_steam_machine_third_order(x, dxdt, gov.steam, mu, hp_idx, ip_idx, lp_idx);
  
#if DEBUG
  cout << "*** GOV debugging data: ***\n";
  cout << "Telec = " << Telec << endl;
  cout << "Pe_ref = " << Pe_ref << endl;
  cout << "Pmech = " << Pmech << endl;
  cout << "delta_omega = " << delta_omega << endl;
  cout << "val_pid_in = " << val_pid_in << endl;
  cout << "div_pid_in = " << dxdt[GP1_idx] << endl;
  cout << "mu = " << mu << endl;
  cout << "mu_div = " << mu_div << endl;
  cout << "GP1 = " << x[GP1_idx] << endl;
  cout << "GP2 = " << x[GP2_idx] << endl;
  cout << "d_GP1_dt = " << dxdt[GP1_idx] << endl;
  cout << "d_GP2_dt = " << dxdt[GP2_idx] << endl;
#endif
  
  //return P_mech;
   return 3.14;
}

/**  */
/** GOV Type 6 is not used so far, thus not implemented yet. */
//void ODE_solver::
//process_EPRI_GOV_TYPE_VI(const vector_type& x, vector_type& dxdt, EPRI_GOV_VI_DATA& gov) {
//
//}

/** 7型GOV是水轮机调速器模型，包括调节系统模型、液压系统模型以及水轮机模型 */
real__t ODE_solver::
process_EPRI_GOV_TYPE_VII(const vector_type& x, vector_type& dxdt, EPRI_GOV_VII_DATA& gov) {
  /** 调节系统模型: input - delta_omega, output - mu */
  /*
  integrate_block(dxdt, x[omega_idx], x[GP1_idx], GP1_idx, 1., 1., gov.TR1); //转速测量
  real__t delta_omega = apply_dead_band(omega_ref - x[GP1_idx], gov.dead_band1);
  delta_omega = apply_limiting(delta_omega, gov.dead_band_Min1, gov.dead_band_Max1);
  real__t val_pid_in = delta_omega * gov.Kw;
  real__t div_pid_in = -apply_dead_band(dxdt[GP1_idx], gov.dead_band1) * gov.Kw;
  */
  /** PID block */
   /*
  real__t val1 = apply_limiting(gov.Kp * val_pid_in, gov.PRO_Min, gov.PRO_Max);
  integrate_block(x, dxdt, val_pid_in, div_pid_in, x[GP2_idx], GP2_idx, 0., gov.Kd, 1., gov.Td);
  real__t val2 = x[GP2_idx];
  dxdt[GP3_idx] = apply_limiting(gov.Ki * val_pid_in, 10. * gov.I_Min, 10. * gov.I_Max);
  real__t val3 = x[GP3_idx];
  real__t val_total = val1 + val2 + val3;
  mu = apply_limiting(val_total, 10. * gov.PID_Min, 10. * gov.PID_Max);
  mu = apply_limiting(mu, gov.Ratelimn, gov.Ratelimp);
  real__t mu_div = div_pid_in * gov.Kp + dxdt[GP2_idx] + dxdt[GP3_idx];
  */
  /** 液压系统 -- not implemented */
    
  /** 水轮机模型: input - mu (gate opening, 调门开度), output - P_mech (mechanical power, 机械功率) */
  //integrate_block(x, dxdt, mu, mu_div, x[GPm_idx], GPm_idx, 1., -gov.Tw, 1., 0.5 * gov.Tw);
  
#if DEBUG
  cout << "*** GOV debugging data: ***\n";
  cout << "Telec = " << Telec << endl;
  cout << "Pe_ref = " << Pe_ref << endl;
  cout << "Pmech = " << Pmech << endl;
  cout << "delta_omega = " << delta_omega << endl;
  cout << "val_pid_in = " << val_pid_in << endl;
  cout << "div_pid_in = " << dxdt[GP1_idx] << endl;
  cout << "mu = " << mu << endl;
  cout << "mu_div = " << mu_div << endl;
  cout << "GP1 = " << x[GP1_idx] << endl;
  cout << "GP2 = " << x[GP2_idx] << endl;
  cout << "GP3 = " << x[GP3_idx] << endl;
  cout << "x[GPm_idx] = " << x[GPm_idx] << endl;
  cout << "d_GP1_dt = " << dxdt[GP1_idx] << endl;
  cout << "d_GP2_dt = " << dxdt[GP2_idx] << endl;
  cout << "d_GP3_dt = " << dxdt[GP3_idx] << endl;
  cout << "dxdt[GPm_idx] = " << dxdt[GPm_idx] << endl << endl;
#endif
  
  return x[GPm_idx];
}

/** 8型GOV是水轮机调速器模型，包括调节系统模型、液压系统模型以及水轮机模型 */
real__t ODE_solver::
process_EPRI_GOV_TYPE_VIII(const vector_type& x, vector_type& dxdt, EPRI_GOV_VIII_DATA& gov) {
  /** 调节系统模型: input - delta_omega, output - delta_mu */
  /*
  integrate_block(dxdt, x[omega_idx], x[GP1_idx], GP1_idx, 1., 1., gov.TR1); //转速测量
  real__t delta_omega = apply_dead_band(omega_ref - x[GP1_idx], gov.dead_band_p);
  delta_omega = apply_limiting(delta_omega, gov.dead_band_Min1, gov.dead_band_Max1);
  delta_omega = apply_limiting(delta_omega, gov.Rtd0, gov.Rti0);
  real__t val_pid_in = delta_omega * gov.Kw;
  real__t div_pid_in = -x[GP1_idx] * gov.Kw;
  */
  /** PID block */
  /*
  real__t val1 = gov.Kp2 * val_pid_in;
  integrate_block(x, dxdt, val_pid_in, div_pid_in, x[GP2_idx], GP2_idx, 0., gov.Kd2, 1., gov.Td2);
  real__t val2 = x[GP2_idx];
  dxdt[GP3_idx] = gov.Ki2 * val_pid_in;
  real__t val3 = apply_limiting(x[GP3_idx], 10. * gov.I_Min2, 10. * gov.I_Max2);
  real__t val_total = val1 + val2 + val3;
  mu = apply_limiting(val_total, 10. * gov.PID_Min2, 10. * gov.PID_Max2);
  mu = apply_limiting(mu, gov.Rtd1, gov.Rti1);
  real__t mu_div = -gov.Kw * gov.Kp2 * dxdt[GP1_idx] + dxdt[GP2_idx] + dxdt[GP3_idx];
  */
  /** 液压系统 -- not implemented */
    
  /** 水轮机模型: input - mu (gate opening, 调门开度), output - P_mech (mechanical power, 机械功率) */
  //integrate_block(x, dxdt, mu, mu_div, x[GPm_idx], GPm_idx, 1., -gov.Tw, 1., 0.5 * gov.Tw);
  
#if DEBUG
  cout << "*** GOV debugging data: ***\n";
  cout << "Telec = " << Telec << endl;
  cout << "Pe_ref = " << Pe_ref << endl;
  cout << "Pmech = " << Pmech << endl;
  cout << "delta_omega = " << delta_omega << endl;
  cout << "val_pid_in = " << val_pid_in << endl;
  cout << "div_pid_in = " << dxdt[GP1_idx] << endl;
  cout << "mu = " << mu << endl;
  cout << "mu_div = " << mu_div << endl;
  cout << "GP1 = " << x[GP1_idx] << endl;
  cout << "GP2 = " << x[GP2_idx] << endl;
  cout << "GP3 = " << x[GP3_idx] << endl;
  cout << "x[GPm_idx] = " << x[GPm_idx] << endl;
  cout << "d_GP1_dt = " << dxdt[GP1_idx] << endl;
  cout << "d_GP2_dt = " << dxdt[GP2_idx] << endl;
  cout << "d_GP3_dt = " << dxdt[GP3_idx] << endl;
  cout << "dxdt[GPm_idx] = " << dxdt[GPm_idx] << endl << endl;
#endif
  
  return x[GPm_idx];
}

/** 9型GOV属于汽轮机电气液压式调速系统模型,在4型调速器的基础上改造而来 */
real__t ODE_solver::
process_EPRI_GOV_TYPE_IX(const vector_type& x, vector_type& dxdt, EPRI_GOV_IX_DATA& gov) {
    /** 电液调节系统2: input - delta_omega, output - mu */
  /*
  real__t delta_omega = omega_ref - x[omega_idx];
  integrate_block(dxdt, delta_omega, x[GP1_idx], GP1_idx, 1., 1., gov.T1); //转速测量
  
  real__t val_pid_in = - gov.K * apply_dead_band(x[GP1_idx], 0.5 * gov.dead_band_tol);
  real__t div_pid_in = - gov.K * apply_dead_band(dxdt[GP1_idx], 0.5 * gov.dead_band_tol);
  mu = 0.;
  real__t mu_div = 0.;
  
  if (gov.control_option == 1) {          //0,调节器压力反馈控制
    val_pid_in += Pm_ref - Pmech;
    mu = process_PID_block(x, dxdt, val_pid_in, div_pid_in, GP2_idx, gov.pid);
    mu_div = dxdt[GP2_idx] - gov.pid.Kp * dxdt[GP1_idx];
  } else if (gov.control_option == 2) {   //1,DEH开环控制
    mu = mu0 + val_pid_in;
    mu_div = div_pid_in;
  } else if (gov.control_option == 3) {   //2,负荷反馈控制
    mu = gov.K2 * val_pid_in;
    mu_div = gov.K2 * div_pid_in;
    val_pid_in += Pe_ref - Telec;
    mu += process_PID_block(x, dxdt, val_pid_in, div_pid_in, GP2_idx, gov.pid);
    mu_div += dxdt[GP2_idx] - gov.pid.Kp * dxdt[GP1_idx];
  }
  */
  /** 电液伺服机构 -- not implemented */
  
  /** 主汽压力变化模型 */
  /** Note: 在现阶段使用的调速器模型中，由于主汽压力变化模型的响应时间很长，
   *        因此在一般的暂态稳定分析计算中很少考虑，此时认为主汽压力为1
   */
  //real__t P_pressure = 1.;
  
  /** 汽轮机模型: input - mu (gate opening, 调门开度), output - P_mech (mechanical power, 机械功率) */
  //mu *= P_pressure;
  //mu_div *= P_pressure;
  //real__t P_mech = process_steam_machine_third_order(x, dxdt, gov.steam, mu, hp_idx, ip_idx, lp_idx);
  
#if DEBUG
  cout << "mu = " << mu << endl;
  cout << "Pmech = " << P_mech << endl;
  cout << "d_GP1_dt = " << dxdt[GP1_idx] << endl;
  cout << "d_GP2_dt = " << dxdt[GP2_idx] << endl;
#endif
  
  //return P_mech;
  return 3.14;
}

real__t ODE_solver::
process_GOV_BPAGG(const vector_type& x, vector_type& dxdt, GOV_BPAGG_DATA& gov) {
  /*
  real__t val_delta_omega = gov.R * (x[omega_idx] - omega_ref);
  real__t div_delta_omega = gov.R * dxdt[omega_idx];
  integrate_block(x, dxdt, val_delta_omega, div_delta_omega, x[GP1_idx], GP1_idx, 1., gov.T2, 1., gov.T1);

  real__t delta_Pmech = min(-x[GP1_idx] + Pm_ref, 20 * gov.Pmax); // Pm_ref???
  integrate_block(dxdt, delta_Pmech, x[GP2_idx], GP2_idx, 1., 1., gov.T3);
  integrate_block(dxdt, x[GP2_idx], x[GP3_idx], GP3_idx, 1., 1., gov.T4);
  integrate_block(x, dxdt, GP3_idx, GPm_idx, 1., gov.F * gov.T5, 1., gov.T5);
  
  return x[GPm_idx];
  */
  return 3.14;
}

}  // namespace transient_analysis
