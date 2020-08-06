/******************************************************************************
 * Copyright (c): (2019-2020) GEIRINA
 * All rights reserved.
 * Project: Power system transient simulation
 *
 * - Author:      Peng Wei, peng.wei@geirina.net
 * - Created on:  Jan. 29, 2020
 * - Last Update: Jan. 29, 2020
 *
 * - This library is part of the ODE_solver class. All excitor (AVR) models
 *   should be put in this file.
 *
*******************************************************************************/

#include "ode_solver.hpp"

namespace transient_analysis {

/** IEEE EXCITOR TYPE 1*/
real__t ODE_solver::
process_EXC_IEEE_I(const vector_type& x, vector_type& dxdt, EXC_IEEE_I_DATA& exc) {
  /** 量测 */
  integrate_block(dxdt, Vm, x[EVa_idx], EVa_idx, 1., 1., exc.TR);

  /** 放大 */
  real__t V_diff = Vm_ref + VS - x[EVa_idx] - x[EVf_idx];
  integrate_block(dxdt, V_diff, x[EVr_idx], EVr_idx, exc.KA, 1., exc.TA);

  /** 励磁机 */
  real__t Vr_val = apply_limiting(x[EVr_idx], 10. * exc.VR_Min, 10. * exc.VR_Max);
  real__t sat_Efd = exc.b * pow(x[EVe_idx] - exc.a, 2);
  real__t V_adj = exc.KE * x[EVe_idx] + sat_Efd;
  integrate_block(dxdt, Vr_val - V_adj, x[EVe_idx], EVe_idx, 1., 0., exc.TE);

  /** 反馈 */
  integrate_block(x, dxdt, EVe_idx, EVf_idx, 0., exc.KF, 1., exc.TF);
  
  return x[EVe_idx];
}

/** 1型AVR为它励式常规励磁系统或采用可控硅调节器的它励式快速励磁系统，即通常具有励磁机的励磁调节系统 */
real__t ODE_solver::
process_EPRI_EXC_TYPE_I(const vector_type& x, vector_type& dxdt, EPRI_EXC_I_DATA& exc) {
  /** 量测环节 */
  real__t V_diff = -Vm + Vm_ref;
  integrate_block(dxdt, V_diff, x[EVr_idx], EVr_idx, exc.Kr, 1., exc.Tr);

  /** 放大环节 */
  V_diff = VS - x[EVf_idx] + x[EVr_idx];
  integrate_block(dxdt, V_diff, x[EVa_idx], EVa_idx, exc.Ka, 1., exc.Ta);

  /** 励磁机 */
  real__t Va_sum = apply_limiting(x[EVa_idx] + Efd0, -10. * exc.Efd_Max, 10. * exc.Efd_Max);
  integrate_block(dxdt, Va_sum, x[EVe_idx], EVe_idx, 1., 1., exc.Te);

  /** 反馈环节 */
  integrate_block(x, dxdt, EVe_idx, EVf_idx, 0., exc.Kf, 1., exc.Tf);
  
  real__t Efd_val = apply_limiting(x[EVe_idx], 10. * exc.Efd_Min, 10. * exc.Efd_Max);
  
#if DEBUG
  cout << "*** AVR debugging data: ***\n";
  cout << "Vm = " << Vm << endl;
  cout << "Vm_ref = " << Vm_ref << endl;
  cout << "VS = " << VS << endl;
  cout << "V_diff = " << V_diff << endl;
  cout << "Va = " << x[EVa_idx] << endl;
  cout << "Vr = " << x[EVr_idx] << endl;
  cout << "Vf = " << x[EVf_idx] << endl;
  cout << "Ve = " << x[EVe_idx] << endl;
  cout << "dVa_dt = " << dxdt[EVa_idx] << endl;
  cout << "dVr_dt = " << dxdt[EVr_idx] << endl;
  cout << "dVf_dt = " << dxdt[EVf_idx] << endl;
  cout << "dVe_dt = " << dxdt[EVe_idx] << endl;
  cout << "Efd0 = " << Efd0 << endl;
  cout << "Efd = " << Efd_val << endl << endl;
#endif
  
  return Efd_val;
}

/** 2型AVR为采用可控硅调节器的自并励和自复励快速励磁系统 */
real__t ODE_solver::
process_EPRI_EXC_TYPE_II(const vector_type& x, vector_type& dxdt, EPRI_EXC_II_DATA& exc) {
  assert(exc.K2 == 0 || exc.K2 == 1);
  
  /** 量测环节 */
  real__t V_diff = -Vm + Vm_ref + VS;
  integrate_block(dxdt, V_diff, x[EV1_idx], EV1_idx, exc.Kr, 1., exc.Tr);

  /** 中间环节 */
  integrate_block(x, dxdt, EV1_idx, EV2_idx, 1., exc.T1, exc.K2, exc.T2);
  integrate_block(x, dxdt, EV2_idx, EV3_idx, 1., exc.T3, 1., exc.T4);

  /** 放大环节 */
  integrate_block(dxdt, x[EV3_idx], x[EVa_idx], EVa_idx, exc.Ka, 1., exc.Ta);
  
  real__t Efd_val = x[EVa_idx] + Efd0;

#if DEBUG
  cout << "*** AVR debugging data: ***\n";
  cout << "Vm = " << Vm << endl;
  cout << "Vm_ref = " << Vm_ref << endl;
  cout << "V_diff = " << V_diff << endl;
  cout << "VS = " << VS << endl;
  cout << "V1 = " << x[EV1_idx] << endl;
  cout << "V2 = " << x[EV2_idx] << endl;
  cout << "V3 = " << x[EV3_idx] << endl;
  cout << "Va = " << x[EVa_idx] << endl;
  cout << "d_V1_dt = " << dxdt[EV1_idx] << endl;
  cout << "d_V2_dt = " << dxdt[EV2_idx] << endl;
  cout << "d_V3_dt = " << dxdt[EV3_idx] << endl;
  cout << "d_Va_dt = " << dxdt[EVa_idx] << endl;
  cout << "Efd = " << Efd_val << endl;
  cout << "Efd_ref = " << Efd0 << endl << endl;
#endif
  
  return apply_limiting(Efd_val, exc.Efd_Min, exc.Efd_Max);
}

/** c.f. P.23 */
real__t ODE_solver::get_exc_input(real__t Xc) {
  real__t angle_diff = -atan2(Vq, Vd) + atan2(Iq, Id); // Voltage angle is the reference
  real__t re_val = Vd + Id * Xc * sin(angle_diff);
  real__t im_val = Vq + Iq * Xc * sin(angle_diff);
  return sqrt(re_val * re_val + im_val * im_val);
}

/** c.f. https://www.neplan.ch/wp-content/uploads/2015/08/Nep_EXCITERS1.pdf P140 */
real__t ODE_solver::get_rectified_regulator(real__t src) {
  return 1.;  // nonlinear function is not used anymore
  
  real__t FEX = 0.;
  if (src <= 0)
    FEX = 1.;
  else if (src <= 0.433)
    FEX = 1. - 0.577 * src;
  else if (src < 0.75)
    FEX = sqrt(0.75 - src * src);
  else if (src <= 1)
    FEX = 1.732 * (1 - src);
  else
    FEX = 0.;
  return FEX;
}

void ODE_solver::
print_for_debug_EXC_III_X(const vector_type& x, vector_type& dxdt, real__t V_in, real__t val_Vdiff,
                          real__t div_Vdiff, real__t Vfe_val, real__t Efd_val) {
#if DEBUG
  cout << "\n*** AVR debugging data: ***\n";
  cout << "Vin = " << V_in << endl;
  cout << "Vm_ref = " << Vm_ref << endl;
  cout << "VS = " << VS << endl;
  cout << "val_Vdiff = " << val_Vdiff << endl;
  cout << "div_Vdiff = " << div_Vdiff << endl;
  cout << "V1 = " << x[EV1_idx] << endl;
  cout << "V2 = " << x[EV2_idx] << endl;
  cout << "V3 = " << x[EV3_idx] << endl;
  cout << "Va = " << x[EVa_idx] << endl;
  cout << "Vr = " << x[EVr_idx] << endl;
  cout << "Ve = " << x[EVe_idx] << endl;
  cout << "Vf = " << x[EVf_idx] << endl;
  cout << "Vfe = " << Vfe_val << endl;
  cout << "d_V1_dt = " << dxdt[EV1_idx] << endl;
  cout << "d_V2_dt = " << dxdt[EV2_idx] << endl;
  cout << "d_V3_dt = " << dxdt[EV3_idx] << endl;
  cout << "d_Va_dt = " << dxdt[EVa_idx] << endl;
  cout << "d_Vr_dt = " << dxdt[EVr_idx] << endl;
  cout << "d_Ve_dt = " << dxdt[EVe_idx] << endl;
  cout << "d_Vf_dt = " << dxdt[EVf_idx] << endl;
  cout << "Efd = " << Efd_val << endl;
  cout << "Efd0 = " << Efd0 << endl << endl;
#endif
}

/** 3型AVR用来模拟由副励磁机向调节器供电的不可控整流交流励磁机励磁系统,主要应用于无刷励磁系统 */
real__t ODE_solver::
process_EPRI_EXC_TYPE_III(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc) {
  real__t Ve_val = apply_limiting(x[EVe_idx], -10000., exc.Ve_Max);
  real__t FEX = get_rectified_regulator(exc.Kc * Ifd / Ve_val);  //Rectified Regulation
  real__t Efd_val = apply_limiting(Ve_val * FEX, -10000., exc.Efd_Max);
  real__t Se = exc.C1 * exp(exc.C2 * Efd_val); /** 励磁机饱和系数 */
  real__t V_adj = Se * Ve_val + Ifd * exc.Kd;
  
  /** 量测环节 */
  real__t V_in = get_exc_input(exc.Xc);
  integrate_block(dxdt, V_in, x[EV1_idx], EV1_idx, 1., 1., exc.Tr);

  /** 串联校正环节 */
  real__t val_Vdiff = -x[EV1_idx] + Vm_ref + VS - x[EVf_idx];
  real__t div_Vdiff = -dxdt[EV1_idx] + dxdt[VS_idx] - dxdt[EVf_idx];
  integrate_block(x, dxdt, val_Vdiff, div_Vdiff, x[EV2_idx], EV2_idx, exc.K, exc.K * exc.T1, exc.Kv, exc.T2);
  integrate_block(x, dxdt, EV2_idx, EV3_idx, 1., exc.T3, 1., exc.T4);

  /** 并联校正环节 */
  real__t Vfe_val = V_adj + exc.Ke * Ve_val;
  real__t Vfe_div = (exc.Ke + Se) * dxdt[EVe_idx] + exc.Kd * dIfd_dt;
  integrate_block(x, dxdt, Vfe_val, Vfe_div, x[EVf_idx], EVf_idx, 0., exc.Kf, 1., exc.Tf);

  /** 功率放大环节 */
  integrate_block(dxdt, x[EV3_idx], x[EVa_idx], EVa_idx, exc.Ka, 1., exc.Ta);

  /** 第二级调节器 */
  real__t Va_val = apply_limiting(x[EVa_idx], exc.Va_Min, exc.Va_Max);
  Va_val -= Vfe_val * exc.KH1;
  integrate_block(dxdt, Va_val, x[EVr_idx], EVr_idx, exc.KB, 1., exc.T5);

  /** 励磁机 */
  real__t Vr_val = apply_limiting(x[EVr_idx], exc.Vr_Min, exc.Vr_Max);
  integrate_block(dxdt, Vr_val - V_adj, x[EVe_idx], EVe_idx, 1., exc.Ke, exc.Te);

  print_for_debug_EXC_III_X(x, dxdt, V_in, val_Vdiff, div_Vdiff, Vfe_val, Efd_val);
  
  return Efd_val;
}

/** 4型AVR用来模拟由副励磁机向调节器供电的不可控整流交流励磁机励磁系统,应用于有刷励磁系统 */
real__t ODE_solver::
process_EPRI_EXC_TYPE_IV(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc) {
  /** 右端项，显式表达 */
  real__t Ve_val = apply_limiting(x[EVe_idx], -10000., exc.Ve_Max);
  real__t FEX = get_rectified_regulator(exc.Kc * Ifd / Ve_val);  //Rectified Regulation
  real__t Efd_val = apply_limiting(Ve_val * FEX, -10000., exc.Efd_Max);
  
  /** 量测环节 */
  real__t V_in = get_exc_input(exc.Xc);
  integrate_block(dxdt, V_in, x[EV1_idx], EV1_idx, 1., 1., exc.Tr);

  /** 串联校正环节 */
  real__t val_Vdiff = -x[EV1_idx] + Vm_ref + VS - x[EVf_idx];
  real__t div_Vdiff = -dxdt[EV1_idx] + dxdt[VS_idx] - dxdt[EVf_idx];
  integrate_block(x, dxdt, val_Vdiff, div_Vdiff, x[EV2_idx], EV2_idx, exc.K, exc.K * exc.T1, exc.Kv, exc.T2);
  integrate_block(x, dxdt, EV2_idx, EV3_idx, 1., exc.T3, 1., exc.T4);
  
  /** 并联校正环节 */
   integrate_block(x, dxdt, EVe_idx, EVf_idx, 0., exc.Kf, 1., exc.Tf);
   
  /** 功率放大环节 */
  integrate_block(dxdt, x[EV3_idx], x[EVa_idx], EVa_idx, exc.Ka, 1., exc.Ta);
  
  /** 第二级调节器 */
  real__t Va_val = apply_limiting(x[EVa_idx], exc.Va_Min, exc.Va_Max);
  Va_val -= Efd_val * exc.KH1;
  integrate_block(dxdt, Va_val, x[EVr_idx], EVr_idx, exc.KB, 1., exc.T5);
  
  /** 励磁机 */
  real__t Se = exc.C1 * exp(exc.C2 * Efd_val); /** 励磁机饱和系数 */
  real__t V_adj = Se * Ve_val + Ifd * exc.Kd;
  real__t Vr_val = apply_limiting(x[EVr_idx], exc.Vr_Min, exc.Vr_Max);
  integrate_block(dxdt, Vr_val - V_adj, x[EVe_idx], EVe_idx, 1., exc.Ke, exc.Te);

  print_for_debug_EXC_III_X(x, dxdt, V_in, val_Vdiff, div_Vdiff, Efd_val, Efd_val);
  
  return Efd_val;
}

/** 5型AVR用来模拟没有副励磁机的交流励磁机不可控整流器励磁系统,应用于无刷励磁系统 */
real__t ODE_solver::
process_EPRI_EXC_TYPE_V(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc) {
  real__t Ve_val = apply_limiting(x[EVe_idx], -10000., exc.Ve_Max);
  real__t FEX = get_rectified_regulator(exc.Kc * Ifd / Ve_val);  //Rectified Regulation
  real__t Efd_val = apply_limiting(Ve_val * FEX, -10000., exc.Efd_Max);
  real__t Se = exc.C1 * exp(exc.C2 * Efd_val); /** 励磁机饱和系数 */
  real__t V_adj = Se * Ve_val + Ifd * exc.Kd;
  
  /** 量测环节 */
  real__t V_in = get_exc_input(exc.Xc);
  integrate_block(dxdt, V_in, x[EV1_idx], EV1_idx, 1., 1., exc.Tr);
  
  /** 串联校正环节 */
  real__t val_Vdiff = -x[EV1_idx] + Vm_ref + VS - x[EVf_idx];
  real__t div_Vdiff = -dxdt[EV1_idx] + dxdt[VS_idx] - dxdt[EVf_idx];
  integrate_block(x, dxdt, val_Vdiff, div_Vdiff, x[EV2_idx], EV2_idx, exc.K, exc.K * exc.T1, exc.Kv, exc.T2);
  integrate_block(x, dxdt, EV2_idx, EV3_idx, 1., exc.T3, 1., exc.T4);

  /** 并联校正环节 */
  real__t Vfe_val = V_adj + exc.Ke * Ve_val;
  real__t Vfe_div = (exc.Ke + Se) * dxdt[EVe_idx] + exc.Kd * dIfd_dt;
  integrate_block(x, dxdt, Vfe_val, Vfe_div, x[EVf_idx], EVf_idx, 0., exc.Kf, 1., exc.Tf);
  
  /** 功率放大环节 */
  integrate_block(dxdt, x[EV3_idx], x[EVa_idx], EVa_idx, exc.Ka, 1., exc.Ta);
  
  /** 第二级调节器 */
  real__t Va_val = apply_limiting(x[EVa_idx], exc.Va_Min, exc.Va_Max);
  Va_val -= Vfe_val * exc.KH1;
  integrate_block(dxdt, Va_val, x[EVr_idx], EVr_idx, exc.KB, 1., exc.T5);
  
  /** 励磁机 */
  real__t Vr_val = apply_limiting(x[EVr_idx], Vm * exc.Vr_Min, Vm * exc.Vr_Max);
  integrate_block(dxdt, Vr_val - V_adj, x[EVe_idx], EVe_idx, 1., exc.Ke, exc.Te);

  print_for_debug_EXC_III_X(x, dxdt, V_in, val_Vdiff, div_Vdiff, Vfe_val, Efd_val);
  
  return Efd_val;
}

/** 6型AVR用来模拟没有副励磁机的交流励磁机不可控整流器励磁系统,应用于有刷励磁系统 */
real__t ODE_solver::
process_EPRI_EXC_TYPE_VI(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc) {
  real__t Ve_val = apply_limiting(x[EVe_idx], -10000., exc.Ve_Max);
  real__t FEX = get_rectified_regulator(exc.Kc * Ifd / Ve_val);  //Rectified Regulation
  real__t Efd_val = apply_limiting(Ve_val * FEX, -10000., exc.Efd_Max);
  
  /** 量测环节 */
  real__t V_in = get_exc_input(exc.Xc);
  integrate_block(dxdt, V_in, x[EV1_idx], EV1_idx, 1., 1., exc.Tr);

  /** 串联校正环节 */
  real__t val_Vdiff = -x[EV1_idx] + Vm_ref + VS - x[EVf_idx];
  real__t div_Vdiff = -dxdt[EV1_idx] + dxdt[VS_idx] - dxdt[EVf_idx];
  integrate_block(x, dxdt, val_Vdiff, div_Vdiff, x[EV2_idx], EV2_idx, exc.K, exc.K * exc.T1, exc.Kv, exc.T2);
  integrate_block(x, dxdt, EV2_idx, EV3_idx, 1., exc.T3, 1., exc.T4);
  
  /** 并联校正环节 */
  integrate_block(x, dxdt, EVe_idx, EVf_idx, 0., exc.Kf, 1., exc.Tf);
  
  /** 功率放大环节 */
  integrate_block(dxdt, x[EV3_idx], x[EVa_idx], EVa_idx, exc.Ka, 1., exc.Ta);
  
  /** 第二级调节器 */
  real__t Va_val = apply_limiting(x[EVa_idx], exc.Va_Min, exc.Va_Max);
  Va_val -= Efd_val * exc.KH1;
  integrate_block(dxdt, Va_val, x[EVr_idx], EVr_idx, exc.KB, 1., exc.T5);
  
  /** 励磁机 */
  real__t Se = exc.C1 * exp(exc.C2 * Efd); /** 励磁机饱和系数 */
  real__t V_adj = Se * Ve_val + Ifd * exc.Kd;
  real__t Vr_val = apply_limiting(x[EVr_idx], Vm * exc.Vr_Min, Vm * exc.Vr_Max);
  integrate_block(dxdt, Vr_val - V_adj, x[EVe_idx], EVe_idx, 1., exc.Ke, exc.Te);

  print_for_debug_EXC_III_X(x, dxdt, V_in, val_Vdiff, div_Vdiff, Efd_val, Efd_val);
  
    return Efd_val;
}

/** 7型AVR用来模拟由副励磁机向调节器供电的不可控整流交流励磁机励磁系统,主要应用于无刷励磁系统
 * 本模型与3型的区别在于并联校正环节加入点位置不同，3型加入点在串联校正环节之前，而7型加入点在串联校正环节之后
 */
real__t ODE_solver::
process_EPRI_EXC_TYPE_VII(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc) {
  real__t Ve_val = apply_limiting(x[EVe_idx], -10000., exc.Ve_Max);
  real__t FEX = get_rectified_regulator(exc.Kc * Ifd / Ve_val);  //Rectified Regulation
  real__t Efd_val = apply_limiting(Ve_val * FEX, -10000., exc.Efd_Max);
  real__t Se = exc.C1 * exp(exc.C2 * Efd_val); /** 励磁机饱和系数 */
  real__t V_adj = Se * Ve_val + Ifd * exc.Kd;
  
  /** 量测环节 */
  real__t V_in = get_exc_input(exc.Xc);
  integrate_block(dxdt, V_in, x[EV1_idx], EV1_idx, 1., 1., exc.Tr);

  /** 串联校正环节 */
  real__t val_Vdiff = -x[EV1_idx] + Vm_ref + VS;
  real__t div_Vdiff = -dxdt[EV1_idx] + dxdt[VS_idx];
  integrate_block(x, dxdt, val_Vdiff, div_Vdiff, x[EV2_idx], EV2_idx, exc.K, exc.K * exc.T1, exc.Kv, exc.T2);
  integrate_block(x, dxdt, EV2_idx, EV3_idx, 1., exc.T3, 1., exc.T4);

  /** 并联校正环节 */
  real__t Vfe_val = V_adj + exc.Ke * Ve_val;
  real__t Vfe_div = (exc.Ke + Se) * dxdt[EVe_idx] + exc.Kd * dIfd_dt;
  integrate_block(x, dxdt, Vfe_val, Vfe_div, x[EVf_idx], EVf_idx, 0., exc.Kf, 1., exc.Tf);

  /** 功率放大环节 */
  integrate_block(dxdt, x[EV3_idx] - x[EVf_idx], x[EVa_idx], EVa_idx, exc.Ka, 1., exc.Ta);

  /** 第二级调节器 */
  real__t Va_val = apply_limiting(x[EVa_idx], exc.Va_Min, exc.Va_Max);
  Va_val -= Vfe_val * exc.KH1;
  integrate_block(dxdt, Va_val, x[EVr_idx], EVr_idx, exc.KB, 1., exc.T5);

  /** 励磁机 */
  real__t Vr_val = apply_limiting(x[EVr_idx], exc.Vr_Min, exc.Vr_Max);
  integrate_block(dxdt, Vr_val - V_adj, x[EVe_idx], EVe_idx, 1., exc.Ke, exc.Te);

  print_for_debug_EXC_III_X(x, dxdt, V_in, val_Vdiff, div_Vdiff, Vfe_val, Efd_val);
  
  return Efd_val;
}

/** 8型AVR用来模拟由副励磁机向调节器供电的不可控整流交流励磁机励磁系统,应用于有刷励磁系统
 * 本模型与4型的区别在于并联校正环节加入点位置不同，4型加入点在串联校正环节之前，而8型加入点在串联校正环节之后
 */
real__t ODE_solver::
process_EPRI_EXC_TYPE_VIII(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc) {
  real__t Ve_val = apply_limiting(x[EVe_idx], -10000., exc.Ve_Max);
  real__t FEX = get_rectified_regulator(exc.Kc * Ifd / Ve_val);  //Rectified Regulation
  real__t Efd_val = apply_limiting(Ve_val * FEX, -10000., exc.Efd_Max);
  
  /** 量测环节 */
  real__t V_in = get_exc_input(exc.Xc);
  integrate_block(dxdt, V_in, x[EV1_idx], EV1_idx, 1., 1., exc.Tr);

  /** 串联校正环节 */
  real__t val_Vdiff = -x[EV1_idx] + Vm_ref + VS;
  real__t div_Vdiff = -dxdt[EV1_idx] + dxdt[VS_idx];
  integrate_block(x, dxdt, val_Vdiff, div_Vdiff, x[EV2_idx], EV2_idx, exc.K, exc.K * exc.T1, exc.Kv, exc.T2);
  integrate_block(x, dxdt, EV2_idx, EV3_idx, 1., exc.T3, 1., exc.T4);
  
  /** 并联校正环节 */
   integrate_block(x, dxdt, EVe_idx, EVf_idx, 0., exc.Kf, 1., exc.Tf);
   
  /** 功率放大环节 */
  integrate_block(dxdt, x[EV3_idx] - x[EVf_idx], x[EVa_idx], EVa_idx, exc.Ka, 1., exc.Ta);
  
  /** 第二级调节器 */
  real__t Va_val = apply_limiting(x[EVa_idx], exc.Va_Min, exc.Va_Max);
  Va_val -= Efd_val * exc.KH1;
  integrate_block(dxdt, Va_val, x[EVr_idx], EVr_idx, exc.KB, 1., exc.T5);
  
  /** 励磁机 */
  real__t Se = exc.C1 * exp(exc.C2 * Efd_val); /** 励磁机饱和系数 */
  real__t V_adj = Se * Ve_val + Ifd * exc.Kd;
  real__t Vr_val = apply_limiting(x[EVr_idx], exc.Vr_Min, exc.Vr_Max);
  integrate_block(dxdt, Vr_val - V_adj, x[EVe_idx], EVe_idx, 1., exc.Ke, exc.Te);
  
  print_for_debug_EXC_III_X(x, dxdt, V_in, val_Vdiff, div_Vdiff, Efd_val, Efd_val);

  return Efd_val;
}

/** 9型AVR用来模拟没有副励磁机的交流励磁机不可控整流器励磁系统,应用于无刷励磁系统
 * 本模型与5型的区别在于并联校正环节加入点位置不同，5型加入点在串联校正环节之前，而9型加入点在串联校正环节之后
 */
real__t ODE_solver::
process_EPRI_EXC_TYPE_IX(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc) {
  real__t Ve_val = apply_limiting(x[EVe_idx], -10000., exc.Ve_Max);
  real__t FEX = get_rectified_regulator(exc.Kc * Ifd / Ve_val);  //Rectified Regulation
  real__t Efd_val = apply_limiting(Ve_val * FEX, -10000., exc.Efd_Max);
  real__t Se = exc.C1 * exp(exc.C2 * Efd_val); /** 励磁机饱和系数 */
  real__t V_adj = Se * Ve_val + Ifd * exc.Kd;

  /** 量测环节 */
  real__t V_in = get_exc_input(exc.Xc);
  integrate_block(dxdt, V_in, x[EV1_idx], EV1_idx, 1., 1., exc.Tr);

  /** 串联校正环节 */
  real__t val_Vdiff = -x[EV1_idx] + Vm_ref + VS;
  real__t div_Vdiff = -dxdt[EV1_idx] + dxdt[VS_idx];
  integrate_block(x, dxdt, val_Vdiff, div_Vdiff, x[EV2_idx], EV2_idx, exc.K, exc.K * exc.T1, exc.Kv, exc.T2);
  integrate_block(x, dxdt, EV2_idx, EV3_idx, 1., exc.T3, 1., exc.T4);

  /** 并联校正环节 */
  real__t Vfe_val = V_adj + exc.Ke * Ve_val;
  real__t Vfe_div = (exc.Ke + Se) * dxdt[EVe_idx] + exc.Kd * dIfd_dt;
  integrate_block(x, dxdt, Vfe_val, Vfe_div, x[EVf_idx], EVf_idx, 0., exc.Kf, 1., exc.Tf);

  /** 功率放大环节 */
  integrate_block(dxdt, x[EV3_idx] - x[EVf_idx], x[EVa_idx], EVa_idx, exc.Ka, 1., exc.Ta);

  /** 第二级调节器 */
  real__t Va_val = apply_limiting(x[EVa_idx], exc.Va_Min, exc.Va_Max);
  Va_val -= Vfe_val * exc.KH1;
  integrate_block(dxdt, Va_val, x[EVr_idx], EVr_idx, exc.KB, 1., exc.T5);

  /** 励磁机 */
  real__t Vr_val = apply_limiting(x[EVr_idx], Vm * exc.Vr_Min, Vm * exc.Vr_Max);
  integrate_block(dxdt, Vr_val - V_adj, x[EVe_idx], EVe_idx, 1., exc.Ke, exc.Te);

  print_for_debug_EXC_III_X(x, dxdt, V_in, val_Vdiff, div_Vdiff, Vfe_val, Efd_val);

  return Efd_val;
}

/** 10型AVR用来模拟没有副励磁机的交流励磁机不可控整流器励磁系统,应用于有刷励磁系统
 * 本模型与6型的区别在于并联校正环节加入点位置不同，6型加入点在串联校正环节之前，而10型加入点在串联校正环节之后
 */
real__t ODE_solver::
process_EPRI_EXC_TYPE_X(const vector_type& x, vector_type& dxdt, EPRI_EXC_III_TO_X_DATA& exc) {
  real__t Ve_val = apply_limiting(x[EVe_idx], -10000., exc.Ve_Max);
  real__t FEX = get_rectified_regulator(exc.Kc * Ifd / Ve_val);  //Rectified Regulation
  real__t Efd_val = apply_limiting(Ve_val * FEX, -10000., exc.Efd_Max);

  /** 量测环节 */
  real__t V_in = get_exc_input(exc.Xc);
  integrate_block(dxdt, V_in, x[EV1_idx], EV1_idx, 1., 1., exc.Tr);

  /** 串联校正环节 */
  real__t val_Vdiff = -x[EV1_idx] + Vm_ref + VS;
  real__t div_Vdiff = -dxdt[EV1_idx] + dxdt[VS_idx];
  integrate_block(x, dxdt, val_Vdiff, div_Vdiff, x[EV2_idx], EV2_idx, exc.K, exc.K * exc.T1, exc.Kv, exc.T2);
  integrate_block(x, dxdt, EV2_idx, EV3_idx, 1., exc.T3, 1., exc.T4);

  /** 并联校正环节, Efd ~ Ve */
  integrate_block(x, dxdt, EVe_idx, EVf_idx, 0., exc.Kf, 1., exc.Tf);

  /** 功率放大环节 */
  integrate_block(dxdt, x[EV3_idx] - x[EVf_idx], x[EVa_idx], EVa_idx, exc.Ka, 1., exc.Ta);

  /** 第二级调节器 */
  real__t Va_val = apply_limiting(x[EVa_idx], exc.Va_Min, exc.Va_Max);
  Va_val -= Efd_val * exc.KH1;
  integrate_block(dxdt, Va_val, x[EVr_idx], EVr_idx, exc.KB, 1., exc.T5);

  /** 励磁机 */
  real__t Se = exc.C1 * exp(exc.C2 * Efd_val); /** 励磁机饱和系数 */
  real__t V_adj = Se * Ve_val + Ifd * exc.Kd;
  real__t Vr_val = apply_limiting(x[EVr_idx], Vm * exc.Vr_Min, Vm * exc.Vr_Max);
  integrate_block(dxdt, Vr_val - V_adj, x[EVe_idx], EVe_idx, 1., exc.Ke, exc.Te);

  print_for_debug_EXC_III_X(x, dxdt, V_in, val_Vdiff, div_Vdiff, Efd_val, Efd_val);

  return Efd_val;
}

/** 11型AVR用来模拟它励可控整流器励磁系统 */
real__t ODE_solver::
process_EPRI_EXC_TYPE_XI(const vector_type& x, vector_type& dxdt, EPRI_EXC_XI_TO_XII_DATA& exc) {
  /** 量测环节 */
  real__t V_in = get_exc_input(exc.Xc);
  integrate_block(dxdt, V_in, x[EV1_idx], EV1_idx, 1., 1., exc.Tr);

  /** 串联校正环节 */
  real__t val_Vdiff = -x[EV1_idx] + Vm_ref + VS;
  real__t div_Vdiff = -dxdt[EV1_idx] + dxdt[VS_idx];
  integrate_block(x, dxdt, val_Vdiff, div_Vdiff, x[EV2_idx], EV2_idx, exc.K, exc.K * exc.T1, exc.Kv, exc.T2);
  integrate_block(x, dxdt, EV2_idx, EV3_idx, 1., exc.T3, 1., exc.T4);
  
  /** 并联校正环节 */
  integrate_block(x, dxdt, EVa_idx, EVf_idx, 0., exc.Kf, 1., exc.Tf);
  
  /** 功率放大环节 */
  integrate_block(dxdt, x[EV3_idx] - x[EVf_idx], x[EVa_idx], EVa_idx, exc.Ka, 1., exc.Ta);
  real__t Va_val = apply_limiting(x[EVa_idx], exc.Va_Min, exc.Va_Max);
  real__t Vr_val = apply_limiting(Va_val, exc.Vr_Min, exc.Vr_Max);
  
  real__t Efd_val = Vr_val - Ifd * exc.Kc;

  print_for_debug_EXC_III_X(x, dxdt, V_in, val_Vdiff, div_Vdiff, Efd_val, Efd_val);
  
  return Efd_val;
}

/** 12型AVR用来模拟自并励励磁系统（电压源可控整流器励磁系统）。
 * 与11型励磁系统不同的是，可控整流器是由发电机端的整流变压器（又称机端变）供电，没有旋转部件，是静止励磁系统的一种。
 */
real__t ODE_solver::
process_EPRI_EXC_TYPE_XII(const vector_type& x, vector_type& dxdt, EPRI_EXC_XI_TO_XII_DATA& exc) {
  /** 量测环节 */
  real__t V_in = get_exc_input(exc.Xc);
  integrate_block(dxdt, V_in, x[EV1_idx], EV1_idx, 1., 1., exc.Tr);

  /** 串联校正环节 */
  real__t val_Vdiff = -x[EV1_idx] + Vm_ref + VS * (exc.Vs_Pos == 0 ? 1. : 0.);
  real__t div_Vdiff = -dxdt[EV1_idx] + dxdt[VS_idx] * (exc.Vs_Pos == 0 ? 1. : 0.);
  integrate_block(x, dxdt, val_Vdiff, div_Vdiff, x[EV2_idx], EV2_idx, exc.K, exc.K * exc.T1, exc.Kv, exc.T2);
  integrate_block(x, dxdt, EV2_idx, EV3_idx, 1., exc.T3, 1., exc.T4);
  
  /** 并联校正环节 */
  integrate_block(x, dxdt, EVa_idx, EVf_idx, 0., exc.Kf, 1., exc.Tf);
  
  /** 功率放大环节 */
  real__t V_diff = x[EV3_idx] - x[EVf_idx] + VS * (exc.Vs_Pos == 0 ? 0. : 1.);
  integrate_block(dxdt, V_diff, x[EVa_idx], EVa_idx, exc.Ka, 1., exc.Ta);
  real__t Va_val = apply_limiting(x[EVa_idx], exc.Va_Min, exc.Va_Max);
  real__t Vr_val = apply_limiting(Va_val, Vm * exc.Vr_Min, Vm * exc.Vr_Max);
  
  real__t Efd_val = Vr_val - Ifd * exc.Kc;

  print_for_debug_EXC_III_X(x, dxdt, V_in, val_Vdiff, div_Vdiff, Efd_val, Efd0);
  
  return Efd_val;
}

}  // namespace transient_analysis
