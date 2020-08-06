/******************************************************************************
 * Copyright (c): (2019-2020) GEIRINA
 * All rights reserved.
 * Project: Power system transient simulation
 *
 * - Author:      Peng Wei, peng.wei@geirina.net
 * - Created on:  Jan. 29, 2020
 * - Last Update: Jan. 29, 2020
 *
 * - This library is part of the ODE_solver class. All power system stabilizer
 *   (PSS) models should be put in this file.
 *
*******************************************************************************/

#include "ode_solver.cuh"

namespace transient_analysis {

/** PSS 1是一种通用的模型，输入信号根据需要可取转速偏差、功率偏差、端电压偏差。模型采用两级移相结构 */
real__t ODE_solver::
process_EPRI_PSS_TYPE_I(const vector_type& x, vector_type& dxdt, EPRI_PSS_I_DATA &pss) {
  assert(pss.K == 0 || pss.K == 1);
  /*
  real__t VP0 = (pss.Kq1 * (x[omega_idx] - omega_ref) -
                 pss.Kq2 * (Telec - Pe_ref) -
                 pss.Kq3 * (Vm - Vm_ref));
  real__t div_VP0 = pss.Kq1 * dxdt[omega_idx];

  integrate_block(x, dxdt, VP0, div_VP0, x[PP1_idx], PP1_idx, 1. * pss.K, 1., 1., pss.Tq);  //隔直环节
  integrate_block(x, dxdt, PP1_idx, PP2_idx, 1., pss.T1e, 1., pss.T2e);                     //移相环节1
  integrate_block(x, dxdt, PP2_idx, VS_idx, 1., pss.T3e, 1., pss.T4e);                      //移相环节2

#if DEBUG
  cout << "*** PSS debugging data: ***\n";
  cout << "d_PP1_dt = " << dxdt[PP1_idx] << endl;
  cout << "d_PP2_dt = " << dxdt[PP2_idx] << endl;
  cout << "d_VS_dt = " << dxdt[VS_idx] << endl << endl;
#endif
  */
  return apply_limiting(x[VS_idx], pss.VS_Min, pss.VS_Max);
}

/** PSS 2也是一种通用的模型，输入信号根据需要可取转速偏差、电磁功率偏差、机械功率偏差。
 * 与1型PSS相比，2型PSS增加了1级隔直环节和1级移相环节，另外体现了测量环节的惯性延迟作用。
 */
real__t ODE_solver::
process_EPRI_PSS_TYPE_II(const vector_type& x, vector_type& dxdt, EPRI_PSS_II_DATA &pss) {
  /*
  real__t VP0 = (pss.Kw * (x[omega_idx] - omega_ref) -
                 pss.Kp * (Telec - Pe_ref) -
                 pss.Kt * (Vm - Vm_ref));
  real__t div_VP0 = pss.Kw * dxdt[omega_idx];

  integrate_block(x, dxdt, VP0, div_VP0, x[PP1_idx], PP1_idx, 0., pss.TW1, 1., pss.TW1);  //隔直环节1
  integrate_block(x, dxdt, PP1_idx, PP2_idx, 0., pss.TW2, 1., pss.TW2);                   //隔直环节2
  integrate_block(x, dxdt, PP2_idx, PP3_idx, 1., pss.T1, 1., pss.T2);                     //移相环节1
  integrate_block(x, dxdt, PP3_idx, PP4_idx, 1., pss.T3, 1., pss.T4);                     //移相环节2
  integrate_block(x, dxdt, PP4_idx, VS_idx, 1., pss.T5, 1., pss.T6);                      //移相环节3
  
#if DEBUG
  cout << "*** PSS debugging data: ***\n";
  cout << "d_PP1_dt = " << dxdt[PP1_idx] << endl;
  cout << "d_PP2_dt = " << dxdt[PP2_idx] << endl;
  cout << "d_PP3_dt = " << dxdt[PP3_idx] << endl;
  cout << "d_PP4_dt = " << dxdt[PP4_idx] << endl;
  cout << "d_VS_dt = " << dxdt[VS_idx] << endl;
#endif
  */
  return apply_limiting(x[VS_idx], pss.VS_Min, pss.VS_Max);
}

/** PSS 4模型输入信号为delta omega与delta Pe，delta omega输入信号产生等值的机械功率信号，
 * 使总信号对机械功率变化不敏感，在系统振荡频率范围内delta Pe输入维持稳定。
 * Remark: 陷波器环节被省略, 实际上与PSS 6相同
 */
real__t ODE_solver::
process_EPRI_PSS_TYPE_IV(const vector_type& x, vector_type& dxdt, EPRI_PSS_IV_VI_DATA &pss) {
  real__t delta_omega = x[omega_idx] - omega_ref;
  integrate_block(dxdt, delta_omega, x[PP1_idx], PP1_idx, pss.Kw, 1., pss.Trw);  //转速测量
  integrate_block(x, dxdt, PP1_idx, PP2_idx, 0., pss.T5, 1., pss.T6);            //转速隔直1
  integrate_block(x, dxdt, PP2_idx, PP3_idx, 0., pss.T7, 1., pss.T7);            //转速隔直2

  real__t delta_Pe = Telec - Pe_ref;
  integrate_block(dxdt, delta_Pe, x[PP4_idx], PP4_idx, pss.Kr, 1., pss.Trp);    //功率测量
  integrate_block(dxdt, x[PP4_idx], x[PP5_idx], PP5_idx, pss.Tw, 1., pss.Tw1);  //功率隔直1
  integrate_block(x, dxdt, PP5_idx, PP6_idx, 0., pss.Tw2, 1., pss.Tw2);         //功率隔直2

  real__t temp_val = pss.Kp * (x[PP3_idx] + (pss.Ks - 1.) * x[PP6_idx]);
  real__t temp_div = pss.Kp * (dxdt[PP3_idx] + (pss.Ks - 1.) * dxdt[PP6_idx]);
  integrate_block(x, dxdt, temp_val, temp_div, x[PP7_idx], PP7_idx, 1., pss.T1, 1., pss.T2);  //移相环节1
  integrate_block(x, dxdt, PP7_idx, PP8_idx, 1., pss.T13, 1., pss.T14);                       //移相环节2
  integrate_block(x, dxdt, PP8_idx, VS_idx, 1., pss.T3, 1., pss.T4);                          //移相环节3

#if DEBUG
  cout << "*** PSS debugging data: ***\n";
  cout << "d_VP1_dt = " << dxdt[PP1_idx] << endl;
  cout << "d_VP2_dt = " << dxdt[PP2_idx] << endl;
  cout << "d_PP3_dt = " << dxdt[PP3_idx] << endl;
  cout << "d_PP4_dt = " << dxdt[PP4_idx] << endl;
  cout << "d_PP5_dt = " << dxdt[PP5_idx] << endl;
  cout << "d_PP6_dt = " << dxdt[PP6_idx] << endl;
  cout << "d_PP7_dt = " << dxdt[PP7_idx] << endl;
  cout << "d_PP8_dt = " << dxdt[PP8_idx] << endl;
  cout << "d_VS_dt = " << dxdt[VS_idx] << endl;
#endif
  
  return apply_limiting(x[VS_idx], pss.VS_Min, pss.VS_Max);
}

/** PSS 5输入信号取为功率偏差信号。该模型包括两级测量环节、三级隔直环节和移相环节。 */
real__t ODE_solver::
process_EPRI_PSS_TYPE_V(const vector_type& x, vector_type& dxdt, EPRI_PSS_V_DATA &pss) {
  real__t delta_Pe = Telec - Pe_ref;
  integrate_block(dxdt, delta_Pe, x[PP1_idx], PP1_idx, 1., 1., pss.T1);        //功率测量1
  integrate_block(dxdt, x[PP1_idx], x[PP2_idx], PP2_idx, 1., 1., pss.T2);      //功率测量2
  integrate_block(x, dxdt, PP2_idx, PP3_idx, 0., pss.T3, 1., pss.T3);          //功率隔直1
  integrate_block(x, dxdt, PP3_idx, PP4_idx, 0., pss.T4, 1., pss.T4);          //功率隔直2
  integrate_block(x, dxdt, PP4_idx, PP5_idx, 0., pss.T5, 1., pss.T5);          //功率隔直3
  integrate_block(dxdt, x[PP5_idx], x[PP6_idx], PP6_idx, pss.K1, 1., pss.T6);  //速度信号放大

  real__t VS_val = (x[PP5_idx] * pss.K2 * pss.a * (1. - pss.p) +
                    x[PP6_idx] * pss.a * (1. - abs(1. - pss.p)));

#if DEBUG
  cout << "*** PSS debugging data: ***\n";
  cout << "d_VP1_dt = " << dxdt[PP1_idx] << endl;
  cout << "d_VP2_dt = " << dxdt[PP2_idx] << endl;
  cout << "d_PP3_dt = " << dxdt[PP3_idx] << endl;
  cout << "d_PP4_dt = " << dxdt[PP4_idx] << endl;
  cout << "d_PP5_dt = " << dxdt[PP5_idx] << endl;
  cout << "d_PP6_dt = " << dxdt[PP6_idx] << endl;
#endif
  
  return apply_limiting(pss.K * VS_val, pss.VS_Min, pss.VS_Max);
}

/** PSS 6型模型结构与PSS 4类似所不同的是，6型PSS模型无陷波器环节 */
real__t ODE_solver::
process_EPRI_PSS_TYPE_VI(const vector_type& x, vector_type& dxdt, EPRI_PSS_IV_VI_DATA &pss) {
  real__t delta_omega = x[omega_idx] - omega_ref;
  integrate_block(dxdt, delta_omega, x[PP1_idx], PP1_idx, pss.Kw, 1., pss.Trw);  //转速测量
  integrate_block(x, dxdt, PP1_idx, PP2_idx, 0., pss.T5, 1., pss.T6);            //转速隔直1
  integrate_block(x, dxdt, PP2_idx, PP3_idx, 0., pss.T7, 1., pss.T7);            //转速隔直2

  real__t delta_Pe = Telec - Pe_ref;
  integrate_block(dxdt, delta_Pe, x[PP4_idx], PP4_idx, pss.Kr, 1., pss.Trp);    //功率测量
  integrate_block(dxdt, x[PP4_idx], x[PP5_idx], PP5_idx, pss.Tw, 1., pss.Tw1);  //功率隔直1
  integrate_block(x, dxdt, PP5_idx, PP6_idx, 0., pss.Tw2, 1., pss.Tw2);         //功率隔直2

  real__t temp_val = pss.Kp * (x[PP3_idx] + (pss.Ks - 1.) * x[PP6_idx]);
  real__t temp_div = pss.Kp * (dxdt[PP3_idx] + (pss.Ks - 1.) * dxdt[PP6_idx]);
  integrate_block(x, dxdt, temp_val, temp_div, x[PP7_idx], PP7_idx, 1., pss.T1, 1., pss.T2);  //移相环节1
  integrate_block(x, dxdt, PP7_idx, PP8_idx, 1., pss.T13, 1., pss.T14);                       //移相环节2
  integrate_block(x, dxdt, PP8_idx, VS_idx, 1., pss.T3, 1., pss.T4);                          //移相环节3

#if DEBUG
  cout << "*** PSS debugging data: ***\n";
  cout << "d_PP1_dt = " << dxdt[PP1_idx] << endl;
  cout << "d_PP2_dt = " << dxdt[PP2_idx] << endl;
  cout << "d_VS_dt = " << dxdt[VS_idx] << endl;
#endif
  
  return apply_limiting(x[VS_idx], pss.VS_Min, pss.VS_Max);
}

/** PSS 8模型的输入信号取自发电机上网变压器的高压侧母线。 */
real__t ODE_solver::
process_EPRI_PSS_TYPE_VIII(const vector_type& x, vector_type& dxdt, EPRI_PSS_VIII_DATA &pss) {
   /*
  real__t V_diff = Vm_ref - Vm;

  integrate_block(dxdt, V_diff, x[PP1_idx], PP1_idx, pss.Kqv, 1., pss.Tqv);  //电压测量
  integrate_block(x, dxdt, PP1_idx, PP2_idx, 1., pss.Tq1p, 1., pss.Tq1);     //移相环节1
  integrate_block(x, dxdt, PP2_idx, PP3_idx, 1., pss.Tq2p, 1., pss.Tq2);     //移相环节2
  integrate_block(x, dxdt, PP3_idx, VS_idx, 1., pss.Tq3p, 1., pss.Tq3);      //移相环节3
  
#if DEBUG
  cout << "*** PSS debugging data: ***\n";
  cout << "d_PP1_dt = " << dxdt[PP1_idx] << endl;
  cout << "d_PP2_dt = " << dxdt[PP2_idx] << endl;
  cout << "d_PP3_dt = " << dxdt[PP3_idx] << endl;
  cout << "d_VS_dt = " << dxdt[VS_idx] << endl;
#endif
  */
  return apply_limiting(x[VS_idx], pss.VS_Min, pss.VS_Max);
}

}  // namespace transient_analysis
