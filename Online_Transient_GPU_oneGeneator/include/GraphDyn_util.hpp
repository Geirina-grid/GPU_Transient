/******************************************************************************
 * Copyright (c): (2019-2020) GEIRINA
 * All rights reserved.
 * Project: Power system transient simulation
 *
 * - Author:      Peng Wei, peng.wei@geirina.net
 * - Created on:  Jan. 29, 2020
 * - Last Update: Jan. 29, 2020
 *
 * - This library contains all the data structures needed in the simulation.
 *
*******************************************************************************/

#ifndef GRAPHDYN_UTIL_HPP_
#define GRAPHDYN_UTIL_HPP_
 
// #include <mpi.h>
#include <boost/numeric/odeint.hpp>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/** ****** Header files needed for NISCLU ****** */
#include <math.h> /* for using atan2 */
#include <sys/time.h>
#include <unistd.h>
#include <algorithm> /* for using std::sort */
#include <iostream>
#include "ExprUtil.hpp"
#include "TextTable.hpp" /* for printing tables */
#include "graphlu.h"   /* for solving sparse linear equations */
#include "graphlu_util.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/set_operations.h>
#include <thrust/execution_policy.h>

using namespace boost::numeric::odeint;
using vector_type = vector<real__t>;
using h2_vector_type = vector< vector<real__t>>;
//using d2_vector_type = thrust::device_vector<thrust::device_vector< value_type >>;
using matrix_type = boost::numeric::ublas::matrix<real__t>;

//typedef vector<double> value_type;
typedef double value_type;
typedef thrust::device_vector< value_type > d_vector_type;
//typedef thrust::device_vector<thrust::device_vector< value_type >> d2_vector_type;
//typedef thrust::device_vector<vector< value_type >> d2_vector_type;

#define PI 3.141592653589793
#define EPS 1.e-14 /* a small number << 1 to justify nonzero numbers */
#define INF 1.e8  /* a large number >> 1 to represent infinity */
#define DEBUG 0    /* DEBUG = 1, debug mode; DEBUG = 0, release mode */

namespace transient_analysis {

enum IDX {
  omega_idx, delta_idx, Eqp_idx, Edp_idx, Eqpp_idx, Edpp_idx,                         // GEN
  delta_mu_idx, GP1_idx, GP2_idx, GP3_idx, GPm_idx, GPf_idx, hp_idx, ip_idx, lp_idx,  // governor
  EV1_idx, EV2_idx, EV3_idx, EVf_idx, EVa_idx, EVr_idx, EVe_idx,                      // excitor
  PP1_idx, PP2_idx, PP3_idx, PP4_idx, PP5_idx, PP6_idx, PP7_idx, PP8_idx, VS_idx,     // PSS
  mu_output_idx, PT_output_idx, Efd_output_idx, VS_output_idx                         // output
};

struct BUS {
  string bus_name;
  int bus_id;
  int bus_type;
  real__t IBase, ZBase;
  real__t Vm, Va;
  real__t Pgen, Qgen;
  real__t Pload, Qload;
  real__t Gshunt, Bshunt, invX;
};

struct GENERATOR {
  string bus_name;
  int bus_id;
  int Gen_Model, Gen_Par; // model number, parameter set number
  int AVR_Model, AVR_Par; // model number, parameter set number
  int GOV_Model, GOV_Par; // model number, parameter set number
  int PSS_Model, PSS_Par; // model number, parameter set number
  real__t Rate_MVA, Rate_MW;
  real__t Xdp, Xdpp, X2, TJ;
  real__t omega_ref, freq_ref;
  real__t Vt_ref, Pe_ref;
  real__t Efd0, mu0;
  real__t Ep0, Epp0;
};

struct LINE {
  string from, to;
  real__t G, B, B_half, Gm, Bm, invX, tap;
  int tap_count;
};

struct GENROU_IEEE_DATA {
  uint__t gen_id;
  uint__t bus_id;
  real__t Xl, Ra;
  real__t Xd, Xdp, Xdpp, Xq, Xqp, Xqpp;
  real__t Td0p, Td0pp, Tq0p, Tq0pp, TJ;
  real__t D;
  real__t S1, S12;
  uint__t gen_type;
  uint__t order;
};

struct EPRI_GEN_DATA {
  string bus_name;
  uint__t gen_id;
  uint__t bus_id;
  real__t Xd, Xdp, Xdpp, Xq, Xqp, Xqpp;
  real__t X2, Ra;
  real__t Td0p, Td0pp, Tq0p, Tq0pp, TJ;
  real__t a, b, n;
  real__t D;
};

struct EPRI_GEN_DATA_D {
  string bus_name;
  uint__t gen_id;
  uint__t bus_id;
  real__t Xd, Xdp, Xdpp, Xq, Xqp, Xqpp;
  real__t X2, Ra;
  real__t Td0p, Td0pp, Tq0p, Tq0pp, TJ;
  real__t a, b, n;
  real__t D;
};

struct EXC_IEEE_I_DATA {
  uint__t exc_type;
  real__t TR;
  real__t KA, TA;
  real__t VR_Max, VR_Min;
  real__t KE, TE;
  real__t KF, TF;
  real__t a, b;
};

struct EPRI_EXC_I_DATA {
  real__t Kr;                //量测环节放大倍数
  real__t Ka;                //放大环节（可控硅传递函数）放大倍数
  real__t Kf;                //反馈环节放大倍数
  real__t Efd_Max, Efd_Min;  //励磁电压上下限，标幺值(p.u.)
  real__t Tr;                //量测环节时间常数，单位为秒(s)
  real__t Ta;                //放大环节（可控硅传递函数）时间常数，单位为秒(s)
  real__t Tf;                //反馈环节时间常数，单位为秒(s)
  real__t Te;                //励磁机时间常数，单位为秒(s)
};

struct EPRI_EXC_II_DATA {
  real__t Kr;                //量测环节放大倍数
  real__t Ka;                //放大环节放大倍数
  real__t K2;                //变换环节类型的参数。K2=0为比例积分，K2=1为移相环节
  real__t Vta;               //强励电压达Efd_Max时的端电压，标幺值,以初始稳态电压Vt0为基准
  real__t Vtb;               //强励电压达Efd_Min时的端电压，标幺值,以初始稳态电压Vt0为基准
  real__t Kpt;               //自励电压系数
  real__t Kit;               //自励电流系数
  real__t Ke;                //换弧压降系数
  real__t Efd_Max, Efd_Min;  //励磁电压上下限，标幺值(p.u.)
  real__t Tr;                //量测环节时间常数，单位为秒(s)
  real__t Ta;                //放大环节时间常数，单位为秒(s)
  real__t T1, T2, T3, T4;    //中间环节时间常数，单位为秒(s)
};

struct EPRI_EXC_III_TO_X_DATA {
  real__t Xc;                      //调差电抗，标幺值(p.u.)
  real__t Tr;                      //量测环节时间常数，单位为秒(s)
  real__t K;                       //串联校正环节的直流增益
  real__t Kv;                      //积分校正选择因子，Kv=0时为纯积分型校正，Kv=1时为比例积分型校正
  real__t T1, T2, T3, T4;          //串联校正环节时间常数，单位为秒(s)，一般有T4<T3<T1<T2
  real__t Ka, Ta, Va_Max, Va_Min;  //功率放大环节
  real__t Kf, Tf;                  //并联校正环节
  real__t KH1;                     //用来调节时间常数补偿度的比例反馈系数
  real__t KB, T5;                  //第二级调节器增益,时间常数
  real__t Vr_Max, Vr_Min;          //调节器输出上下限
  real__t Ke, Te, Ve_Max;  //励磁机
  real__t C1, C2; //用于求取励磁机饱和系数Se的参数
  real__t Kc;                      //与换相电抗相关的整流器负荷系数
  real__t Kd;                      //交流励磁机负载电流电枢反应的去磁系数，是交流励磁机电抗的函数
  real__t Efd_Max;                 //励磁电压上限，标幺值(p.u.)
  
  real__t MQL;                     //瞬时强励限制启用标志，MQL=1时启用瞬时强励限制，MQL=0时不启用
  real__t KL1;                     //励磁机励磁电流限制增益（瞬时强励限制），MQL=0 时不需要
  real__t VL1R;                    //励磁机电流限制（瞬时强励限制），标幺值(p.u.)，MQL=0 时不需要
  real__t MGL;                     //过励限制启用标志，MGL=1 时启用过励限制，MGL=0 时不启用
  real__t Ifd_Inf, Vfd_Inf;        //发电机磁场长期允许电流/电压，标幺值(p.u.)，MGL=0 时不需要
  real__t B;                       //过励发热允许值，MGL=0 时不需要
  real__t C;                       //过励恢复系数，标幺值(p.u.)，MGL=0 时不需要
  real__t T;                       //允许的过励时间，单位为秒(s)，MGL=0 时不需要
  real__t KL2;                     //过励限制回路增益，MGL=0 时不需要
  real__t TL1, TL2;                //过励限制回路时间常数，单位为秒(s)，MGL=0 时不需要
  real__t TYPEDL;                  //低励限制类型：TYPEDL=0 时不启用低励限制（缺省）
                                   //TYPEDL=1 时为直线型低励限制，TYPEDL=2时为圆周型低励限制
  real__t P1, Q1, P2, Q2;          //低励限制曲线上第一(二)点有功(无功)功率，标幺值，TYPEDL=0不需要
  real__t KH2;                     //低励限制回路增益，TYPEDL=0 不需要
  real__t TH1, TH2;                //低励限制回路时间常数，单位为秒(s)，TYPEDL=0 不需要
};

struct EPRI_EXC_XI_TO_XII_DATA {
  real__t Xc;                      //调差电抗，标幺值(p.u.)
  real__t Tr;                      //量测环节时间常数，单位为秒(s)
  real__t K;                       //串联校正环节的直流增益
  real__t Kv;                      //积分校正选择因子，Kv=0时为纯积分型校正，Kv=1时为比例积分型校正
  real__t T1, T2, T3, T4;          //串联校正环节时间常数，单位为秒(s)，一般有T4<T3<T1<T2
  real__t Ka, Ta, Va_Max, Va_Min;  //功率放大环节
  real__t Kf, Tf;                  //并联校正环节
  real__t Vr_Max, Vr_Min;          //调节器输出上下限
  real__t Kc;                      //与换相电抗相关的整流器负荷系数
  uint__t Vs_Pos;                  //Vs输入位置，可以选择Vs输入在串联校正环节前或后 (仅EXC 12中使用)
  
  real__t MGL;                     //过励限制启用标志，MGL=1 时启用过励限制，MGL=0 时不启用
  real__t Ifd_Inf, Vfd_Inf;        //发电机磁场长期允许电流/电压，标幺值(p.u.)，MGL=0 时不需要
  real__t B;                       //过励发热允许值，MGL=0 时不需要
  real__t C;                       //过励恢复系数，标幺值(p.u.)，MGL=0 时不需要
  real__t T;                       //允许的过励时间，单位为秒(s)，MGL=0 时不需要
  real__t KL2;                     //过励限制回路增益，MGL=0 时不需要
  real__t TL1, TL2;                //过励限制回路时间常数，单位为秒(s)，MGL=0 时不需要
  real__t TYPEDL;                  //低励限制类型：TYPEDL=0 时不启用低励限制（缺省）
                                   //TYPEDL=1 时为直线型低励限制，TYPEDL=2时为圆周型低励限制
  real__t P1, Q1, P2, Q2;          //低励限制曲线上第一(二)点有功(无功)功率，标幺值，TYPEDL=0不需要
  real__t KH2;                     //低励限制回路增益，TYPEDL=0 不需要
  real__t TH1, TH2;                //低励限制回路时间常数，单位为秒(s)，TYPEDL=0 不需要
};

struct PID_DATA {
  real__t Kp, Kd, Ki;        //比例环节增益,积分环节增益,微分环节增益
  real__t I_Max, I_Min;      //积分环节上/下限
  real__t PID_Max, PID_Min;  //PID 调节器输出上/下限
};

struct STEAM_DATA {
  real__t Fhp, Fip, Flp;  //高/中/低压缸功率比例
  real__t lambda;         //高压缸功率自然过调系数
  real__t Tch, Trh, Tco;  //蒸汽容积时间常数，再热蒸汽容积时间常数，交叉管时间常数，单位为秒(s)
};

struct EHS_DATA {
  PID_DATA pid;                 //电液转换PID模块
  real__t VEL_Open, VEL_Close;  //油动机最大开启/关闭速度,标幺值
  real__t TO, TC;               //油动机开启/关闭时间常数
  real__t T2;                   //反馈环节滤波时间常数
  real__t P_Max, P_Min;         //汽门/导水页开度上/下限
  int is_open;
};

struct EPRI_GOV_I_DATA {
  uint__t gov_type;              // 1 - steam, 2 - hydro
  real__t K_delta;               //量测环节放大倍数
  real__t Ki;                    //硬负反馈放大倍数,对于汽轮机,Ki=1；
  real__t Kbeta;                 //软负反馈放大倍数,对于汽轮机,Kbeta=0
  real__t dead_band_tol;         //调速器死区，标幺值(p.u.)
  real__t sigma_Max, sigma_Min;  //配压阀行程上/下限，标幺值(p.u.)
  real__t mu_Max, mu_Min;        //导水叶(汽门)开度上/下限，标幺值(p.u.)
  real__t alpha;                 //汽轮机过热系数。若无中间过热，alpha=1；对于水轮机，alpha=1
  real__t Ti;                    //水轮机软反馈时间常数，单位为秒(s)。对于汽轮机无此参数，可取Ti=100000
  real__t TS;                    //伺服机构时间常数，单位为秒(s)
  real__t TW;                    //水锤效应时间常数，单位为秒(s);若无水锤效应，TW=0
  real__t T0;                    //蒸汽容积时间常数，单位为秒(s);对于水轮机，T0=TW/2
  real__t TRH;                   //汽轮机中间过热时间常数，单位为秒(s)。
                                 //对无中间过热的汽轮机和水轮机，无TRH 参数，可令TRH=100000和0
  real__t mu_zero, mu_zero_base;
  real__t KmH;
};

struct EPRI_GOV_II_DATA {
  real__t dead_band_tol;        //调节死区，计算时取其值的1/2作为正负死区
  real__t K;                    //转速放大倍数
  real__t Tr;                   //调速器滑阀组时间常数
  real__t Tb;                   //中间滑阀组时间常数
  real__t VEL_Open, VEL_Close;  //油动机最大开启/关闭速度,标幺值
  real__t TO, TC;               //油动机开启/关闭时间常数
  real__t U_Max, U_Min;         //原动机最大/小汽门开度，标幺值
  STEAM_DATA steam;             //汽轮机模型
};

struct EPRI_GOV_III_DATA{
  real__t dead_band_tol;    //调节死区，计算时取其值的1/2作为正负死区
  real__t K1, T1;           //转速变换
  PID_DATA pid_load;        //负荷控制器PID
  PID_DATA pid_pressure;    //调节级压力控制器PID
  real__t CON_Max, CON_Min; //电液调节系统输出上下限
  EHS_DATA ehs;             //电液伺服系统页
  STEAM_DATA steam;         //汽轮机模型
  int load_feed_control, load_pid_control, pressure_pid_control;
};

struct EPRI_GOV_IV_DATA {
  real__t dead_band_tol;  //调节死区，计算时取其值的1/2作为正负死区
  real__t K;              //转速放大倍数
  real__t T1;             //转速测量时间常数
  real__t K2;             //负荷前馈系数
  PID_DATA pid;           //负荷控制器PID or 调节级压力控制器PID
  EHS_DATA ehs;           //电液伺服系统页
  STEAM_DATA steam;       //汽轮机模型
  int control_option;     //控制方式选择: 1,调节器压力反馈控制、2,DEH开环控制、3,负荷反馈控制
};

struct EPRI_GOV_V_DATA {
  real__t dead_band_tol;  //调节死区，计算时取其值的1/2作为正负死区
  real__t K1;             //转速放大倍数
  real__t T1;             //转速测量时间常数，单位为秒
  real__t K2;             //主汽压力偏差放大倍数
  PID_DATA pid;           //负荷控制器PID
  EHS_DATA ehs;           //电液伺服系统页
  STEAM_DATA steam;       //汽轮机模型
  int control_speed;
  int control_option;     //控制方式选择: 0,CCS自动控制、1,负荷开环控制、2,带主汽压力修正负荷控制
};

/** 7型GOV -- 水轮机调速器模型，包括调节系统模型、液压系统模型以及水轮机模型 */
struct EPRI_GOV_VII_DATA {
  int auto_manual_mode;                    //控制方式: 0-自动，1-手动
  int additional_control;                  //附加调节方式: 0-功率调节， 1-开度调节，2-无调节
  int Y_control_input_option;              //开度调节输入信号选择：0-开度，1-Y PID输出值
  int additional_control_input_position;   //功率/开度偏差信号接入点选择：0-PID前，1-PID内积分项前
  real__t YC;                              //手动方式开度给定值，当控制方式选为手动调节时有效
  real__t T1, TR1;                         //转速延迟时间,转速测量环节时间常数，单位为秒(s)
  real__t dead_band1;                      //转速偏差死区，标幺值(p.u.，相对于额定频率)
  real__t dead_band_Max1, dead_band_Min1;  //频率调节限幅最大值/最小值
  real__t Kw;                              //转速偏差放大倍数

  real__t T2, TR2;                         //功率延迟时间,转速测量环节时间常数，单位为秒(s)
  real__t dead_band2;                      //功率偏差死区，标幺值(p.u.，相对于额定功率)
  real__t dead_band_Max2, dead_band_Min2;  //功率调节限幅最大值/最小值
  real__t Ep;                              //功率偏差放大倍数

  real__t T3, TR3;                         //开度延迟时间,转速测量环节时间常数，单位为秒(s)
  real__t dead_band3;                      //开度偏差死区，标幺值(p.u.，相对于额定功率)
  real__t dead_band_Max3, dead_band_Min3;  //开度调节限幅最大值/最小值
  real__t Bp;                              //开度偏差放大倍数

  /** 调节系统PID */
  real__t Kp, Kd, Ki;          //比例环节增益,积分环节增益,微分环节增益
  real__t Td;                  //微分时间常数
  real__t I_Max, I_Min;        //积分环节上/下限
  real__t PID_Max, PID_Min;    //PID 调节器输出上/下限
  real__t PRO_Max, PRO_Min;    //PID 调节器输出上/下限
  real__t Ratelimp, Ratelimn;  //开度调节速率上/下限

  EHS_DATA ehs;  //液压系统（含执行机构）
  real__t T4;    //液压系统功率输出纯延迟时间，单位为秒(s)

  /** 水轮机相关参数 */
  int hydro_type;                            //0-不考虑水击模型，1-刚性水击模型，2-近似弹性水击模型，3-弹性水击模型
  real__t Tf;                                //油动机反馈时间
  real__t Tw;                                //水轮机水锤效应时间常数
  real__t Tr;                                //管道反射时间，也称弹性水击时间常数
  real__t At;                                //水轮机增益，标幺值
  real__t YNL;                               //空载开度
  real__t dead_band_Max_w, dead_band_Min_w;  //水轮机正向/反向死区
};

/** 8型调速器 -- 水轮机模型,仅实现 A：开度模式 */
struct EPRI_GOV_VIII_DATA {
  /** A模式参数*/
  int additional_control;                  //附加调节方式: 1-功率调节， 2-开度调节
  int Y_control_input_option;              //开度调节输入信号选择：1-开度，2-Y PID输出值
  int additional_control_input_position;   //开度偏差信号接入点选择：1-PID前，2-PID内积分项前
  real__t T1, TR1;                         //转速延迟时间,转速测量环节时间常数，单位为秒(s)
  real__t dead_band_p, dead_band_n;        //转速偏差死区，标幺值(p.u.，相对于额定频率)
  real__t dead_band_Max1, dead_band_Min1;  //频率调节限幅最大值/最小值
  real__t Kw;                              //转速偏差放大倍数
  real__t Rti0, Rtd0;                      //正向/负向速率限制
  real__t T3, TR3;                         //开度延迟时间,转速测量环节时间常数，单位为秒(s)
  real__t dead_band3;                      //开度偏差死区，标幺值(p.u.，相对于额定功率)
  real__t dead_band_Max3, dead_band_Min3;  //开度调节限幅最大值/最小值
  real__t Bp;                              //开度偏差放大倍数
  /** 调节系统PID */
  real__t Kp2, Kd2, Ki2;        //比例环节增益,积分环节增益,微分环节增益
  real__t Td2;                //微分时间常数
  real__t I_Max2, I_Min2;      //积分环节上/下限
  real__t PID_Max2, PID_Min2;  //PID 调节器输出上/下限
  
  real__t Rti1, Rtd1;                      //正向/负向速率限制

  //  /** B,功率模式参数*/
  //  real__t T2, TR2;                         //功率延迟时间,转速测量环节时间常数，单位为秒(s)
  //  real__t dead_band2;                      //功率偏差死区，标幺值(p.u.，相对于额定功率)
  //  real__t dead_band_Max2, dead_band_Min2;  //功率调节限幅最大值/最小值
  //  real__t Ep;                              //功率偏差放大倍数
  //  real__t Rti2, Rtd2; //正向/负向速率限制
  //  PID_DATA pid_B;      //控制器PID

  EHS_DATA ehs;  //液压系统（含执行机构）
  real__t T4;    //液压系统功率输出纯延迟时间，单位为秒(s)
  real__t Tf;    //油动机反馈时间

  /** 水轮机相关参数 */
  real__t Tw;  //水轮机水锤效应时间常数
};

struct EPRI_GOV_IX_DATA {
  real__t dead_band_tol;              //调节死区，计算时取其值的1/2作为正负死区
  real__t K;                          //转速放大倍数
  real__t T1;                         //转速测量时间常数
  real__t K2;                         //负荷前馈系数
  real__t Tr_pm;                      //机械功率量测时间常数，单位为秒
  real__t Tr_pe;                      //电磁功率量测时间常数，单位为秒(s)
  real__t Tdelay1, Tdelay2, Tdelay3;  //延时,取值0~5s
  real__t Rr, Rf;                     //上升/下降速率，功率标幺值/s
  PID_DATA pid;                       //负荷控制器PID or 调节级压力控制器PID
  EHS_DATA ehs;                       //电液伺服系统页
  STEAM_DATA steam;                   //汽轮机模型
  int control_option;                 //控制方式选择: 0,调节器压力反馈控制、1,DEH开环控制、2,负荷反馈控制
};

struct EPRI_PSS_I_DATA {
  uint__t pss_type;
  real__t Kq1, Kq2, Kq3;       //转速偏差放大倍数,电磁功率偏差放大倍数,电压偏差放大倍数
  real__t K;                   //类型系数: K=0为隔直环节（惯性微分环节）；K=1为移相环节
  real__t Tq;                  //隔直环节时间常数
  real__t T1e, T2e, T3e, T4e;  //移相环节时间常数
  real__t VS_Max, VS_Min;      //PSS输出上/下限
};

struct EPRI_PSS_II_DATA {
  uint__t pss_type;
  real__t Kw, Kp, Kt;              //转速偏差放大倍数,电磁功率偏差放大倍数,电压偏差放大倍数
  real__t TW1, TW2;                //隔直环节时间常数,单位为秒(s)，通常取较大值(~10s)
  real__t T1, T2, T3, T4, T5, T6;  //移相环节时间常数
  real__t VS_Max, VS_Min;          //PSS输出上/下限
};

struct EPRI_PSS_IV_VI_DATA {
  uint__t pss_type;
  real__t Kw, Trw;                   //转速偏差放大倍数,转速测量时间常数
  real__t T5, T6, T7;                //转速隔直环节时间常数,单位为秒(s)
  real__t Kr, Trp;                   //功率偏差放大倍数,功率测量时间常数
  real__t Ks;                        //功率偏差补偿系数
  real__t Tw, Tw1, Tw2;              //功率隔直环节时间常数，单位为秒(s)
  real__t T9, T10, T12;              //陷波器时间常数，单位为秒(s)
  real__t T1, T2, T3, T4, T13, T14;  //移相环节时间常数，单位为秒(s)
  real__t Kp;                        //PSS比例放大倍数
  real__t VS_Max, VS_Min;            //PSS输出上/下限
};

struct EPRI_PSS_V_DATA {
  uint__t pss_type;
  real__t T1, T2;          //测量环节时间常数
  real__t T3, T4, T5;      //隔直环节时间常数,单位为秒(s)
  real__t K1;              //速度信号放大倍数
  real__t T6;              //速度信号时间常数
  real__t K2;              //加速度信号放大倍数
  real__t a;               //PSS 可调幅值
  real__t p;               //PSS 可调相位
  real__t K;               //系数，一般为1
  real__t VS_Max, VS_Min;  //PSS输出上/下限
};

struct EPRI_PSS_VI_DATA {
  uint__t pss_type;
  real__t Kw, Trw;                   //转速偏差放大倍数,转速测量时间常数
  real__t T5, T6, T7;                //转速隔直环节时间常数,单位为秒(s)
  real__t Kr, Trp;                   //功率偏差放大倍数,功率测量时间常数
  real__t Ks;                        //功率偏差补偿系数
  real__t Tw, Tw1, Tw2;              //功率隔直环节时间常数，单位为秒(s)
  real__t T1, T2, T3, T4, T13, T14;  //移相环节时间常数，单位为秒(s)
  real__t Kp;                        //PSS比例放大倍数
  real__t VS_Max, VS_Min;            //PSS输出上/下限
};

struct EPRI_PSS_VIII_DATA {
  uint__t pss_type;
  real__t Kqv, Tqv;                         //电压控制放大倍数,电压测量时间常数
  real__t Tq1, Tq1p, Tq2, Tq2p, Tq3, Tq3p;  //移相环节时间常数，单位为秒(s)
  real__t VS_Max, VS_Min;                   //PSS输出上/下限
};

struct EXC_DATA {
  uint__t exc_type;
  real__t KR, TR;
  uint__t K2 /* = 0 or 1 */;
  real__t T1, T2, T3, T4;
  real__t KA, TA;
  real__t K_pt, K_it, K_e;
  real__t KF, TF;
  real__t Te;
  real__t Efd_Max, Efd_Min;
  real__t Vt_Max, Vt_Min;
};

struct GOV_I_DATA {
  uint__t gov_type /* 1 - steam, 2 - hydro */;
  real__t K_delta;
  real__t TS;
  real__t dead_band_tol;
  real__t sigma_Max, sigma_Min;
  real__t mu_zero, mu_zero_base;
  real__t mu_Max, mu_Min;
  real__t T0, TW, alpha, TRH;
  real__t KmH;
  real__t Ki, Kbeta, Ti;
};

struct GOV_BPAGG_DATA {
  uint__t gov_type /* 1 - steam, 2 - hydro */;
  real__t R;
  real__t F;
  real__t Pmax;
  real__t T1, T2, T3, T4, T5;
};

struct PSS_I_DATA {
  uint__t pss_type;
  real__t Kq1, Kq2, Kq3;
  uint__t K /* = 0 or 1 */;
  real__t Tq;
  real__t KW, TW;
  real__t T1e, T2e, T3e, T4e;
  real__t VS_Max, VS_Min;
};

struct GENERATOR_DATA_PACK {
  uint__t EXC_type, GOV_type, PSS_type, GEN_type;
  
  GENROU_IEEE_DATA genrou;
  EPRI_GEN_DATA gen;
  
  EXC_IEEE_I_DATA exc_ieee_1;
  EPRI_EXC_I_DATA exc_1;
  EPRI_EXC_II_DATA exc_2;
  EPRI_EXC_III_TO_X_DATA exc_3_10;
  EPRI_EXC_XI_TO_XII_DATA exc_11_12;
  
  GOV_BPAGG_DATA gov_bpagg;
  EPRI_GOV_I_DATA gov_1;
  EPRI_GOV_II_DATA gov_2;
  EPRI_GOV_III_DATA gov_3;
  EPRI_GOV_IV_DATA gov_4;
  EPRI_GOV_V_DATA gov_5;
  EPRI_GOV_VII_DATA gov_7;
  EPRI_GOV_VIII_DATA gov_8;
  EPRI_GOV_IX_DATA gov_9;
  
  EPRI_PSS_I_DATA pss_1;
  EPRI_PSS_II_DATA pss_2;
  EPRI_PSS_IV_VI_DATA pss_4_6;
  EPRI_PSS_V_DATA pss_5;
  EPRI_PSS_VIII_DATA pss_8;

  real__t omega_ref;
  real__t Vt_ref;
  real__t Pe_ref;
  real__t freq_ref;
  real__t Efd0;
  real__t mu0;
  real__t Rate_MW;
};

struct GENERATOR_DATA_PACK_D {
  value_type EXC_type, GOV_type, PSS_type, GEN_type;

  //GENROU_IEEE_DATA genrou;
  EPRI_GEN_DATA gen;

  //EXC_IEEE_I_DATA exc_ieee_1;
  //EPRI_EXC_I_DATA exc_1;
  //EPRI_EXC_II_DATA exc_2;
  //EPRI_EXC_III_TO_X_DATA exc_3_10;
  //EPRI_EXC_XI_TO_XII_DATA exc_11_12;

  //GOV_BPAGG_DATA gov_bpagg;
  //EPRI_GOV_I_DATA gov_1;
  //EPRI_GOV_II_DATA gov_2;
  //EPRI_GOV_III_DATA gov_3;
  //EPRI_GOV_IV_DATA gov_4;
  //EPRI_GOV_V_DATA gov_5;
  //EPRI_GOV_VII_DATA gov_7;
  //EPRI_GOV_VIII_DATA gov_8;
  //EPRI_GOV_IX_DATA gov_9;

  //EPRI_PSS_I_DATA pss_1;
  //EPRI_PSS_II_DATA pss_2;
  //EPRI_PSS_IV_VI_DATA pss_4_6;
  //EPRI_PSS_V_DATA pss_5;
  //EPRI_PSS_VIII_DATA pss_8;

  value_type omega_ref;
  value_type Vt_ref;
  value_type Pe_ref;
  value_type freq_ref;
  value_type Efd0;
  value_type mu0;
  value_type Rate_MW;
};

template <typename T>
void sortrows(std::vector<std::vector<T>>& matrix);

template <typename T>
void print_matrix(const vector<vector<T>>& mtx, const string& info);

template <typename T>
void print_vector(const vector<T>& vec, const string& info);

template <typename T>
void print_array(T* vec, const int size, const string& info);

void printLogo();
void print_line();

}  // namespace transient_analysis

#endif
