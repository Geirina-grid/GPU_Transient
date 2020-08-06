/******************************************************************************
 * Copyright (c): (2019-2020) GEIRINA
 * All rights reserved.
 * Project: Power system transient simulation
 *
 * - Author:      Peng Wei, peng.wei@geirina.net
 * - Created on:  Jan. 29, 2020
 * - Last Update: Jan. 29, 2020
 *
 * - This library is part of the Transient_Solver class. The initial conditions
 *   for all the variables in the generator differential equations are calculated
 *   here.
 *
*******************************************************************************/

#include "transient.hpp"

namespace transient_analysis {

real__t Transient_Solver::
get_exc_input(real__t Vd, real__t Vq, real__t Id, real__t Iq, real__t Xc) {
  real__t rad_diff = -atan2(Vq, Vd) + atan2(Iq, Id); // Voltage angle is the reference
  real__t re_val = Vd + Id * Xc * sin(rad_diff);
  real__t im_val = Vq + Iq * Xc * sin(rad_diff);
  return sqrt(re_val * re_val + im_val * im_val);
}

void Transient_Solver::compute_ODE_initial_values() {
  for (auto& g_hldr : generators) {
    auto & bus_name=g_hldr.first;
    auto & gen=g_hldr.second;
#if DEBUG
    printf("Transient simulation computing generator ODE intial values for bus %s:\n", bus_name.c_str());
#endif
    uint__t bus_id = bus_name_to_id[bus_name] - 1;
    assert(bus_id >= 0 && bus_id < nBuses);
    
    BUS &bus = buses[bus_name];
    EPRI_GEN_DATA &gen_data = all_gen[gen.Gen_Par];
        
    real__t Xd   = (gen.Gen_Par == 0) ? gen.Xdp  : gen_data.Xd;
    real__t Xdp  = (gen.Gen_Par == 0) ? gen.Xdp  : gen_data.Xdp;
    real__t Xdpp = (gen.Gen_Par == 0) ? gen.Xdpp : gen_data.Xdpp;
    real__t Xq   = (gen.Gen_Par == 0) ? gen.Xdp  : gen_data.Xq;
    real__t Xqp  = (gen.Gen_Par == 0) ? gen.Xdp  : gen_data.Xqp;
    real__t Xqpp = (gen.Gen_Par == 0) ? gen.Xdpp : gen_data.Xqpp;
    
    real__t Vm = bus.Vm;
    real__t Va = bus.Va * PI / 180;  // degree to rad
    real__t Pg = bus.Pgen;  // gen.Rate_MVA;
    real__t Qg = bus.Qgen;  // gen.Rate_MW;
    
    real__t Vx  = Vm * cos(Va);
    real__t Vy  = Vm * sin(Va);
    real__t Ix  = (Pg * cos(Va) + Qg * sin(Va)) / Vm;
    real__t Iy  = (Pg * sin(Va) - Qg * cos(Va)) / Vm;
    real__t EQx = Vx + gen_data.Ra * Ix - Xq * Iy;
    real__t EQy = Vy + gen_data.Ra * Iy + Xq * Ix;
    
    real__t delta0 = atan2(EQy, EQx);
    
    real__t Vd = Vx * sin(delta0) - Vy * cos(delta0);
    real__t Vq = Vx * cos(delta0) + Vy * sin(delta0);
    real__t Id = Ix * sin(delta0) - Iy * cos(delta0);
    real__t Iq = Ix * cos(delta0) + Iy * sin(delta0);
    
    real__t Pmech0 = Pg;
    
    real__t Efd0  = sqrt(EQx * EQx + EQy * EQy) + (Xd - Xq) * Id;
    real__t Edp0  = Vd + gen_data.Ra * Id - Xqp * Iq;
    real__t Eqp0  = Vq + gen_data.Ra * Iq + Xdp * Id;
    real__t Edpp0 = Vd + gen_data.Ra * Id - Xqpp * Iq;
    real__t Eqpp0 = Vq + gen_data.Ra * Iq + Xdpp * Id;
    
//    // for comparison:
//    real__t Edp0_new = (Xq - Xqp) * Iq;
//    real__t Eqp0_new = Efd0 - (Xd - Xdp) *Id;
//    cout  << "two method difference: " << Edp0 - Edp0_new << ", " << Eqp0 - Eqp0_new << endl;
    
    real__t Ifd0   = gen_data.a * Eqp0 + gen_data.b * pow(Eqp0, gen_data.n);
    real__t omega0 = 1.0;
    
    gen.omega_ref = omega0;
    gen.freq_ref  = 50.;
    gen.Vt_ref    = Vm;
    gen.Pe_ref    = Pmech0;
    gen.Efd0      = Efd0;
    gen.mu0       = 0;
    gen.Ep0       = sqrt(Edp0 * Edp0 + Eqp0 * Eqp0);
    gen.Epp0      = sqrt(Edpp0 * Edpp0 + Eqpp0 * Eqpp0);
    
    /** GEN */
    gen_solution[bus_name][omega_idx] = omega0;
    gen_solution[bus_name][delta_idx] = delta0;
    gen_solution[bus_name][Edpp_idx]  = Edpp0;
    gen_solution[bus_name][Eqpp_idx]  = Eqpp0;
    gen_solution[bus_name][Edp_idx]   = Edp0;
    gen_solution[bus_name][Eqp_idx]   = Eqp0;
    
    /** GOV */
    switch (gen.GOV_Model) {
      case 0: break;
      case 1: {
        real__t alpha = all_gov_1[gen.GOV_Par].alpha;
        real__t KmH   = gen.Rate_MW / 100.;
        gen_solution[bus_name][delta_mu_idx] = 0;
        gen_solution[bus_name][GPf_idx]      = 0;
        gen_solution[bus_name][GP1_idx]      = Pmech0;
        gen_solution[bus_name][GP2_idx]      = Pmech0 * (1. - alpha);
        gen.mu0 = Pmech0 / KmH;
        break;
      }
      case 3: {
        real__t mu0 = Pmech0;
        gen_solution[bus_name][GP1_idx] = 0.;
        gen_solution[bus_name][GP2_idx] = mu0;
        gen_solution[bus_name][hp_idx]  = Pmech0;
        gen_solution[bus_name][ip_idx]  = Pmech0;
        gen_solution[bus_name][lp_idx]  = Pmech0;
        break;
      }
      case 4:
      case 9: { // need to check with the control_option.
        real__t mu0 = Pmech0;
        gen_solution[bus_name][GP1_idx] = 0.;
        gen_solution[bus_name][GP2_idx] = mu0;
        gen_solution[bus_name][hp_idx]  = Pmech0;
        gen_solution[bus_name][ip_idx]  = Pmech0;
        gen_solution[bus_name][lp_idx]  = Pmech0;
        gen.mu0 = Pmech0;
        break;
      }
      case 5: { // need to check with the control_option.
        real__t mu0 = Pmech0;
        real__t control_option = all_gov_5[gen.GOV_Par].control_option;
        gen_solution[bus_name][GP1_idx] = 0.;
        gen_solution[bus_name][GP2_idx] = (control_option == 3) ? mu0 : 0.;
        gen_solution[bus_name][hp_idx]  = Pmech0;
        gen_solution[bus_name][ip_idx]  = Pmech0;
        gen_solution[bus_name][lp_idx]  = Pmech0;
        gen.mu0 = Pmech0;
        break;
      }
      case 7: {
        real__t mu0 = Pmech0;
        gen_solution[bus_name][GP1_idx] = 1.;
        gen_solution[bus_name][GP2_idx] = 0.;
        gen_solution[bus_name][GP3_idx] = mu0;
        gen_solution[bus_name][GPm_idx] = Pmech0;
        gen.mu0 = Pmech0;
        break;
      }
      case 8: {
        real__t mu0 = Pmech0;
        gen_solution[bus_name][GP1_idx] = 1.;
        gen_solution[bus_name][GP2_idx] = 0.;
        gen_solution[bus_name][GP3_idx] = mu0;
        gen_solution[bus_name][GPm_idx] = Pmech0;
        gen.mu0 = Pmech0;
        break;
      }
      default: {std::cerr << "Error: unsupported governor (GOV) type...\n"; std::terminate(); break;}
    }
    
    /** AVR */
    switch (gen.AVR_Model) {
      case 0: break;
      case 1:
        gen_solution[bus_name][EVr_idx] = 0;
        gen_solution[bus_name][EVa_idx] = 0;
        gen_solution[bus_name][EVe_idx] = Efd0;
        gen_solution[bus_name][EVf_idx] = 0;
        break;
      case 2:
        gen_solution[bus_name][EV1_idx] = 0;
        gen_solution[bus_name][EV2_idx] = 0;
        gen_solution[bus_name][EV3_idx] = 0;
        gen_solution[bus_name][EVa_idx] = 0;
        break;
      case 3: case 5: case 7: case 9: {
        real__t Vin0  = get_exc_input(Vd, Vq, Id, Iq, all_exc_3_10[gen.AVR_Par].Xc);
        real__t Se0   = all_exc_3_10[gen.AVR_Par].C1 * exp(all_exc_3_10[gen.AVR_Par].C2 * Efd0);
        real__t Vadj0 = Se0 * Efd0 + Ifd0 * all_exc_3_10[gen.AVR_Par].Kd;
        real__t Vfe0  = Vadj0 + Efd0 * all_exc_3_10[gen.AVR_Par].Ke;
        real__t Vr0   = Efd0 * all_exc_3_10[gen.AVR_Par].Ke + Vadj0;
        real__t Va0   = Vr0 / all_exc_3_10[gen.AVR_Par].KB + Vfe0 * all_exc_3_10[gen.AVR_Par].KH1;
        real__t V30   = Va0 / all_exc_3_10[gen.AVR_Par].Ka;
        gen_solution[bus_name][EV1_idx] = Vin0;
        gen_solution[bus_name][EV2_idx] = V30;
        gen_solution[bus_name][EV3_idx] = V30;
        gen_solution[bus_name][EVa_idx] = Va0;
        gen_solution[bus_name][EVr_idx] = Vr0;
        gen_solution[bus_name][EVe_idx] = Efd0;
        gen_solution[bus_name][EVf_idx] = 0;
        gen.Vt_ref = V30 * all_exc_3_10[gen.AVR_Par].Kv / all_exc_3_10[gen.AVR_Par].K + Vin0;
        break;
      }
      case 4: case 6: case 8: case 10: {
        real__t Vin0  = get_exc_input(Vd, Vq, Id, Iq, all_exc_3_10[gen.AVR_Par].Xc);
        real__t Se    = all_exc_3_10[gen.AVR_Par].C1 * exp(all_exc_3_10[gen.AVR_Par].C2 * Efd0);
        real__t Vadj0 = Se * Efd0 + Ifd0 * all_exc_3_10[gen.AVR_Par].Kd;
        real__t Vr0   = Efd0 * all_exc_3_10[gen.AVR_Par].Ke + Vadj0;
        real__t Va0   = Vr0 / all_exc_3_10[gen.AVR_Par].KB + Efd0 * all_exc_3_10[gen.AVR_Par].KH1;
        real__t V30   = Va0 / all_exc_3_10[gen.AVR_Par].Ka;
        gen_solution[bus_name][EV1_idx] = Vin0;
        gen_solution[bus_name][EV2_idx] = V30;
        gen_solution[bus_name][EV3_idx] = V30;
        gen_solution[bus_name][EVa_idx] = Va0;
        gen_solution[bus_name][EVr_idx] = Vr0;
        gen_solution[bus_name][EVe_idx] = Efd0;
        gen_solution[bus_name][EVf_idx] = 0;
        gen.Vt_ref = V30 * all_exc_3_10[gen.AVR_Par].Kv / all_exc_3_10[gen.AVR_Par].K + Vin0;
        break;
      }
      case 11: case 12: {
        real__t Vin0 = get_exc_input(Vd, Vq, Id, Iq, all_exc_11_12[gen.AVR_Par].Xc);
        real__t Vr0  = Efd0 + Ifd0 * all_exc_11_12[gen.AVR_Par].Kc;
        real__t V30  = Vr0 / all_exc_11_12[gen.AVR_Par].Ka;
        gen_solution[bus_name][EV1_idx] = Vin0;
        gen_solution[bus_name][EV2_idx] = V30;
        gen_solution[bus_name][EV3_idx] = V30;
        gen_solution[bus_name][EVa_idx] = Vr0;
        gen_solution[bus_name][EVf_idx] = 0;
        gen.Vt_ref = V30 * all_exc_11_12[gen.AVR_Par].Kv / all_exc_11_12[gen.AVR_Par].K + Vin0;
        break;
      }
      default: {std::cerr << "Error: unsupported excitor (AVR) type...\n"; std::terminate(); break;}
    }
    
    /** PSS */
    gen_solution[bus_name][PP1_idx] = 0;
    gen_solution[bus_name][PP2_idx] = 0;
    gen_solution[bus_name][PP3_idx] = 0;
    gen_solution[bus_name][PP4_idx] = 0;
    gen_solution[bus_name][PP5_idx] = 0;
    gen_solution[bus_name][PP6_idx] = 0;
    gen_solution[bus_name][PP7_idx] = 0;
    gen_solution[bus_name][PP8_idx] = 0;
    gen_solution[bus_name][VS_idx]  = 0;
    
    /** for output */
    gen_solution[bus_name][mu_output_idx]  = Pmech0;
    gen_solution[bus_name][PT_output_idx]  = Pmech0;
    gen_solution[bus_name][Efd_output_idx] = Efd0;
    gen_solution[bus_name][VS_output_idx]  = 0;
    
#if DEBUG
    cout << "bus_name = " << bus_name << endl;
    cout << "EQx = " << EQx << endl;
    cout << "EQy = " << EQy << endl;
    cout << "Edp0 = " << Edp0 << endl;
    cout << "Eqp0 = " << Eqp0 << endl;
    cout << "Edpp0 = " << Edpp0 << endl;
    cout << "Eqpp0 = " << Eqpp0 << endl;
    cout << "Efd0 = " << Efd0 << endl;
    cout << "mu0 = " << gen.mu0 << endl;
    
    cout << "Vm = " << Vm << endl;
    cout << "Va (degree) = " << Va * 180 / PI << endl;
    cout << "Pg = " << Pg << endl;
    cout << "Qg = " << Qg << endl;
    cout << "Vx = " << Vx << endl;
    cout << "Vy = " << Vy << endl;
    cout << "Ix = " << Ix << endl;
    cout << "Iy = " << Iy << endl;
    cout << "Im = " << sqrt(Ix * Ix + Iy * Iy) << endl;
    
    cout << "Vd = " << Vd << endl;
    cout << "Vq = " << Vq << endl;
    cout << "Id = " << Id << endl;
    cout << "Iq = " << Iq << endl << endl;
#endif
  }
  printf("\nInitial values of differential variables:\n");
  print_gen_solution();
}

}
