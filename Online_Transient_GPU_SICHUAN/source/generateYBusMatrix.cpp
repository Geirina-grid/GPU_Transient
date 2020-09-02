/******************************************************************************
 * Copyright (c): (2020) GEIRINA
 * All rights reserved.
 * Project: Power system transient simulation
 *
 * - Author:      Peng Wei, peng.wei@geirina.net
 * - Created on:  Jan. 29, 2020
 * - Last Update: Feb. 15, 2020
 *
 * - This library is part of the Transient_Solver class. The Y bus matrix are
 *   established here.
 *
*******************************************************************************/

#include "transient.hpp"

namespace transient_analysis {

void Transient_Solver::add_to_gEdge_all() {
  /** initialize all diagonal entries of the Y bus matrix */
  gEdge_all.clear();
  for (int i = 0; i < nBuses; ++i) {
    gEdge_all.push_back({1. * i, 1. * i, 0., 0.});
  }
  
  for (auto &g_hldr : line) {
    auto & id=g_hldr.first;
    auto & e=g_hldr.second;
    int s = bus_name_to_id[e.from] - 1;
    int t = bus_name_to_id[e.to] - 1;
    real__t tap = e.tap / e.tap_count;
    tap = (abs(tap) < EPS) ? 1. : tap;
    real__t tap_abs = abs(tap);
    
    /** the off-diagonal entry (s, t) */
    gEdge_all.push_back({1. * s, 1. * t, -tap_abs * e.G, -tap_abs * e.B});
    
    /** the diagonal entry (s, s) */
    real__t tap_sq = tap > 0 ? tap * tap : 1.;
    gEdge_all[s][2] += tap_sq * e.G;
    gEdge_all[s][3] += tap_sq * e.B + e.B_half;
    
    /** add Gm, Bm */
    gEdge_all[t][2] += e.Gm;
    gEdge_all[t][3] += e.Bm;
    
    /** store the number of neighbors, needed for the CSR pointer */
    gVertex_all[s][0]++;
  }
}

void Transient_Solver::process_shunts_and_loads() {
  for (int i = 0; i < nBuses; ++i) {
    string bus_name = bus_id_to_name[i + 1];
    BUS& bus = buses[bus_name];
    real__t Vm_square = bus.Vm * bus.Vm;
    assert(Vm_square > EPS);
    gEdge_all[i][2] += bus.Gshunt;  // G shunt
    gEdge_all[i][3] += bus.Bshunt;  // B shunt
    
    /** NOTE: The load can either be added to the diagonal entries of the Y
     * bus matrix, or be converted into current and then added to the right
     * hand side. In our case, we choose adding to the diagonal entries of
     * the Y bus matrix.
     */
    assert(bus.bus_type <= 3 && bus.bus_type >= 1);
    if (bus.bus_type == 1) {  // 1 - load bus
      /** ATTENTION: the load buses sometimes have work stored on Pgen and Qgen
       * these values need to be substracted from Pload and Qload
       */
      gEdge_all[i][2] += +(bus.Pload - bus.Pgen) / Vm_square;  // Gii
      gEdge_all[i][3] += -(bus.Qload - bus.Qgen) / Vm_square;  // Bii
    } else if (bus.bus_type <= 3 && bus.bus_type >= 2) {              // gen bus
      /** ATTENTION: the gen buses may also have nonzero Pload and Qload */
      gEdge_all[i][2] += +bus.Pload / Vm_square;  // Gii
      gEdge_all[i][3] += -bus.Qload / Vm_square;  // Bii
    } else {
      std::cerr << "Error: unknown bus type... Program terminated!\n";
      std::terminate();
    }
  }
}

void Transient_Solver::generate_edges_matrix() {
  /** reset values in gVertex_all */
   for (int i = 0; i < gVertex_all.size(); ++i) {
     gVertex_all[i][0] = 0;
   }
   
  add_to_gEdge_all();
  process_shunts_and_loads();
  
#if DEBUG
//  print_matrix(gEdge_all, "gEdge_all matrix");
  cout << "Transient simulation generating edge matrix success!" << endl;
#endif
}

void Transient_Solver::convert_nodes() {
  /**
   *     bus_type -- 3, slack bus
   *              -- 2, generator bus (PV bus)
   *              -- 1, load bus (PQ bus)
   */
  
  for (auto& g_hldr : generators) {
    auto & bus_name=g_hldr.first;
    auto & gen=g_hldr.second;
    uint__t bus_id   = bus_name_to_id[bus_name] - 1;
    real__t rate_mva = 1.; //gen.Rate_MVA;
    real__t delta    = gen_solution[bus_name][delta_idx];
    real__t Ra       = all_gen[gen.Gen_Par].Ra;
    real__t Xdpp, Xqpp, Edpp, Eqpp;

    switch (gen.Gen_Model) {
      case 6: case 3: case -6: case -3:
        Xdpp = all_gen[gen.Gen_Par].Xdpp;
        Xqpp = all_gen[gen.Gen_Par].Xqpp;
        Edpp = gen_solution[bus_name][Edpp_idx];
        Eqpp = gen_solution[bus_name][Eqpp_idx];
        break;
      case 5: case -5:
        Xdpp = all_gen[gen.Gen_Par].Xdp;
        Xqpp = all_gen[gen.Gen_Par].Xqp;
        Edpp = gen_solution[bus_name][Edp_idx];
        Eqpp = gen_solution[bus_name][Eqp_idx];
        break;
      case 4: case -4:
        Xdpp = all_gen[gen.Gen_Par].Xdpp;
        Xqpp = all_gen[gen.Gen_Par].Xqpp;
        Edpp = 0;
        Eqpp = gen.Epp0;
        break;
      case 2: case 1: case -2: case -1:
        Xdpp = all_gen[gen.Gen_Par].Xdp;
        Xqpp = all_gen[gen.Gen_Par].Xqp;
        Edpp = 0.;
        Eqpp = gen_solution[bus_name][Eqp_idx];
        break;
      case 0:
        Xdpp = gen.Gen_Par == 0 ? gen.Xdp : all_gen[gen.Gen_Par].Xdp;
        Xqpp = gen.Gen_Par == 0 ? gen.Xdp : all_gen[gen.Gen_Par].Xdp;
        Edpp = 0.;
        Eqpp = gen.Ep0;
        break;
      default:
        std::cerr << "Error: unsupported generator (GEN) type...\n";
        std::terminate();
        break;
    }
    
    real__t denom = Ra * Ra + Xdpp * Xqpp;
    assert(abs(denom) > EPS);
    
//    real__t Vx = gVoltage[bus_id];
//    real__t Vy = gVoltage[bus_id + nBuses];
//    real__t Vd = Vx * sin(delta) - Vy * cos(delta);
//    real__t Vq = Vx * cos(delta) + Vy * sin(delta);
//
//    gen_dq_current[bus_name][0] = (+Ra * (Edpp - Vd) + Xqpp * (Eqpp - Vq)) / denom;
//    gen_dq_current[bus_name][1] = (-Xdpp * (Edpp - Vd) + Ra * (Eqpp - Vq)) / denom;
//
//    cout << "dq currents in generateYBusMatrix:" << endl;
//    cout << "Id, Iq = " << gen_dq_current[bus_name][0] << ", " << gen_dq_current[bus_name][1] << endl;
    
    real__t t1 = (Xdpp - Xqpp) / 2.;
    real__t t2 = (Xdpp + Xqpp) / 2.;
    
    real__t Gx = rate_mva * (+Ra - t1 * sin(2 * delta)) / denom;
    real__t Bx = rate_mva * (+t2 + t1 * cos(2 * delta)) / denom;
    real__t Gy = rate_mva * (+Ra + t1 * sin(2 * delta)) / denom;
    real__t By = rate_mva * (-t2 + t1 * cos(2 * delta)) / denom;
    
    real__t Expp = +Edpp * sin(delta) + Eqpp * cos(delta);
    real__t Eypp = -Edpp * cos(delta) + Eqpp * sin(delta);
    
    /* the modified current on the gen bus */
    gCurrent[bus_id]          = (Gx * Expp + Bx * Eypp);
    gCurrent[bus_id + nBuses] = (By * Expp + Gy * Eypp);
    
    /* store the values on the vertices for future use */
    gVertex_all[bus_id][1] = Gx;
    gVertex_all[bus_id][2] = Bx;
    gVertex_all[bus_id][3] = By;
    gVertex_all[bus_id][4] = Gy;
    
#if DEBUG
    if (bus_name == "10010") {
      cout << "bus_id = " << bus_id << endl;
      cout << "delta = " << delta << endl;
      cout << "Ix = " << gCurrent[bus_id] << endl;
      cout << "Iy = " << gCurrent[bus_id + nBuses] << endl;
      cout << "Gx = " << Gx << endl;
      cout << "Bx = " << Bx << endl;
      cout << "Gy = " << Gy << endl;
      cout << "By = " << By << endl;
    }
#endif
  }
  
#if DEBUG
//  print_matrix(gVertex_all, "gVertex_all matrix");
  cout << "Transient simulation convert buses success!" << endl;
#endif
}

void Transient_Solver::convert_to_CSR_arrays() {
  /* edges need to be sorted */
  sortrows(gEdge_all);

  /* row turnning pointer (ep) */
  ep[0] = 0;
  uint__t idx = 1;
  /* the real part */
  for (int i = 0; i < nBuses; ++i, ++idx) {
    ep[idx] = ep[idx - 1] + 2 * (1 + (int)gVertex_all[i][0]);
  }
  /* the imaginary part */
  for (int i = 0; i < nBuses; ++i, ++idx) {
    ep[idx] = ep[idx - 1] + 2 * (1 + (int)gVertex_all[i][0]);
  }

  /* row index (ei) and value (eY) */
  idx = 0;
  /* the real part */
  for (uint__t i = 0; i < gEdge_all.size(); ++i, idx += 2) {
    int from = (int)gEdge_all[i][0];
    int to   = (int)gEdge_all[i][1];
    ei[idx] = to;
    ei[idx + 1] = to + nBuses;
    
    eY[idx] = gEdge_all[i][2];              // = Gik
    eY[idx + 1] = -gEdge_all[i][3];         // = -Bik
    if (from == to) {                       // diagonal entries
      eY[idx] += gVertex_all[from][1];      // += Gx
      eY[idx + 1] += gVertex_all[from][2];  // += Bx
    }
  }
  /* the imaginary part */
  for (uint__t i = 0; i < gEdge_all.size(); ++i, idx += 2) {
    int from = gEdge_all[i][0];
    int to   = gEdge_all[i][1];

    ei[idx] = to;
    ei[idx + 1] = to + nBuses;
    
    eY[idx] = gEdge_all[i][3];              // = Bik
    eY[idx + 1] = gEdge_all[i][2];          // = Gik
    if (from == to) {                       // diagonal entries
      eY[idx] += gVertex_all[from][3];      // += By
      eY[idx + 1] += gVertex_all[from][4];  // += Gy
    }
  }
#if DEBUG
  printf("Transient simulation generating CSR arrays success!\n");
#endif
}

void Transient_Solver::generate_Y_Bus_matrix() {
#if DEBUG
  print_line();
  cout << "    Y bus Matrix number of rows: " << n << endl
       << "    Total number of nonzeros: " << nnz << endl;
  print_line();
#endif

  assert(eY && ei && ep);

//#if DEBUG
//  std::cout << "n = " << n << "\tnnz = " << nnz << std::endl;
//  printf(" i\teY\tei\tep:\n");
//  for (int i = 0; i < nnz; ++i) {
//    if (i <= n)
//      printf("%2d\t%.4g\t%6d\t%6d\n", i, eY[i], ei[i], ep[i]);
//    else
//      printf("%2d\t%.4g\t%6d\n", i, eY[i], ei[i]);
//  }
//  std::cout << std::endl;
//#endif

  /** initialization */
  int ret_p = GraphLU_Initialize(Ybus_matrix);
  if (ret_p < 0) {
    printf("Initialize error: %d\n", ret_p);
  }

  /** sparse matrix created from eY, ei, ep */
  ret_p = GraphLU_CreateMatrix(Ybus_matrix, n, nnz, eY, ei, ep);
  if (ret_p < 0) {
    printf("CreateMatrix error: %d\n", ret_p);
  }

  /** analyze the sparse matrix */
  ret_p = GraphLU_Analyze(Ybus_matrix);
  if (ret_p < 0) {
    printf("GraphLU_Analyze error: %d\n", ret_p);
  }

  /** set control parameters */
  Ybus_matrix->cfgi[9] = 1;

  /** symbolic factorize */
  ret_p = Graphlu_Util::symbolic_factorize(Ybus_matrix);
  if (ret_p < 0) {
    printf("GraphLU Symbolic factorize error: %d\n", ret_p);
  }

  /** numerical factorize */
  ret_p = Graphlu_Util::numerical_factorize(Ybus_matrix);
  if (ret_p < 0) {
    printf("GraphLU Numerical factorize error: %d\n", ret_p);
  }
}

}
