/******************************************************************************
 * Copyright (c): (2019-2020) GEIRINA
 * All rights reserved.
 * Project: Power system transient simulation
 *
 * - Author:      Peng Wei, peng.wei@geirina.net
 * - Created on:  Jan. 29, 2020
 * - Last Update: Jan. 29, 2020
 *
 * - This library is part of the Transient_Solver class. The functions are
 *   responsible for outputing calcuated result data into csv files.
 *
*******************************************************************************/

#include "transient.hpp"

namespace transient_analysis {

void Transient_Solver::
write_gen_info(const string path, const vector_type& x, const real__t t) {
  std::ofstream outfile;
  outfile.open(path, std::ios_base::app);
  outfile << t << ",\t";
  for (uint__t i = 0; i < x.size(); ++i)
    outfile << ((abs(x[i]) < EPS) ? 0. : x[i]) << ",\t";
  outfile << endl;
}

void Transient_Solver::
write_bus_info(const string path, real__t* gVoltage, const int n, const real__t t) {
  std::ofstream outfile;
  outfile.open(path, std::ios_base::app);

  outfile << t << ",\t";
  for (uint__t i = 0; i < n; ++i)
    outfile << ((abs(gVoltage[i]) < EPS) ? 0. : gVoltage[i]) << ",\t";
  outfile << endl;
}

void Transient_Solver::output_data_to_csv_one_step() {
  /** output: omega, mu, Pt, delta, Efd, VS, Vm(voltage at current bus) */
  int nOutput = 11, idx = 0;
  vector_type output_list(nOutput * nGen, 0);
  for (auto &g_hldr : generators) {
    auto & bus_name=g_hldr.first;
    auto & gen=g_hldr.second;
    uint__t bus_id = bus_name_to_id[bus_name] - 1;
    real__t Vx = gVoltage[bus_id];
    real__t Vy = gVoltage[bus_id + nBuses];
    
    output_list[idx + 0] = gen_solution[bus_name][omega_idx];
    output_list[idx + 1] = gen_solution[bus_name][mu_output_idx];
    output_list[idx + 2] = gen_solution[bus_name][PT_output_idx];
    output_list[idx + 3] = gen_solution[bus_name][delta_idx];
    output_list[idx + 4] = gen_solution[bus_name][Efd_output_idx];
    output_list[idx + 5] = gen_solution[bus_name][VS_output_idx];
    output_list[idx + 6] = sqrt(Vx * Vx + Vy * Vy);
    output_list[idx + 7] = gen_solution[bus_name][Edpp_idx];
    output_list[idx + 8] = gen_solution[bus_name][Eqpp_idx];
    output_list[idx + 9] = gen_solution[bus_name][Edp_idx];
    output_list[idx + 10] = gen_solution[bus_name][Eqp_idx];
    
    idx += nOutput;
  }
  write_gen_info("../output/genBusResults.csv", output_list, current_time);
  write_bus_info("../output/allBusResults.csv", gVoltage, nBuses * 2, current_time);
#if DEBUG
  cout << "Transient simulation output data to csv success!" << endl;
#endif
}

}
