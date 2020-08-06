/******************************************************************************
 * Copyright (c): (2020) GEIRINA
 * All rights reserved.
 * Project: Power system transient simulation
 *
 * - Author:      Peng Wei, peng.wei@geirina.net
 * - Created on:  Jan. 29, 2020
 * - Last Update: Feb. 15, 2020
 *
 * - This library is part of the Transient_Solver class. The functions are
 *   responsible for printing calcuated result data to the console.
 *
*******************************************************************************/

#include <transient.hpp>

namespace transient_analysis {

  void Transient_Solver::print_system_summary() {
    printf("Transient simulation system summary:\n");
    TextTable t('-', '|', '+');
    t.addRow(vector<string>{"Component", "count"});
    t.addRow(vector<string>{"Generator", to_string(generators.size())});
    t.addRow(vector<string>{"Bus", to_string(buses.size())});
    t.addRow(vector<string>{"Gen Par", to_string(all_gen.size())});
    t.addRow(vector<string>{"Gov 1 Par", to_string(all_gov_1.size())});
    t.addRow(vector<string>{"Gov 2 Par", to_string(all_gov_2.size())});
    t.addRow(vector<string>{"Gov 3 Par", to_string(all_gov_3.size())});
    t.addRow(vector<string>{"Gov 4 Par", to_string(all_gov_4.size())});
    t.addRow(vector<string>{"Gov 5 Par", to_string(all_gov_5.size())});
    t.addRow(vector<string>{"Gov 7 Par", to_string(all_gov_7.size())});
    t.addRow(vector<string>{"Gov 8 Par", to_string(all_gov_8.size())});
    t.addRow(vector<string>{"Gov 9 Par", to_string(all_gov_9.size())});
    t.addRow(vector<string>{"Exc 1 Par", to_string(all_exc_1.size())});
    t.addRow(vector<string>{"Exc 2 Par", to_string(all_exc_2.size())});
    t.addRow(vector<string>{"Exc 3-10 Par", to_string(all_exc_3_10.size())});
    t.addRow(vector<string>{"Exc 11-12 Par", to_string(all_exc_11_12.size())});
    t.addRow(vector<string>{"PSS 1 Par", to_string(all_pss_1.size())});
    t.addRow(vector<string>{"PSS 2 Par", to_string(all_pss_2.size())});
    t.addRow(vector<string>{"PSS 4-6 Par", to_string(all_pss_4_6.size())});
    t.addRow(vector<string>{"PSS 5 Par", to_string(all_pss_5.size())});
    t.addRow(vector<string>{"PSS 8 Par", to_string(all_pss_8.size())});

    std::cout << t;
  }

  void Transient_Solver::print_bus_data() {
    TextTable t('-', '|', '+');
    t.addRow(vector<string>{"Bus ID", "Bus Name", "Bus Type", "Vm (p.u.)", "Va (rad)",
                            "Va (deg)", "Pgen", "Qgen", "Pload", "Qload"});
    vector<int> typeCount(4, 0);
    for (int i = 0; i < nBuses; ++i) {
      string bus_name = bus_id_to_name[i + 1];
      typeCount[buses[bus_name].bus_type]++;
      
      if (nBuses >= 5000 && i % 400 != 0) continue;
      
      real__t Vx = gVoltage[i], Vy = gVoltage[nBuses + i];
      real__t Vm = sqrt(Vx * Vx + Vy * Vy);
      real__t Va = atan2(Vy, Vx);
      
      t.add(to_string(i + 1));
      t.add(bus_name);
      t.add(bus_types[buses[bus_name].bus_type]);
      t.add(to_string(Vm));
      t.add(to_string(Va));
      t.add(to_string(Va * 180 / PI));
      t.add(to_string(buses[bus_name].Pgen));
      t.add(to_string(buses[bus_name].Qgen));
      t.add(to_string(buses[bus_name].Pload));
      t.add(to_string(buses[bus_name].Qload));
      t.endOfRow();
      
      cout
      << i << ", "
      << bus_name << ", "
      << bus_types[buses[bus_name].bus_type] << ", "
      << Vm << ", "
      << Va * 180 / PI << ", "
      << buses[bus_name].Pgen << ", "
      << buses[bus_name].Qgen << ", "
      << buses[bus_name].Pload << ", "
      << buses[bus_name].Qload
      << endl;
    }
    cout << "Total Slack Bus: " << typeCount[3] << endl;
    cout << "Total PV Bus: " << typeCount[2] << endl;
    cout << "Total PQ Bus: " << typeCount[1] << endl;
    cout << "Total Bus: " << typeCount[1] + typeCount[2] + typeCount[3] << endl;
    std::cout << t;
  }

//  void Transient_Solver::print_branch_data() {
//    TextTable t('-', '|', '+');
//
//    t.addRow(vector<string>{"From Name", " To Name", "From ID", "To ID", "G", "B", "B half", "tap", "tap count"});
//    t.addRow(vector<string>{" ", " ", " ", " ", "(p.u.)", "(p.u.)", "(p.u.)", " ", " "});
//
//    for (auto& p : lines) {
//      BRANCH line = p.second;
//      t.addRow(vector<string>{line.from, line.to, to_string(line.G), to_string(line.B),
//                              to_string(bus_namd_to_id[line.from]), to_string(bus_namd_to_id[line.to]),
//                              to_string(line.B_half), to_string(line.tap), to_string(line.tap_count)});
//    }
//    std::cout << t;
//  }

  void Transient_Solver::print_gen_solution() {
    TextTable t('-', '|', '+');

    t.addRow(vector<string>{"Bus ID", "Bus Name", "Bus Type", "delta", "omega",
                            "Edp", "Eqp", "Edpp", "Eqpp", "Efd",
                            "mu", "Pmech", "V_PSS", "Vt_ref"});
    t.addRow(vector<string>{" ", " ", " ", "(deg)", "(p.u.)",
                            "(p.u.)", "(p.u.)", "(p.u.)", "(p.u.)", "(p.u.)",
                            "(p.u.)", "(p.u.)", "(p.u.)", "(p.u.)"});
    int i = 0;
    for (auto &g_hldr : gen_solution) {
      auto & bus_name=g_hldr.first;
      auto & sol=g_hldr.second;
      if (gen_solution.size() >= 100 && ++i % 40 != 0) continue;

      uint__t bus_id = bus_name_to_id[bus_name];
      t.add("Bus " + to_string(bus_id));
      t.add(bus_name);
      t.add(bus_types[buses[bus_name].bus_type]);
      t.add(to_string(sol[delta_idx] * 180 / PI));
      t.add(to_string(sol[omega_idx]));
      t.add(to_string(sol[Edp_idx]));
      t.add(to_string(sol[Eqp_idx]));
      t.add(to_string(sol[Edpp_idx]));
      t.add(to_string(sol[Eqpp_idx]));
      t.add(to_string(sol[Efd_output_idx]));
      t.add(to_string(sol[mu_output_idx]));
      t.add(to_string(sol[PT_output_idx]));
      t.add(to_string(sol[VS_output_idx]));
      t.add(to_string(generators[bus_name].Vt_ref));
      t.endOfRow();
    }
    std::cout << t;
  }

  void Transient_Solver::print_gen_parameters() {
    printf("Generator parameters:\n");
    TextTable t('-', '|', '+');
    t.addRow(vector<string>{"Par_No", "XD", "XDP", "XDPP", "Xq",
                            "XqP", "XqPP", "X2", "Ra", "TD0P",
                            "TD0PP", "Tq0P", "Tq0PP", "TJ", "a",
                            "b", "n", "D_Sec"});
    t.addRow(vector<string>{" ", "(p.u.)", "(p.u.)", "(p.u.)", "(p.u.)",
                            "(p.u.)", "(p.u.)", "(p.u.)", "(p.u.)", "(s)",
                            "(s)", "(s)", "(s)", "(s)", "--",
                            "--", "--", "--"});
    for (auto &g_hldr : all_gen) {
      auto & par_no=g_hldr.first;
      auto & gen=g_hldr.second;
      t.add(to_string(par_no));
      t.add(to_string(gen.Xd));
      t.add(to_string(gen.Xdp));
      t.add(to_string(gen.Xdpp));
      t.add(to_string(gen.Xq));
      t.add(to_string(gen.Xqp));
      t.add(to_string(gen.Xqpp));
      t.add(to_string(gen.X2));
      t.add(to_string(gen.Ra));
      t.add(to_string(gen.Td0p));
      t.add(to_string(gen.Td0pp));
      t.add(to_string(gen.Tq0p));
      t.add(to_string(gen.Tq0pp));
      t.add(to_string(gen.TJ));
      t.add(to_string(gen.a));
      t.add(to_string(gen.b));
      t.add(to_string(gen.n));
      t.add(to_string(gen.D));
      t.endOfRow();
    }
    std::cout << t;
  }

}
