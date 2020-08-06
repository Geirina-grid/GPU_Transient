/******************************************************************************
 * Copyright (c): (2019-2020) GEIRINA
 * All rights reserved.
 * Project: Power system transient simulation
 *
 * - Author:      Peng Wei, peng.wei@geirina.net
 * - Created on:  June 21, 2019
 * - Last Update: Jan. 29, 2020
 *
 * The equations are stacked in order corresponding to the structure:
 *
 *          +--------------------------------+
 *          |                                |
 *          V                                |
 *      Governor --> Turbine --> Generator ==+------> Load
 *                                  ^        |
 *                                  |        |
 *                              Excitation <-+
 *                                  ^        |
 *                                  |        |
 *                                 PSS <-----+
 *
 * Remarks 7/28/2019: 1. Xl was not considered in all cases
 *                    2. Shunt in 10Gen39Bus was not considered
 *
 *
 *******************************************************************************/

#ifndef TRANSIENT_HPP_
#define TRANSIENT_HPP_

#include "GraphDyn_util.hpp"
#include "ode_solver.hpp"
#include "readData.hpp"
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>
#include <iostream>
#include <cmath>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <thrust/device_vector.h>

//#include "printData.hpp"
//#include "outputData.hpp"
//#include "computeInitialValue.hpp"
//#include "generateYBusMatrix.hpp"

using namespace std;
using namespace boost::numeric::odeint;
#ifndef MAX_VERTICES
#define MAX_VERTICES 15000
#endif

namespace transient_analysis {

class Transient_Solver {
public:
  Transient_Solver(real__t start_time, real__t end_time, real__t time_stepping, string main_path);
  ~Transient_Solver();
  void run(int argc, char** argv);

 private:
  void load_system_data(string folder_path);
  
  void get_bus_name_id_mapping();
  void dump_Voltages();
  void setup_system();
  void add_to_gEdge_all();
  void process_shunts_and_loads();
  
  void generate_edges_matrix();
  void convert_nodes();
  void convert_to_CSR_arrays();
  void generate_Y_Bus_matrix();
  string solve_algebraic_equations_one_step(SGraphLU* matrix, real__t* rhs);
  void compute_ODE_initial_values();
  real__t get_exc_input(real__t Vd, real__t Vq, real__t Id, real__t Iq, real__t Xc);
  void output_data_to_csv_one_step();
  void update_step_size();
  void apply_fault(const uint__t busID, const real__t t, const real__t btime, const real__t etime);
  void matrix_vector_mult(const SGraphLU* matrix, const real__t* src, real__t* dst);
  void write_gen_info(const string path, const vector_type& x, const real__t t);
  void write_bus_info(const string path, real__t* gVoltage, const int n, const real__t t);
  void print_system_summary();
  void print_gen_solution();
  void print_bus_data();
  void print_gen_parameters();
  void print_branch_data();
  void get_generator_current();
  void assign_ode_solver_data(string bus_name, GENERATOR& gen);

  map<string, GENERATOR> generators;
  map<string, BUS> buses;
  unordered_map<string, int> bus_name_to_id;
  unordered_map<string, LINE> line;
  unordered_map<int, EPRI_GEN_DATA> all_gen;
  unordered_map<int, EPRI_EXC_I_DATA> all_exc_1;
  unordered_map<int, EPRI_EXC_II_DATA> all_exc_2;
  unordered_map<int, EPRI_EXC_III_TO_X_DATA> all_exc_3_10;
  unordered_map<int, EPRI_EXC_XI_TO_XII_DATA> all_exc_11_12;
  unordered_map<int, EPRI_GOV_I_DATA> all_gov_1;
  unordered_map<int, EPRI_GOV_II_DATA> all_gov_2;
  unordered_map<int, EPRI_GOV_III_DATA> all_gov_3;
  unordered_map<int, EPRI_GOV_IV_DATA> all_gov_4;
  unordered_map<int, EPRI_GOV_V_DATA> all_gov_5;
  unordered_map<int, EPRI_GOV_VII_DATA> all_gov_7;
  unordered_map<int, EPRI_GOV_VIII_DATA> all_gov_8;
  unordered_map<int, EPRI_GOV_IX_DATA> all_gov_9;
  unordered_map<int, EPRI_PSS_I_DATA> all_pss_1;
  unordered_map<int, EPRI_PSS_II_DATA> all_pss_2;
  unordered_map<int, EPRI_PSS_IV_VI_DATA> all_pss_4_6;
  unordered_map<int, EPRI_PSS_V_DATA> all_pss_5;
  unordered_map<int, EPRI_PSS_VIII_DATA> all_pss_8;
  
  map<string, vector_type> gen_solution;
  map<string, vector_type> gen_error;
  map<string, vector_type> gen_dq_current;

  vector<vector_type> bus_voltage;

  ODE_solver system;

  //unordered_map<string, int> bus_name_to_id;
  unordered_map<int, string> bus_id_to_name;
  unordered_map<int, string> bus_types;
  unordered_map<int, string>  bus_name_map;

  vector<vector<real__t>> gVertex_all;
  vector<vector<real__t>> gEdge_all;
  real__t* gCurrent;
  real__t* gVoltage;
  real__t* eY;
  uint__t* ei;
  uint__t* ep;
  SGraphLU* Ybus_matrix;
  uint__t n;
  uint__t nnz;

  uint__t nBuses;
  uint__t nGen;
  uint__t nGenUnknowns;

  real__t current_time;
  real__t start_time;
  real__t end_time;
  real__t time_stepping;
  uint__t num_of_steps;
  real__t tol;
  real__t step_size;
  real__t max_step_size;

  string main_path;

  bool is_modify_Y_bus_matrix;
  bool is_fault_treated;
};
}
#endif /* TRANSIENT_HPP_ */
