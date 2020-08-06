#include <transient.hpp>
//#include <trans_cuda.hpp>
#include <map>
#include <thrust/for_each.h>
#include <string>
#include <utility>
//#include <boost/functional/hash.hpp>
//#include <boost/numeric/odeint/integrate/integrate_const.hpp>
//#include <boost/numeric/odeint/external/thrust/thrust.hpp>
//#include <iostream>
//#include <cmath>
//#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
//#include <thrust/device_vector.h>
//#include <thrust/iterator/permutation_iterator.h>
//#include <thrust/iterator/counting_iterator.h>

//#include <boost/numeric/odeint/external/thrust/thrust.hpp>
//#include <thrust/for_each.h>
//#include <thrust/device_vector.h>
//#include <thrust/execution_policy.h>
using namespace std;
using namespace boost::numeric::odeint;
namespace transient_analysis {
/*
typedef pair<string, string> strpair;
struct printf_functor
{
  //template <string T1, GENERATOR T2>;
  __host__ __device__
  void operator()(const pair<string, string> mp)
  {
    // note that using printf in a __device__ function requires
    // code compiled for a GPU with compute capability 2.0 or
    // higher (nvcc --arch=sm_20)
    //printf("%d\n", x);
    //printf("%d\n", y);
    printf("mp first=%s\n", mp.first);
    printf("success in printf_functor\n");
    //uint__t bus_id = bus_name_to_id[bus_name] - 1;
  }
};
*/
//extern "C" void test_cuda(unordered_map<int, int> generators);

//extern "C" void test_cuda(map<string, GENERATOR> generators);
//void test_cuda(map<string, GENERATOR> &generators);

//extern __global__ void test_cuda();

Transient_Solver::
Transient_Solver(real__t start_time, real__t end_time, real__t time_stepping, string main_path)
    : start_time(start_time), end_time(end_time), time_stepping(time_stepping), main_path(main_path) {}

Transient_Solver::~Transient_Solver() {
  GraphLU_Destroy(Ybus_matrix);
  free(Ybus_matrix);
  free(eY);
  free(ei);
  free(ep);
  free(gCurrent);
  free(gVoltage);
}

void Transient_Solver::setup_system() {
  load_system_data(main_path);
  print_system_summary();
  get_bus_name_id_mapping();
  bus_types[1] = "PQ Bus";
  bus_types[2] = "PV Bus";
  bus_types[3] = "Slack Bus";
  
  string output_path = "../output";
  remove((output_path + "/genBusResults.csv").c_str());
  remove((output_path + "/allBusResults.csv").c_str());
  
  nBuses       = buses.size();
  nGen         = generators.size();
  nGenUnknowns = VS_output_idx + 1;
  n            = 2 * nBuses;
  nnz          = 4 * (line.size() + nBuses);
    
  for (auto &g_hldr : generators) {
    auto & bus_name=g_hldr.first;
    auto & gen=g_hldr.second;
    gen_solution  [bus_name] = vector<real__t>(nGenUnknowns, 0.);
    gen_error     [bus_name] = vector<real__t>(nGenUnknowns, 0.);
    gen_dq_current[bus_name] = vector<real__t>(2, 0.);
  }
  bus_voltage = vector<vector<real__t>> (nBuses, vector<real__t>(2, 0.));
  //  dump_Voltages();

  print_gen_parameters();
//  print_branch_data();
  
  gCurrent = (real__t*)calloc(n, sizeof(real__t));
  gVoltage = (real__t*)calloc(n, sizeof(real__t));

  /** allocate memory for arrays used in CSR form */
  eY = (real__t*)malloc(sizeof(real__t) * (nnz));    // nonzero values in Ybus matrix
  ei = (uint__t*)malloc(sizeof(uint__t) * (nnz));    // column idx of Ybus matrix
  ep = (uint__t*)malloc(sizeof(uint__t) * (n + 1));  // initial row pointers

  Ybus_matrix = (SGraphLU*)malloc(sizeof(SGraphLU));
  gVertex_all = vector<vector<real__t>>(nBuses, vector<real__t>(5, 0.));
  
  current_time  = start_time;
  num_of_steps  = 0;
  tol           = 1.e-4;
  max_step_size = 0.01;

  is_modify_Y_bus_matrix = true;
  is_fault_treated = false;
}

void Transient_Solver::load_system_data(string folder_path) {
  read_bus_data(buses, folder_path + "/system/Node.csv");
  read_generator_node_data(buses, generators, folder_path + "/system/Generator.csv");
  read_fdpf_data(buses, "../output/fdpf.csv");
//  read_load_data(buses, folder_path + "/system/Load.csv");
  read_compensator_P_data(buses, folder_path + "/system/Compensator_P.csv");
//  read_DC_line_data(buses, line, folder_path + "/system/DC_Line.csv");
  read_AC_line_data(buses, line, folder_path + "/system/AC_Line.csv");
  read_two_winding_transformer_data(buses, line, folder_path + "/system/Two_winding_transformer.csv");
  read_three_winding_transformer_data(buses, line, folder_path + "/system/Three_winding_transformer.csv");
  
  read_EPRI_GEN_data(all_gen, folder_path + "/parameter/Synchronous_Machine.csv");
  read_EPRI_GOV_I_data(all_gov_1, folder_path + "/parameter/Governor_1.csv");
  read_EPRI_GOV_II_data(all_gov_2, folder_path + "/parameter/Governor_2.csv");
  read_EPRI_GOV_III_data(all_gov_3, folder_path + "/parameter/Governor_3.csv");
  read_EPRI_GOV_IV_data(all_gov_4, folder_path + "/parameter/Governor_4.csv");
  read_EPRI_GOV_V_data(all_gov_5, folder_path + "/parameter/Governor_5.csv");
  read_EPRI_GOV_VII_data(all_gov_7, folder_path + "/parameter/Governor_7.csv");
  read_EPRI_GOV_VIII_data(all_gov_8, folder_path + "/parameter/Governor_8.csv");
  read_EPRI_GOV_IX_data(all_gov_9, folder_path + "/parameter/Governor_9.csv");
  read_EPRI_EXC_I_data(all_exc_1, folder_path + "/parameter/AVR_1.csv");
  read_EPRI_EXC_II_data(all_exc_2, folder_path + "/parameter/AVR_2.csv");
  read_EPRI_EXC_III_TO_X_data(all_exc_3_10, folder_path + "/parameter/AVR_3_to_10.csv");
  read_EPRI_EXC_XI_TO_XII_data(all_exc_11_12, folder_path + "/parameter/AVR_11_to_12.csv");
  read_EPRI_PSS_I_data(all_pss_1, folder_path + "/parameter/PSS_1.csv");
  read_EPRI_PSS_II_data(all_pss_2, folder_path + "/parameter/PSS_2.csv");
  read_EPRI_PSS_IV_VI_data(all_pss_4_6, folder_path + "/parameter/PSS_4_6.csv");
  read_EPRI_PSS_V_data(all_pss_5, folder_path + "/parameter/PSS_5.csv");
  read_EPRI_PSS_VIII_data(all_pss_8, folder_path + "/parameter/PSS_8.csv");

#if DEBUG
  cout << "Transient simulation loading data success!\n";
#endif
}

void Transient_Solver::get_bus_name_id_mapping() {
  int idx = 1; // id begins with 1
  for (auto &b : buses) {
    bus_name_to_id[b.first] = idx;
    bus_id_to_name[idx] = b.first;
    ++idx;
  }
#if DEBUG
  printf("Transient simulation getting bus_name to bus_id mapping success!\n");
#endif
}

void Transient_Solver::dump_Voltages() {
  if (current_time == start_time) {
    for (int i = 0; i < buses.size(); ++i) {
      string bus_name = bus_id_to_name[i];
      bus_voltage[i][0] = buses[bus_name].Vm;
      bus_voltage[i][1] = buses[bus_name].Va;
    }
  } else {
    for (int i = 0; i < buses.size(); ++i) {
      real__t Vx = gVoltage[i];
      real__t Vy = gVoltage[i + nBuses];
      bus_voltage[i][0] = sqrt(Vx * Vx + Vy * Vy);
      bus_voltage[i][1] = atan2(Vy, Vx);
    }
  }
#if DEBUG
  printf("Transient simulation dumping voltages success!\n");
#endif
}

string Transient_Solver::
solve_algebraic_equations_one_step(SGraphLU* matrix, real__t* rhs) {
  /* call the graphlu solver to solve the sparse linear system
   * notice that GraphLU_Solve_Singular allows the matrix to be *numerically*
   * singular */
  int ret_solve = GraphLU_Solve_Singular(matrix, rhs, 0);
  if (ret_solve < 0) {
    printf("Error: solve_algebraic_equations_one_step: %d\n", ret_solve);
    return "FAILED";
  }
  return "SUCCESS";
}

/* a matrix-vector multiplier (not used in this code) */
void Transient_Solver::
matrix_vector_mult(const SGraphLU* matrix, const real__t* src, real__t* dst) {
  uint__t  n  = matrix->n;
  real__t* ax = matrix->ax;
  uint__t* ai = matrix->ai;
  uint__t* ap = matrix->ap;

  for (int i = 0; i < n; ++i) {
    dst[i] = 0;
    for (int k = ap[i]; k < ap[i + 1]; ++k) {
      dst[i] += ax[k] * src[ai[k]];
    }
  }
}

void Transient_Solver::update_step_size() {
  real__t err = 0;
  for (auto &g : gen_error) {
    for (auto &v : g.second) err += abs(v * v);
  }
  err = max(sqrt(err), EPS);
  printf("current error = %2.12f\n", err);
  
  real__t q = pow(tol / err, 1. / 4.) * 0.8;
  if (err < tol)
    step_size = min(max(q, 0.1), 4.) * step_size;
  else
    step_size = min(max(q, 0.1), 1.) * step_size;

  step_size = min(max_step_size, step_size);
}

void Transient_Solver::get_generator_current() {
  for (auto &g_hldr : generators) {
    auto & bus_name=g_hldr.first;
    auto & gen=g_hldr.second;
    uint__t bus_id = bus_name_to_id[bus_name] - 1;
    uint__t Gen_Par = gen.Gen_Par;
    if (Gen_Par == 0) continue;
    
    real__t Vx    = gVoltage[bus_id];
    real__t Vy    = gVoltage[bus_id + nBuses];
    real__t delta = gen_solution[bus_name][delta_idx];
    
    real__t Vd = Vx * sin(delta) - Vy * cos(delta);
    real__t Vq = Vx * cos(delta) + Vy * sin(delta);
    real__t Ra = all_gen[Gen_Par].Ra;
    
    real__t Xdpp = all_gen[Gen_Par].Xdpp;
    real__t Xqpp = all_gen[Gen_Par].Xqpp;
    real__t Edpp = gen_solution[bus_name][Edpp_idx];
    real__t Eqpp = gen_solution[bus_name][Eqpp_idx];
    
    real__t denom = Ra * Ra + Xdpp * Xqpp;
    assert(denom > EPS);
    
    gen_dq_current[bus_name][0] = (+Ra * (Edpp - Vd) + Xqpp * (Eqpp - Vq)) / denom;
    gen_dq_current[bus_name][1] = (-Xdpp * (Edpp - Vd) + Ra * (Eqpp - Vq)) / denom;
    
//    cout << "dq currents in get_generator_current:" << endl;
//    cout << "Id, Iq = " << gen_dq_current[bus_name][0] << ", " << gen_dq_current[bus_name][1] << endl;

    #if DEBUG
      printf("Transient simulation getting generator current success for bus %s!\n", bus_name.c_str());
    #endif
  }
}

void Transient_Solver::
apply_fault(const uint__t busID, const real__t t, const real__t fault_btime, const real__t fault_etime) {
  /* busID is indexed from 1 */

  if (t >= fault_btime && t < fault_etime - time_stepping) {  // fault happens
    if (is_fault_treated)
      return;
    for (int i = 0; i < gEdge_all.size(); ++i) {
      if (gEdge_all[i][0] == busID - 1 && gEdge_all[i][1] == busID - 1) {
        gEdge_all[i][2] += INF;
        is_fault_treated = true;
        break;
      }
    }
  } else if (abs(t - fault_etime) <= EPS) {  // fault cleared
    is_fault_treated = false;
  } else if (t > fault_etime) {
    if (is_fault_treated)
      return;
    for (int i = 0; i < gEdge_all.size(); ++i) {
      if (gEdge_all[i][0] == busID - 1 && gEdge_all[i][1] == busID - 1) {
        gEdge_all[i][2] -= INF;
        is_fault_treated = true;
        break;
      }
    }
  }
}

/** this function should be called after the compute_ODE_initial_values */
void Transient_Solver::assign_ode_solver_data(string bus_name, GENERATOR& gen) {
  system.parameters.EXC_type = gen.AVR_Model;
  switch (gen.AVR_Model) {
    case 1:
      system.parameters.exc_1 = all_exc_1[gen.AVR_Par]; break;
    case 2:
      system.parameters.exc_2 = all_exc_2[gen.AVR_Par]; break;
    case 3: case 4: case 5: case 6: case 7: case 8: case 9: case 10:
      system.parameters.exc_3_10 = all_exc_3_10[gen.AVR_Par]; break;
    case 11: case 12:
      system.parameters.exc_11_12 = all_exc_11_12[gen.AVR_Par]; break;
    default:
      break;
  }
  
  system.parameters.GOV_type = gen.GOV_Model;
  switch (gen.GOV_Model) {
    case 1: system.parameters.gov_1 = all_gov_1[gen.GOV_Par]; break;
    case 2: system.parameters.gov_2 = all_gov_2[gen.GOV_Par]; break;
    case 3: system.parameters.gov_3 = all_gov_3[gen.GOV_Par]; break;
    case 4: system.parameters.gov_4 = all_gov_4[gen.GOV_Par]; break;
    case 5: system.parameters.gov_5 = all_gov_5[gen.GOV_Par]; break;
    case 7: system.parameters.gov_7 = all_gov_7[gen.GOV_Par]; break;
    case 8: system.parameters.gov_8 = all_gov_8[gen.GOV_Par]; break;
    case 9: system.parameters.gov_9 = all_gov_9[gen.GOV_Par]; break;
  }

  system.parameters.PSS_type = gen.PSS_Model;
  switch (gen.PSS_Model) {
    case 1: system.parameters.pss_1   = all_pss_1[gen.PSS_Par]; break;
    case 2: system.parameters.pss_2   = all_pss_2[gen.PSS_Par]; break;
    case 4: system.parameters.pss_4_6 = all_pss_4_6[gen.PSS_Par]; break;
    case 5: system.parameters.pss_5   = all_pss_5[gen.PSS_Par]; break;
    case 6: system.parameters.pss_4_6 = all_pss_4_6[gen.PSS_Par]; break;
    case 8: system.parameters.pss_8   = all_pss_8[gen.PSS_Par]; break;
  }
  
  system.parameters.GEN_type = gen.Gen_Model;
  if (gen.Gen_Par == 0) {
    all_gen[0].Xdp  = gen.Xdp < EPS  ? 0.0001     : gen.Xdp;
    all_gen[0].Xdpp = gen.Xdpp < EPS ? 0.0001     : gen.Xdpp;
    all_gen[0].TJ   = gen.TJ < EPS   ? 999999.875 : gen.TJ;
    all_gen[0].X2   = gen.X2;
    all_gen[0].Ra   = 0.;
  }
  system.parameters.gen = all_gen[gen.Gen_Par];
  system.parameters.gen.bus_id = bus_name_to_id[bus_name] - 1;
  
  system.parameters.omega_ref = gen.omega_ref;
  system.parameters.freq_ref  = gen.freq_ref;
  system.parameters.Pe_ref    = gen.Pe_ref;
  system.parameters.Vt_ref    = gen.Vt_ref;
  system.parameters.Efd0      = gen.Efd0;
  system.parameters.mu0       = gen.mu0;
  system.parameters.Rate_MW   = gen.Rate_MW;
}

void Transient_Solver::run(int argc, char** argv) {
  /** initialize simulation settings */
  setup_system();
  compute_ODE_initial_values();
  runge_kutta_dopri5<d_vector_type, value_type , d_vector_type , value_type> dopri5_stepper_type;
  print_bus_data();
  
  /** convert the data on the edges, notice that the conversion is
  once for all, unless the topology of the network changes */
  generate_edges_matrix();
  const value_type dt = 0.1;
  
  while (current_time <= end_time) {

    /* output every 10 steps - modify if necessary */
    if ((num_of_steps % 5 == 0 && current_time > 0) || num_of_steps == 1) {
      output_data_to_csv_one_step();
    }
    
    /* make sure that current_time + dt does not exceed end_time */
    time_stepping = min(time_stepping, end_time - current_time);
    if (time_stepping == 0) break;
    
    /** add a fault to a LOAD bus */
    int fault_bus_id = 28;
    assert(fault_bus_id <= nBuses);
    apply_fault(fault_bus_id, current_time, 2., 2.2);

    /** update the buses info, then generate the new Y bus matrix */
    convert_nodes();
    convert_to_CSR_arrays();
    generate_Y_Bus_matrix();

    /** As graphlu returns the solution in-place, we copy gCurrent to gVoltage */
    memcpy(gVoltage, gCurrent, sizeof(*gCurrent) * n);
  
    /** solve one step of algebraic equations */
    string result = solve_algebraic_equations_one_step(Ybus_matrix, gVoltage);
    if (result == "FAILED") {
      std::cerr << "Solving algebraic equations FAILED, please check!!!\n";
      std::terminate();
    }
    
    //get_generator_current();
    //unordered_map<string, string>  mp = {{"1","12"}, {"2", "23"}, {"3", "34"}};   
   
    //mp['1'] = '12';
    //mp['2'] = '23';
    //mp['3'] = '34';
    /** solve one step of ordinary differential equations */
    //thrust::for_each(thrust::device, generators.begin(), generators.end(), printf_functor());
    //for (auto& g_hldr : mp) {
    //  auto & bus_name=g_hldr.first;
    //  auto & gen=g_hldr.second;
    //  std::cout << "bus_name address:" << &bus_name << std::endl;
    //  std::cout << "bus_name:" << bus_name << std::endl;
    //}
    //test_cuda(mp);
    //test_cuda(generators);
    //printf("Begin thrust for_each...");
    //thrust::for_each(mp.begin(), mp.end(), printf_functor());
    //thrust::for_each(generators.begin(), generators.end(), printf_functor());
    //printf("End thruest for_each **********");
    for (auto& g_hldr : generators) {
      auto & bus_name=g_hldr.first;
      auto & gen=g_hldr.second;

      uint__t bus_id = bus_name_to_id[bus_name] - 1;

      /** update the Vx and Vy values at each gen node */
      system.Vx = gVoltage[bus_id];
      system.Vy = gVoltage[bus_id + nBuses];
      system.Id = gen_dq_current[bus_name][0];
      system.Iq = gen_dq_current[bus_name][1];

      //cout << "dq currents in solving ODE:" << endl;
      //cout << "Vx, Vy = " << system.Vx << ", " << system.Vy << endl;
      //cout << "Id, Iq = " << system.Id << ", " << system.Iq << endl;

      assign_ode_solver_data(bus_name, gen);

      d_vector_type d_gen_solution_set = gen_solution[bus_name];
      d_vector_type d_gen_error_set = gen_error[bus_name];
      //integrate_const( dopri5_stepper_type , system , d_gen_solution_set , 0.0 , 1.0 , dt );
      dopri5_stepper_type.do_step(system, d_gen_solution_set, current_time,
                                  time_stepping, d_gen_error_set);
      //dopri5_stepper_type.do_step(system, gen_solution[bus_name], current_time,
      //                            time_stepping, gen_error[bus_name]);
      //integrate_adaptive( make_controlled( 1.0e-6 , 1.0e-6 , dopri5_stepper_type() ) , system , x , t , t + 1.0 , 0.1 );
      //cout << "End of dp_step===========================" << endl;
      /** post-process: prepare for outputs */
      gen_solution[bus_name][mu_output_idx]  = system.get_mu();
      gen_solution[bus_name][PT_output_idx]  = system.get_Pmech();
      gen_solution[bus_name][Efd_output_idx] = system.get_Efd();
      gen_solution[bus_name][VS_output_idx]  = system.get_VS();
      /* update the time and number of steps */
      //current_time += time_stepping;
      //num_of_steps++;
    }
    //cout << "\n+++Now number of steps: " << num_of_steps << endl;

    
#if DEBUG
    if (num_of_steps % 1 == 0) {
//      print_matrix<real__t>(gEdge_all, "System Matrix");
//      print_array<real__t>(gCurrent, buses.size(), "RHS: ");
//      print_array<real__t>(gVoltage, buses.size(), "Solution: ");

      print_bus_data();
      print_gen_solution();
    }
#endif

    /* update the time and number of steps */
    current_time += time_stepping;
    num_of_steps++;
  }

  printf("\n\nTransient Simulation Result:\n");
  print_bus_data();
  print_gen_solution();

  /* Done with the simulation! Congratulations! */
  cout << "\n+++Total number of steps: " << num_of_steps << endl;
}
}  // namespace transient_analysis
