/******************************************************************************
 * 
 * Project: Fast Decoupled Power Flow Computation
 *
 * - Author:      Peng Wei
 * - Created on:  Aug 23, 2019
 * - Last Update: Jan 14, 2019
 * 
 ** Update History:
 *  01/14/2020[Peng] - rewriting to comply with the new EPRI model data
 *        1. new csv parser is used
 ******************************************************************************
 */

#ifndef POWERFLOW_HPP_
#define POWERFLOW_HPP_

#include "GraphDyn_util.hpp"
#include "readData.hpp"

#ifndef MAX_VERTICES
#define MAX_VERTICES 15000
#endif

namespace transient_analysis {
class fast_decoupled_power_flow {
public:
  fast_decoupled_power_flow(string main_path);
  ~fast_decoupled_power_flow();
  void run();

private:
  void compute_deltaP_deltaQ(uint64_t n, real__t &maxDeltaP, real__t &maxDeltaQ,
                             real__t *eG, real__t *eB, uint__t *ei, uint__t *ep,
                             real__t *deltaP, real__t *deltaQ,
                             real__t *Vm, real__t *Va, real__t *Pn, real__t *Qn,
                             uint__t *btype) ;
  void compute_P_Q(uint64_t n,
                   real__t *eG, real__t *eB,
                   uint__t *ei, uint__t *ep,
                   real__t *Vm, real__t *Va,
                   real__t *Pn, real__t *Qn,
                   uint__t *btype);
  
  void apply_forward_backward_sub(uint64_t n, real__t *src,
                                  real__t *lx, uint__t *li, size_t *lp,
                                  real__t *ux, uint__t *ui, size_t *up,
                                  uint__t *rp, uint__t *cpi,
                                  real__t *rows, real__t *cols) ;
  int perform_LU_factorization(SGraphLU *matrix, uint64_t n, uint64_t nnz,
                               real__t *ax, uint__t *ai, uint__t *ap) ;
  void print_bus_data(int n, const real__t *Vm, const real__t *Va,
                      const real__t *Pn, const real__t *Qn,
                      const uint__t *btype) ;
  void output_powerflow_to_csv(string path, uint64_t n, real__t *Vm, real__t *Va, real__t *Pn, real__t *Qn);
  string fdpf_LU_factorize (int64_t nonZerosBp, int64_t nonZerosBpp,
                            vector<vector<real__t>> &gVertex_all,
                            vector<vector<real__t>> &gEdge_all,
                            real__t max_change_P, real__t max_change_Q,
                            int64_t maxiter, int64_t iteration,
                            string output_path) ;
  
  void load_system_data(string path);
  void get_bus_name_id_mapping();
  void setup_system();
  void process_shunts();
  void generate_edges_matrix();
  void generate_vertices_matrix();
  void add_to_gEdge_all();
  void remove_ZBR();
  
  map<string, GENERATOR> generators;
  map<string, BUS> buses;
  unordered_map<string, LINE> line;
  
  unordered_map<string, int> bus_name_to_id;
  unordered_map<int, string> bus_id_to_name;
  unordered_map<int, string> bus_types;
  
  vector<vector<real__t>> gVertex_all;
  vector<vector<real__t>> gEdge_all;
  vector<int> pqv_degree;
  vector<int> pq_degree;
  vector<int> all_degree;
  
  uint__t nBuses;
  uint__t nonZerosBp;
  uint__t nonZerosBpp;
  
  uint__t maxIter;
  uint__t iterations;
  
  real__t max_change_P;
  real__t max_change_Q;
  
  string main_path;
  string bus_data_path;
  string line_data_path;
  string output_path;
};

}
#endif
