#include "powerflow.hpp"

namespace transient_analysis {


void fast_decoupled_power_flow::
compute_P_Q(uint64_t n,
            real__t *eG, real__t *eB,
            uint__t *ei, uint__t *ep,
            real__t *Vm, real__t *Va,
            real__t *Pn, real__t *Qn,
            uint__t *btype) {

  for (int i = 0; i < n; ++i) {
    Pn[i] = Qn[i] = 0.;
    
      // Calculate network injections
    for (int j = ep[i]; j < ep[i + 1]; ++j) {
      int p = ei[j];
      real__t cos_val = cos(Va[i] - Va[p]);
      real__t sin_val = sin(Va[i] - Va[p]);
      
      Pn[i] += -Vm[i] * Vm[p] * (eG[j] * cos_val + eB[j] * sin_val);
      Qn[i] += -Vm[i] * Vm[p] * (eG[j] * sin_val - eB[j] * cos_val);
    }
  }
}


void fast_decoupled_power_flow::
compute_deltaP_deltaQ(uint64_t n,
                      real__t &maxDeltaP, real__t &maxDeltaQ,
                      real__t *eG, real__t *eB,
                      uint__t *ei, uint__t *ep,
                      real__t *deltaP, real__t *deltaQ,
                      real__t *Vm, real__t *Va,
                      real__t *Pn, real__t *Qn,
                      uint__t *btype) {
  // reset maxdeltaP and maxdeltaQ
  maxDeltaP = maxDeltaQ = 0;
  int maxPosP = -1, maxPosQ = -1;
  // Calculate deltaP and deltaQ
  for (int i = 0; i < n; ++i) {
    deltaP[i] = (btype[i] != 3) ? Pn[i] : 0.;
    deltaQ[i] = (btype[i] == 1) ? Qn[i] : 0.;
    
    // Calculate network injections
    if (btype[i] != 3) {   // non slack bus
      for (int j = ep[i]; j < ep[i + 1]; ++j) {
        int p = ei[j];
        real__t cos_val = cos(Va[i] - Va[p]);
        real__t sin_val = sin(Va[i] - Va[p]);
        
        deltaP[i] += -Vm[i] * Vm[p] * (eG[j] * cos_val + eB[j] * sin_val);
        if (btype[i] < 2)  // calculate Q for PQ buses
          deltaQ[i] += -Vm[i] * Vm[p] * (eG[j] * sin_val - eB[j] * cos_val);
      }
    }
    
    // Get maxDeltaP and maxDeltaQ
    if (maxDeltaP < abs(deltaP[i])) {
      maxPosP = i;
      maxDeltaP = abs(deltaP[i]);
    }
    
    if (maxDeltaQ < abs(deltaQ[i])) {
      maxPosQ = i;
      maxDeltaQ = abs(deltaQ[i]);
    }
    
//    maxDeltaP = max(maxDeltaP, abs(deltaP[i]));
//    maxDeltaQ = max(maxDeltaQ, abs(deltaQ[i]));
    
    
    assert(Vm[i] > EPS);
    deltaP[i] /= Vm[i];
    deltaQ[i] /= Vm[i];
  }
  cout << "maxDeltaP = " << maxDeltaP << ", maxDeltaQ = " << maxDeltaQ << endl;
  cout << "maxPosP = " << maxPosP << ", maxPosQ = " << maxPosQ << endl;
}

void fast_decoupled_power_flow::
apply_forward_backward_sub(uint64_t n, real__t *src,
                           real__t *lx, uint__t *li, size_t *lp,
                           real__t *ux, uint__t *ui, size_t *up,
                           uint__t *rp, uint__t *cpi,
                           real__t *rows, real__t *cols) {
  assert(n > 0);
  assert(src && lx && li && lp && ux && ui && up && rp && cpi && rows && cols);
  
  real__t *b = (real__t *)calloc(n, sizeof(real__t));
  assert(b);
  
  // Multiply src with row-scaling (rows) to get rhs
  for (int i = 0; i < n; ++i) {
    b[i] = src[rp[i]] * rows[rp[i]];
  }
  
  // Full forward substitution
  for (int i = 0; i < n; ++i) {
    real__t sum = 0., diag = 0.;
    for (int p = lp[i]; p < lp[i + 1]; ++p) {
      if (i != li[p]) {  // non-diagnal
        sum += lx[p] * b[li[p]];
      } else {
        diag = lx[p];
      }
    }
    b[i] -= sum;
    if (abs(diag) < EPS)
      cout << "warnning: " << diag << ", " <<  cpi[i]+1 << ", " <<bus_id_to_name[cpi[i]+1] << endl;
    assert(abs(diag) > EPS);
    b[i] /= diag;
  }
  
  // Full backward substitution
  for (int i = n - 1; i >= 0; --i) {
    real__t sum = 0.;
    for (int p = up[i]; p < up[i + 1]; ++p) {
      if (i != ui[p]) {  // non-diagnal
        sum += ux[p] * b[ui[p]];
      }
    }
    b[i] -= sum;
  }
  
  // Scale and permute back, then store the result to src
  for (int i = 0; i < n; ++i) {
    src[i] = b[cpi[i]] * cols[cpi[i]];
  }
}

int fast_decoupled_power_flow::
perform_LU_factorization(SGraphLU *matrix, uint64_t n, uint64_t nnz,
                         real__t *ax, uint__t *ai, uint__t *ap) {
#if DEBUG
  cout << "perform_LU_factorization!\n";
#endif

  int ret_p = GraphLU_Initialize(matrix);
  if (ret_p < 0) {
    printf("Initialize error: %d\n", ret_p);
  }
  
//#if DEBUG
//  std::cout << "n = " << n << "\tnnz = " << nnz << std::endl;
//  printf("i\tax\t\tai\tap:\n");
//  for (int i = 0; i < nnz; ++i) {
//    if (ai[i] != 2363) continue;
////    if (nnz > 10000 && i % 100 != 0) continue;
//    printf("%d\t%+.6f\t%d\t", i, ax[i], ai[i]);
//    if (i <= n) std::cout << ap[i];
//    if (i > 0 && ai[i] == ai[i-1]) cout << "   ===>>> duplicates warnning!!!";
//    std::cout << std::endl;
//  }
//  std::cout << std::endl;
//#endif
  
  // Create Matrix
  ret_p = GraphLU_CreateMatrix(matrix, n, nnz, ax, ai, ap);
  if (ret_p < 0) {
    printf("CreateMatrix error: %d\n", ret_p);
  }
  matrix->cfgf[0] = 0.001;  // partial pivoting tolerance, default = 0.001
  
  // Set control parameters
  matrix->cfgi[1] = 1;  // MC64 scaling indicator: 0 = no; 1 = yes (recommended)
  
  // Analyze (column/row ordering, scaling)
  ret_p = GraphLU_Analyze(matrix);
  if (ret_p < 0) {
    printf("GraphLU_Analyze error: %d\n", ret_p);
  }
#if DEBUG
  printf("Sparse Matrix Analysis Time: %.4g\n", matrix->stat[0]);
#endif
  
  /**
   * GraphLU_CreateScheduler creates the task scheduler for parallel LU factorization.
   * If we want to run parallel factorization or parallel re-factorization, this
   * function should be called after GraphLU Analyze.
  */
  int ret = GraphLU_CreateScheduler(matrix);
#if DEBUG
  printf("\nTime of creating scheduler: %.4g\n", matrix->stat[4]);
  printf("Suggestion: %s\n", ret == 0 ? "parallel" : "sequential");
#endif
  
  int error = 0;
  
  if (!ret) {  // parallel factorization
    
    /**
     * GraphLU_CreateThreads creates threads for parallel computation.
     * The second argument (thread) specifies the number of threads, including
     * the main thread. The last argument (check) specifies whether to check
     * the number of threads or not.
    */
    GraphLU_CreateThreads(matrix, 2, TRUE);
#if DEBUG
    printf("Total cores: %d, threads created: %d\n", (int)(matrix->stat[9]), (int)(matrix->cfgi[5]));
#endif
    
    /**
    * This function binds threads to cores (unbind = FALSE) or unbinds threads
    * from cores (unbind = TRUE).
    */
    GraphLU_BindThreads(matrix, FALSE);
    
    // Numerical LU factorization with partial pivoting, parallel
    error = GraphLU_Factorize_MT(matrix);
#if DEBUG
    printf("Parallel factorization time: %.4g\n", matrix->stat[1]);
#endif
    
    if (error < 0)
      cout << "GraphLU_Factorize_MT error:" << error << endl;
    
  } else {  // Sequential factorization
    error = GraphLU_Factorize(matrix);
#if DEBUG
    printf("Sequential factorization time: %.4g\n\n", matrix->stat[1]);
#endif
    
    if (error < 0)
      cout << "GraphLU_Factorize error:" << error << endl;
  }
  
#if DEBUG
  if (error >= 0) {
    print_line();
    printf("++++GraphLU Factorize Success!\n");
    cout << "    Number of total nonzeros after factorization: "
    << matrix->lu_nnz << "\n";
    cout << "    Number of nonzeros in L: " << matrix->l_nnz << "\n";
    cout << "    Number of nonzeros in U: " << matrix->u_nnz << "\n";
    print_line();
  }
#endif
  
  return error;
}

void fast_decoupled_power_flow::
print_bus_data(int n, const real__t *Vm, const real__t *Va,
               const real__t *Pn, const real__t *Qn,
               const uint__t *btype) {
  TextTable t('-', '|', '+');
  
  t.addRow(vector<string>{"Bus ID", "Bus Name", "Vm (p.u.)", "Va (deg)", "P (p.u.)",
                          "Q (p.u.)", "Bus Type"});
  vector<int> typeCount(4, 0);
  for (int i = 0; i < n; ++i) {
    typeCount[btype[i]]++;

//    if (Vm[i] >= 0) continue;
    if ((n >= 6000 && i % 10 != 0)) continue;
    t.add(to_string(i + 1));
    t.add(bus_id_to_name[i + 1]);
    t.add(to_string(Vm[i]));
    t.add(to_string(Va[i] * 180 / PI));
    t.add(to_string(Pn[i]));
    t.add(to_string(Qn[i]));
    t.add(bus_types[btype[i]]);
    t.endOfRow();
  }
  cout << "Total Slack Bus: " << typeCount[3] << endl;
  cout << "Total PV Bus: " << typeCount[2] << endl;
  cout << "Total PQ Bus: " << typeCount[1] << endl;
  cout << "Total Bus: " << typeCount[1] + typeCount[2] + typeCount[3] << endl;

  std::cout << t;
}

void fast_decoupled_power_flow::
output_powerflow_to_csv(string path, uint64_t n, real__t *Vm, real__t *Va, real__t *Pn, real__t *Qn) {
  remove((path).c_str());

  std::ofstream outfile;
  outfile.open(path, std::ios_base::app);
  
  outfile << "Bus_Name,Vm,Va,P,Q\n";
  for (int i = 0; i < n; ++i)
    outfile << bus_id_to_name[i+1] << ","
    << Vm[i] << "," << Va[i] / PI * 180 << ","
    << Pn[i] << "," << Qn[i]
    << std::endl;
  
  std::cout << "Data written to: " << path << std::endl;
}

string fast_decoupled_power_flow::
fdpf_LU_factorize (int64_t nonZerosBp, int64_t nonZerosBpp,
                   vector<vector<real__t>> &gVertex_all,
                   vector<vector<real__t>> &gEdge_all,
                   real__t max_change_P, real__t max_change_Q,
                   int64_t maxiter, int64_t iteration,
                   string output_path) {

  // -----------------------------------------------------------------
  //        Initialize variables and arrays
  // -----------------------------------------------------------------
  real__t maxDeltaP = 0;
  real__t maxDeltaQ = 0;
  string result = "FAILED";
  
  // Get the dimension and the nnz of the matrix B' and B"
  uint__t n      = gVertex_all.size();
  uint__t nnz    = nonZerosBp;
  uint__t n_pp   = gVertex_all.size();
  uint__t nnz_pp = nonZerosBpp;
  uint__t n_e    = gVertex_all.size();
  uint__t nnz_e  = gEdge_all.size();
  
  #if DEBUG
    printf("\n");
    print_line();
    printf("++++Start Running GRAPHLU_fdpf_LU_factorize function!\n");
    cout << "Bp   Number of rows: " << n
         << ",\tNumber of nonzeros: " << nonZerosBp << endl;
    cout << "Bpp  Number of rows: " << n_pp
         << ",\tNumber of nonzeros: " << nonZerosBpp << endl;
    cout << "YBus Number of rows: " << n_e
         << ",\tNumber of nonzeros: " << gEdge_all.size() << endl;
    print_line();
  #endif
  
  real__t *ax      = (real__t *)calloc(nnz, sizeof(real__t));       // values in B'
  uint__t *ai      = (uint__t *)calloc(nnz, sizeof(uint__t));       // column indices of B'
  uint__t *ap      = (uint__t *)calloc(n + 1, sizeof(uint__t));     // initial row pointers
  real__t *ax_pp   = (real__t *)calloc(nnz_pp, sizeof(real__t));    // values in B"
  uint__t *ai_pp   = (uint__t *)calloc(nnz_pp, sizeof(uint__t));    // column indices of B"
  uint__t *ap_pp   = (uint__t *)calloc(n_pp + 1, sizeof(uint__t));  // initial row pointers
  real__t *eG      = (real__t *)calloc(nnz_e, sizeof(real__t));     // G values in Y
  real__t *eB      = (real__t *)calloc(nnz_e, sizeof(real__t));     // B values in Y
  uint__t *ei      = (uint__t *)calloc(nnz_e, sizeof(uint__t));     // column indices of Y
  uint__t *ep      = (uint__t *)calloc(n_e + 1, sizeof(uint__t));   // initial row pointers
  real__t *deltaP  = (real__t *)calloc(n, sizeof(real__t));         // b in the Ax=b
  real__t *deltaQ  = (real__t *)calloc(n, sizeof(real__t));         // b in the Ax=b
  real__t *Vm      = (real__t *)calloc(n, sizeof(real__t));         // V magnitude
  real__t *Va      = (real__t *)calloc(n, sizeof(real__t));         // V angle
  real__t *Pn      = (real__t *)calloc(n, sizeof(real__t));         // bus active power
  real__t *Qn      = (real__t *)calloc(n, sizeof(real__t));         // bus reactive power
  uint__t *btype   = (uint__t *)calloc(n, sizeof(uint__t));         // bus type
  
  assert(ax && ai && ap && ax_pp && ai_pp && ap_pp && eG && eB && ei &&
         ep && deltaP && deltaQ && Vm && Va && Pn && Qn && btype);
  
  // --------------------------------------------------------------------------
  //                Sort all input HeapAccum
  // --------------------------------------------------------------------------
  // Input gEdge_all is assumed to be unsorted.
  // Sort the rows before assigning values to local arrays.
  sortrows(gEdge_all);
  
  // --------------------------------------------------------------------------
  //                Convert Ap, Ap_pp, ep and the Vertex
  // --------------------------------------------------------------------------
  ap[0] = ap_pp[0] = ep[0] = 0;
  for (int i = 0, k = 1; i < n; i++, k++) {
    ap   [k] = (int)gVertex_all[i][0] + ap[k - 1];
    ap_pp[k] = (int)gVertex_all[i][1] + ap_pp[k - 1];
    ep   [k] = (int)gVertex_all[i][2] + ep[k - 1];
    Vm   [i] = gVertex_all[i][3];
    Va   [i] = gVertex_all[i][4];
    Pn   [i] = gVertex_all[i][5];
    Qn   [i] = gVertex_all[i][6];
    btype[i] = gVertex_all[i][7];
  }
  
  printf("\nPower Flow Initial Values:\n");
  print_bus_data(n, Vm, Va, Pn, Qn, btype);
  
  // --------------------------------------------------------------------------
  //      Convert Ybus (ei, eG and eB), B' (ei and Bp_x) and B" (ei and Bpp_x)
  // --------------------------------------------------------------------------
  int i_p = 0, i_pp = 0;
  for (int i = 0; i < nnz_e; ++i) {
    int idx = (int)gEdge_all[i][4];
    ei[i] = idx;
    eG[i] = gEdge_all[i][2];
    eB[i] = gEdge_all[i][3];
    if (abs(gEdge_all[i][5]) > EPS) {
      ai[i_p] = idx;
      ax[i_p] = gEdge_all[i][5];
      ++i_p;
    }
    
    if (abs(gEdge_all[i][6]) > EPS) {
      ai_pp[i_pp] = idx;
      ax_pp[i_pp] = gEdge_all[i][6];
      ++i_pp;
    }
  }
  cout << "i_p = " << i_p << ", i_pp = " << i_pp << endl;
  assert(i_p == nonZerosBp && i_pp == nonZerosBpp);
  // --------------------------------------------------------------------------
  //                Call GRAPHLU and Factorize B' Matrix
  // --------------------------------------------------------------------------

#if DEBUG
  cout << "\n++++Start factorizing B' Matrix...\n";
#endif
  
  SGraphLU *graphlu_Bp = (SGraphLU *)malloc(sizeof(SGraphLU));
  assert(graphlu_Bp);
  
  int error_Bp = perform_LU_factorization(graphlu_Bp, n, nnz, ax, ai, ap);
  
  real__t *lx_Bp = nullptr, *ux_Bp = nullptr;
  uint__t *li_Bp = nullptr, *ui_Bp = nullptr;
  size_t  *lp_Bp = nullptr, *up_Bp = nullptr;
  
  GraphLU_DumpLU(graphlu_Bp, &lx_Bp, &li_Bp, &lp_Bp, &ux_Bp, &ui_Bp, &up_Bp);
  
  // get the permutation arrays, rp_Bp and cp_Bp may be different
  // if the matrix is not symmetric
  uint__t *rp_Bp   = graphlu_Bp->row_perm;
  uint__t *cpi_Bp  = graphlu_Bp->col_perm_inv;
  real__t *rows_Bp = graphlu_Bp->row_scale;
  real__t *cols_Bp = graphlu_Bp->col_scale_perm;
    
  // --------------------------------------------------------------------------
  //                Call GRAPHLU and Factorize B" Matrix
  // --------------------------------------------------------------------------
#if DEBUG
  cout << "\n++++Start factorizing B\" Matrix...\n";
#endif
  SGraphLU *graphlu_Bpp = (SGraphLU *)malloc(sizeof(SGraphLU));
  assert(graphlu_Bpp);
  
  int error_Bpp = perform_LU_factorization(graphlu_Bpp, n_pp, nnz_pp, ax_pp, ai_pp, ap_pp);
  
  real__t *lx_Bpp = nullptr, *ux_Bpp = nullptr;
  uint__t *li_Bpp = nullptr, *ui_Bpp = nullptr;
  size_t  *lp_Bpp = nullptr, *up_Bpp = nullptr;
    
  GraphLU_DumpLU(graphlu_Bpp, &lx_Bpp, &li_Bpp, &lp_Bpp, &ux_Bpp, &ui_Bpp, &up_Bpp);
  
  // Get the permutation arrays, rp and cp may be different if
  // the matrix is not symmetric.
  uint__t *rp_Bpp   = graphlu_Bpp->row_perm;
  uint__t *cpi_Bpp  = graphlu_Bpp->col_perm_inv;
  real__t *rows_Bpp = graphlu_Bpp->row_scale;
  real__t *cols_Bpp = graphlu_Bpp->col_scale_perm;
  
  // --------------------------------------------------------------------------
  //                Perform Newton's method
  // --------------------------------------------------------------------------

  if (error_Bp >= 0 && error_Bpp >= 0) { // if no error, start the iteration
#if DEBUG
    print_line();
    cout << "++++Start iteratively updating deltaP and deltaQ...\n";
    print_line();
#endif
    
    int iter = 0;
    for (; iter < maxiter; ++iter) {
      
#if DEBUG
      printf("\n\ncurrent status at step %d:\n", iter);
      print_bus_data(n, Vm, Va, Pn, Qn, btype);
#endif
      
      compute_deltaP_deltaQ(n, maxDeltaP, maxDeltaQ, eG, eB, ei, ep,
                            deltaP, deltaQ, Vm, Va, Pn, Qn, btype);
      
      // Check if the program converges
      if (maxDeltaP < max_change_P && maxDeltaQ < max_change_Q) {
        result = "OK";
        break;
      }
      
      // update Va, which is associated to B' matrix
      apply_forward_backward_sub(n, deltaP,
                                 lx_Bp, li_Bp, lp_Bp,
                                 ux_Bp, ui_Bp, up_Bp,
                                 rp_Bp, cpi_Bp, rows_Bp, cols_Bp);
      
      for (int i = 0; i < n; ++i) {
        Va[i] -= deltaP[i];
      }
      
      compute_deltaP_deltaQ(n, maxDeltaP, maxDeltaQ, eG, eB, ei, ep,
                            deltaP, deltaQ, Vm, Va, Pn, Qn, btype);
      
      // Check if the program converges
      if (maxDeltaP < max_change_P && maxDeltaQ < max_change_Q) {
        result = "OK";
        break;
      }
      
      // update Vm, which is associated to B" matrix
      apply_forward_backward_sub(n, deltaQ,
                                 lx_Bpp, li_Bpp, lp_Bpp,
                                 ux_Bpp, ui_Bpp, up_Bpp,
                                 rp_Bpp, cpi_Bpp, rows_Bpp, cols_Bpp);
      
      for (int i = 0; i < n; ++i) {
        Vm[i] -= deltaQ[i];
      }
    }
    
#if DEBUG
    printf( "\n+++ Terminated at iteration %d, with maxDeltaP = %e, maxDeltaQ = %e\n", iter, maxDeltaP, maxDeltaQ);
#endif
    
    // ------------------------------------------------------------------------
    //                Store the Result back to the array
    // ------------------------------------------------------------------------
    compute_P_Q(n, eG, eB, ei, ep, Vm, Va, Pn, Qn, btype);
    
    printf("\n\nPower Flow Results:\n");
    print_bus_data(n, Vm, Va, Pn, Qn, btype);
    
    
    for (int i = 0; i < n; ++i) {
      gVertex_all[i][3] = Vm[i];
      gVertex_all[i][4] = Va[i];  // in radian
    }
    
    iteration = iter;
    output_powerflow_to_csv(output_path + "/fdpf.csv", n, Vm, Va, Pn, Qn);

  } else { //factorization failed
    result = "FAILED";
    print_line();
    cout << "Factorization FAILED...\n";
    print_line();
  }
  
// EXIT_fdpf_factorize:
  GraphLU_Destroy(graphlu_Bp);
  GraphLU_Destroy(graphlu_Bpp);
  
  free(ax);    free(ai);    free(ap);
  free(ax_pp); free(ai_pp); free(ap_pp);
  free(graphlu_Bp); free(graphlu_Bpp);
  
  free(eG); free(eB); free(ei); free(ep);
  
  free(deltaP); free(deltaQ);
  free(Vm); free(Va); free(Pn); free(Qn);
  
  free(btype);
  
  return result;
}

fast_decoupled_power_flow::
fast_decoupled_power_flow(string main_path): main_path(main_path) {}

fast_decoupled_power_flow::~fast_decoupled_power_flow() {}

void fast_decoupled_power_flow::setup_system() {
  output_path = "../output";

  load_system_data(main_path);
  get_bus_name_id_mapping();
  bus_types[1] = "PQ Bus";
  bus_types[2] = "PV Bus";
  bus_types[3] = "Slack Bus";
  
  nBuses = buses.size();
  gVertex_all = vector<vector<real__t>>(nBuses, vector<real__t>(8, 0.));
  all_degree.resize(nBuses, 0);
  pqv_degree.resize(nBuses, 0);
  pq_degree.resize(nBuses, 0);

  nonZerosBp  = 0;
  nonZerosBpp = 0;
  maxIter     = 1000;
  iterations  = 0;
  max_change_P = max_change_Q = 1.e-4;
  
#if DEBUG
  printf("Fast decoupled power flow setting up system success!\n");
#endif
}

void fast_decoupled_power_flow::load_system_data(string folder_path) {
  read_bus_data(buses, folder_path + "/system/Node.csv");
  read_load_data(buses, folder_path + "/system/Load.csv");
  read_compensator_P_data(buses, folder_path + "/system/Compensator_P.csv");
  read_generator_node_data(buses, generators, folder_path + "/system/Generator.csv");
//  read_DC_line_data(buses, line, folder_path + "/system/DC_Line.csv");
  read_AC_line_data(buses, line, folder_path + "/system/AC_Line.csv");
  read_two_winding_transformer_data(buses, line, folder_path + "/system/Two_winding_transformer.csv");
  read_three_winding_transformer_data(buses, line, folder_path + "/system/Three_winding_transformer.csv");
//  remove_ZBR();
#if DEBUG
  printf("Fast decoupled power flow loading data success!\n");
#endif
}

void fast_decoupled_power_flow::get_bus_name_id_mapping() {
  int idx = 1; // id begins with 1
  for (auto &bus : buses) {
    bus_name_to_id[bus.first] = idx;
    bus_id_to_name[idx] = bus.first;
    ++idx;
  }
#if DEBUG
  printf("Fast decoupled power flow getting bus_name to bus_id mapping success!\n");
#endif
}

/**
 * gEdge_all contains:
 *  0     1   2   3   4    5    6
 * from, to, eG, eB, ei, B'x, B"x
 */
void fast_decoupled_power_flow::remove_ZBR() {
  cout << "line count original: " << line.size() << endl;
  int numZBR = 0;
  while (true) {
    ++numZBR;
    unordered_map<string, LINE> new_line = line;
    unordered_set<string> to_remove;
    int count = 0;
    for (auto g_hldr : new_line) {
    auto & key=g_hldr.first;
    auto & e=g_hldr.second;
      if (count >= 1) {
        for (auto &v : to_remove)
          buses.erase(v);
        
        break;
      }
//      if (e.from == e.to) line.erase(key);
      if (e.to == "d-gch51g_525.0000_445761" || e.from == "d-gch51g_525.0000_445761")
        cout << "d-gch51g_525.0000_445761 processed!\n";
      if (abs(e.B) >= 1.e6 - 1.e-2) {
        ++count;
//        cout << "removing line: " << key << endl;
        
        buses[e.from].Pgen    += buses[e.to].Pgen;
        buses[e.from].Qgen    += buses[e.to].Qgen;
        buses[e.from].Pload   += buses[e.to].Pload;
        buses[e.from].Qload   += buses[e.to].Qload;
        buses[e.from].Gshunt  += buses[e.to].Gshunt;
        buses[e.from].Bshunt  += buses[e.to].Bshunt;
        buses[e.from].invX    += buses[e.to].invX;
        buses[e.from].bus_type = buses[e.to].bus_type;
        to_remove.insert(e.to);
        cout << "bus removed: " << e.to << endl;
        line.erase(key);
        cout << "line count now: " << line.size()
             << " bus count now: " << buses.size() << endl;
        
        for (auto g_hldr : new_line) {
          auto & id=g_hldr.first;
          auto & v=g_hldr.second;
          if (key == id) continue;                // already taken care of
          if (v.from == e.to && v.to == e.from) { // the reverse edge
            line.erase(id);
          } else if (v.from == e.to) {
//            cout << "modifying line: " << id << endl;
            string id2 = e.from + "#" + v.to;
            line[id2].from    = e.from;
            line[id2].to      = v.to;
            line[id2].G      += v.G;
            line[id2].B      += v.B;
            line[id2].invX   += v.invX;
            line[id2].B_half += v.B_half;
            line[id2].Gm     += v.Gm;
            line[id2].Bm     += v.Bm;
            line[id2].tap    += v.tap;
            line[id2].tap_count++;

            line.erase(id);
          } else if (v.to == e.to) {
//            cout << "modifying line: " << id << endl;
            string id2 = v.from + "#" + e.from;
            line[id2].from    = v.from;
            line[id2].to      = e.from;
            line[id2].G      += v.G;
            line[id2].B      += v.B;
            line[id2].invX   += v.invX;
            line[id2].B_half += v.B_half;
            line[id2].Gm     += v.Gm;
            line[id2].Bm     += v.Bm;
            line[id2].tap    += v.tap;
            line[id2].tap_count++;

            line.erase(id);
          }
        }
      }
    }
    
    if (count == 0 || numZBR == 32) break;
  }
  
}

void fast_decoupled_power_flow::add_to_gEdge_all() {
  cout << "Total lines: " << line.size() << endl;
  cout << "Total buses: " << buses.size() << endl;
 
  gEdge_all.clear();
  /** diagonal entries */
  for (int i = 0; i < nBuses; ++i) {
    gEdge_all.push_back({1. + i, 1. + i, 0., 0., 1. * i, 0., 0.});
  }
  
  for (auto &g_hldr : line) {
    auto & id=g_hldr.first;
    auto & e=g_hldr.second;
    if (buses.count(e.from) == 0 || buses.count(e.to) == 0) cout << id << endl;
    assert(buses.count(e.from) );
    assert(buses.count(e.to));

    uint__t s     = bus_name_to_id[e.from] - 1;
    uint__t t     = bus_name_to_id[e.to] - 1;
    uint__t sType = buses[e.from].bus_type;
    uint__t tType = buses[e.to].bus_type;
    real__t tap   = e.tap / e.tap_count;
    // if (e.tap_count >= 4) cout << ">>>>>>" <<  e.from << ", " << e.to  << ", " << e.tap_count << endl;
    tap = (abs(tap) < EPS) ? 1. : tap;
    real__t tap_abs = abs(tap);

    /** the off-diagonal entry (s, t) */
    all_degree[s]++;
    if (sType != 3 && tType != 3) {   // both are non slack buses
      nonZerosBp++;  pqv_degree[s]++;
      if (sType <= 1 && tType <= 1) { // both are PQ buses
        nonZerosBpp++;  pq_degree[s]++;
        gEdge_all.push_back({1. + s, 1. + t, -e.G * tap_abs, -e.B * tap_abs, 1. * t, e.invX, -e.B * tap_abs});
      } else {
        gEdge_all.push_back({1. + s, 1. + t, -e.G * tap_abs, -e.B * tap_abs, 1. * t, e.invX, 0.});
      }
    } else {                           // at least one slack bus
      gEdge_all.push_back({1. + s, 1. + t, -e.G * tap_abs, -e.B * tap_abs, 1. * t, 0., 0.});
    }
    
    /** the diagonal entry (s, s) */
    real__t tap_sq = tap > 0 ? tap * tap : 1.;
    gEdge_all[s][2] += tap_sq * e.G;
    gEdge_all[s][3] += tap_sq * e.B + e.B_half;
    gEdge_all[s][5] -= e.invX;

    /** add Gm and Bm to diagonal entries */
    gEdge_all[t][2] += e.Gm;
    gEdge_all[t][3] += e.Bm;
  }
}

void fast_decoupled_power_flow::process_shunts() {
  /** take care of shunts on the buses -- they are added to the diagonal entries */
  for (int i = 0; i < nBuses; ++i) {
    string bus_name = bus_id_to_name[i + 1];
    BUS& bus = buses[bus_name];
    real__t Vm_sq = bus.Vm * bus.Vm;
    assert(Vm_sq > EPS);
    gEdge_all[i][2] += bus.Gshunt;  // G shunt
    gEdge_all[i][3] += bus.Bshunt;  // B shunt
    gEdge_all[i][5] -= bus.invX;
  }
}

void fast_decoupled_power_flow::generate_edges_matrix() {
  add_to_gEdge_all();
  process_shunts();
  
  /** postprocess the diagonal entries */
  for (int i = 0; i < nBuses; ++i) {
    nonZerosBp++;
    nonZerosBpp++;
    string bus_name = bus_id_to_name[i + 1];
    int busType = buses[bus_name].bus_type;
    if (busType == 1) {             // PQ
//      if (abs(gEdge_all[i][3]) < EPS)
//        cout << "warnning: " << bus_name << endl;
//      if (abs(gEdge_all[i][5]) < EPS)
//             cout << "warnning: " << bus_name << endl;
//      gEdge_all[i][5] = gEdge_all[i][5] < EPS ? 0.00001 : gEdge_all[i][5];
      gEdge_all[i][6] = gEdge_all[i][3];
//      assert(abs(gEdge_all[i][3]) > EPS);
    } else if (busType == 2) {      // PV
      gEdge_all[i][6] = 1.;
    } else if (busType == 3) {      // SLACK
      gEdge_all[i][6] = gEdge_all[i][5] = 1.;
    }
  }

#if DEBUG
//  print_matrix<real__t>(gEdge_all, "gEdge_all matrix");
#endif
}

/**
 * gVertex_all contains:
 * 0    1    2   3   4   5   6   7
 * B'p, B"p, ep, Vm, Va, Pn, Qn, busType
 */

void fast_decoupled_power_flow::generate_vertices_matrix() {
  for (int i = 0; i < nBuses; ++i) {
    string  bus_name = bus_id_to_name[i + 1];
    BUS& bus = buses[bus_name];
    uint__t bus_type = bus.bus_type;
    
    gVertex_all[i][2] = (real__t)all_degree[i] + 1.;
    gVertex_all[i][5] = bus.Pgen - bus.Pload;
    gVertex_all[i][6] = bus.Qgen - bus.Qload;
    gVertex_all[i][7] = (real__t)bus_type;
    if (bus_type <= 1) {             // PQ
      gVertex_all[i][0] = (real__t)pqv_degree[i] + 1.;
      gVertex_all[i][1] = (real__t)pq_degree[i] + 1.;
      gVertex_all[i][3] = 1.;
      gVertex_all[i][4] = 0.;
    } else if (bus_type == 2)  {     // PV
      gVertex_all[i][0] = (real__t)pqv_degree[i] + 1.;
      gVertex_all[i][1] = 1.;
      gVertex_all[i][3] = bus.Vm;
      gVertex_all[i][4] = 0.;
    } else if (bus_type == 3) {     // SLACK
      gVertex_all[i][0] = 1.;
      gVertex_all[i][1] = 1.;
      gVertex_all[i][3] = bus.Vm;
      gVertex_all[i][4] = bus.Va * PI / 180;
    }
  }
}

void fast_decoupled_power_flow::run() {
  setup_system();
  generate_edges_matrix();
  generate_vertices_matrix();
  
  for (int i = 0; i < gEdge_all.size(); ++i) {
    if (gEdge_all[i][0] == gEdge_all[i][1] && abs(gEdge_all[i][5]) < EPS )
      cout << bus_id_to_name[i+1] << ", " << gEdge_all[i][0] << ", "
           <<  gEdge_all[i][5] << ", " << gEdge_all[i][6] << endl;
  }
  
  string result = fdpf_LU_factorize(nonZerosBp, nonZerosBpp,
                                    gVertex_all, gEdge_all,
                                    max_change_P, max_change_Q,
                                    maxIter, iterations, output_path);
  
  print_line();
  printf("Fast decoupled power flow result: %s\n", result.c_str());
  print_line();
  printf("\n\n");
  
}

}
