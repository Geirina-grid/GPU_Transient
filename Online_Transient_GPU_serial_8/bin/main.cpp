#include "transient.hpp"
#include "powerflow.hpp"
#include "trans_cuda.hpp"

using namespace transient_analysis;
int main(int argc, char** argv) {
  std::clock_t start = std::clock();
  printLogo();

/** list of data paths:
 *    data_China
 *    data_epri_wscc9
 *    data_epri_wscc9-test
 *    data_epri_36
 *    data_epri_36-test
 *    data_China
 *    data_epri_SC
 *    data_epri_SC_wo_hvdc
 */

  string data_path   = "../data/data_epri_36-test";
  string output_path = "../output";

  fast_decoupled_power_flow fdpf(data_path);
  fdpf.run();

  /** simulation time control parameters */
  real__t start_time = 0.;
  real__t end_time = 20.;
  real__t time_stepping = 0.01;

  Transient_Solver solver(start_time, end_time, time_stepping, data_path);
  solver.run(argc, argv);

  printf("+++Results written to %s\n", output_path.c_str());
  printf("+++End of simulation!\n");
  printf("+++Elapsed time: %.4f seconds\n\n", (std::clock() - start) / (real__t)CLOCKS_PER_SEC);

  return 0;
}
