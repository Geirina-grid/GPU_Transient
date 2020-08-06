#include "ode_solver.cuh"

namespace transient_analysis {

void ODE_solver::apply_perturbation(real__t t, EPRI_GEN_DATA& gen) {
  /*
  if (t >= 10 && t <= 15 && (gen.bus_id == 0)) {
    Vm_ref = parameters.Vt_ref * 1.0;
  } else {
    Vm_ref = parameters.Vt_ref;
  }

  if (t >= 5 && t <= 20 && (gen.bus_id == 2 || gen.bus_id == 3)) {
    omega_ref = parameters.omega_ref * 1.0;
  } else {
    omega_ref = parameters.omega_ref;
  }
   */
   printf("Here in apply_perturbation running on CPU.");

//  if (t >= 5 && t <= 20 && (gen.bus_id == 2 || gen.bus_id == 3)) {
//    gov.mu_zero = gov.mu_zero_base * 1.0;
//  } else {
//    gov.mu_zero = gov.mu_zero_base;
//  }
}

}  // namespace transient_analysis
