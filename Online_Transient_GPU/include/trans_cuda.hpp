#ifndef trans_cuda_hpp
#define trans_cuda_hpp

//#include "GraphDyn_util.hpp"
//#include <unordered_map>
#include <map>
//#include <unistd.h>

using namespace std;
//namespace transient_analysis {
struct GENERATOR {
  string bus_name;
  int bus_id;
  int Gen_Model, Gen_Par; // model number, parameter set number
  int AVR_Model, AVR_Par; // model number, parameter set number
  int GOV_Model, GOV_Par; // model number, parameter set number
  int PSS_Model, PSS_Par; // model number, parameter set number
  double Rate_MVA, Rate_MW;
  double  Xdp, Xdpp, X2, TJ;
  double  omega_ref, freq_ref;
  double  Vt_ref, Pe_ref;
  double  Efd0, mu0;
  double  Ep0, Epp0;
};

//namespace transient_analysis {
//extern "C" void test_cuda(unordered_map<int, int> generators);

//void test_cuda(map<string, GENERATOR> &generators);

extern "C" void test_cuda(map<string, GENERATOR> generators);

int cuda_for();

//class Transient_Cuda{
//    public: Transient_Cuda(real__t start_time, real__t end_time, real__t time_stepping, string main_path);
//    public: ~Transient_Cuda();
//    public: void run(int argc, char** argv);
//};

#endif
