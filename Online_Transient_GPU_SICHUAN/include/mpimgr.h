#ifndef MPIMGR_H
#define MPIMGR_H

#include <mpi.h>
#include <fstream>

using namespace dealii;

class MPIMGR
{    
public:
  MPIMGR(int argc, char **argv);
  ~MPIMGR();
  
  int get_nb_processes(){ return nb_processes;  };
  int get_this_process(){ return this_process;  };
  
  MPI_Comm get_communicator(){ return mpi_communicator; };
  
  void write(std::string str) { pcout<<str<<std::endl; }
  void write_in_line(std::string str) { pcout<<str; }
  
  bool is_this_my_job(int quad_no) const;
  int get_max_job_per_proc(){ return job_per_proc;};
  
  ConditionalOStream get_ostream() {return pcout; };
  
  void distribute_jobs(int nb_quad);
  
private:       
  int nb_processes;
  int this_process;
  
  MPI_Comm mpi_communicator;
  ConditionalOStream pcout;
  
  int job_per_proc;
  
  Utilities::MPI::MPI_InitFinalize mpi_initialzation;
};

#endif
