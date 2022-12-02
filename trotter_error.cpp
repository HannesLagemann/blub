#include "JUSQUACE.hpp"
#include "basis_subroutines.hpp"
#include "chip_library.hpp"

inline double euclidean_distance(double *state0,double *state1,unsigned int size)
{
  double diff_real,diff_img;
  double error=0;
  for(unsigned int i=0;i<size;i++)
  {
    diff_real =state0[2*i] - state1[2*i];
    diff_img  =state0[2*i+1] - state1[2*i+1];
    error=error+diff_real*diff_real+diff_img*diff_img;
  }
  error=sqrt(error);
  return error;
}

inline double error(double tau, unsigned int T,JUSQUACE_PFA &pfa_data,JUSQUACE_FD &fd_data)
{
  ArchitectureII(pfa_data,tau);
  ArchitectureII(fd_data,tau);

  pfa_data.N_ctl0=1;
  pfa_data.N_ctl1=1;
  fd_data.N_ctl0=pfa_data.N_ctl0;
  fd_data.N_ctl1=pfa_data.N_ctl1;

  double err=0;
  double max=0;
  unsigned int size=0;
  if(pfa_data.state_size==fd_data.state_size){size=fd_data.state_size;}
  else{exit(0);}

  double *state_trans_pfa= new double[2*size]();
  double *state_trans_fd= new double[2*size]();

  pfa_data.T_start=T;
  pfa_data.T_end=T+1;
  fd_data.T_start=T;
  fd_data.T_end=T+1;

  for(unsigned int i=0;i<size;i++)
  {

    set_state_zero(state_trans_pfa,size);
    set_state_zero(state_trans_fd,size);

    state_trans_pfa[2*i+0]=1;
    state_trans_fd[2*i+0]=1;

    pfa_data.time_offset_ctl0=0;
    pfa_data.time_offset_ctl1=0;

    pfa_data.pfa_time_evolution(state_trans_pfa);
    fd_data.fd_time_evolution(state_trans_fd);

    err=euclidean_distance(state_trans_pfa,state_trans_fd,size);
    if(err>max){max=err;}
  }

  delete[] state_trans_pfa;
  delete[] state_trans_fd;
  return max;
}

int main(int argc, char *argv[])
{

  cout.precision(15);
  const unsigned int width=25;
  unsigned int T=25;
  double tau=1.0;
  double err=0;

  JUSQUACE_PFA pfa_data;
  JUSQUACE_FD fd_data;

  cout << setw(width) << left << "#Time";
  cout << setw(width) << left << "Time grid parameter" ;
  cout << setw(width) << left << "Iterations" ;
  cout << setw(width) << left << "Euclidean distance (Pfa vs. Exact)";
  cout << "\n";

  unsigned int num_divisions=7;
  for(unsigned int i=0;i<num_divisions;i++)
  {
    tau=tau/10.0;
    T=T*10;
    err=error(tau,T,pfa_data,fd_data);
    cout << setw(width) << left << tau*T ;
    cout << setw(width) << left << tau ;
    cout << setw(width) << left << T ;
    cout << setw(width) << left << err ;
    cout << "\n";
  }
  return 0;
}
