#include "JUSQUACE.hpp"
#include "basis_subroutines.hpp"
#include "chip_library.hpp"

int main(int argc, char *argv[])
{
  cout.precision(10);
  const unsigned int width=30;
  bool info=false,tar_matrix=false,tar_state=false,tar_spec=true;
  string option;
  if(info==true)
  {
    cout << "Please enter the program function you would like to use (input_matrix/target_state/spectroscopy): ";
    cin >> option;
    if(option=="input_matrix"){tar_matrix=true;}
    else if(option=="target_state"){tar_state=true;}
    else if(option=="spectroscopy"){tar_spec=true;}
    else {cout << "The input does not match any of the functions available";exit(1);}
  }

  JUSQUACE_PFA chip;
  ArchitectureII(chip);//load chip data into the JUSQUACE_PFA class and precompute some constants
  chip.N_ctl0=1;//number of flux-drive terms
  chip.N_ctl1=0;//number of charge-drive terms

  clock_t cpu_startTime, cpu_endTime;//elapsed runtime variables
  double cpu_ElapseTime=0;//same as last line
  const unsigned int tar_idx=0;//state vector index of the intialal basis state
  const unsigned int end=4;//number of prob. amps. printed
  const unsigned int dim=1 << chip.num_qubits;//dimension of the qubit-subspace matrix
  const unsigned int offset=2*dim;//offset for for matrix rows (flat row-major storage format)

  double *state=new double[2*chip.state_size]();//state vector data strucutre
  memset(state,0,2*chip.state_size*sizeof(double));
  double *input_matrix= new double[2*dim*dim]();//qubit-subspace matrix
  memset(input_matrix,0,2*dim*dim*sizeof(double));
  vector<uint64_t> TAR_IDX={0,16,32,64};//state vector indices for the prob. amp. we intend to store

  if(info==1)
  {

    cout << "#State vector " << chip.state_size << "[Elements]" << endl;

    cout << "#State vector " << (double)(chip.state_size*16)/((double) pow(2,10)) << "[KiB]" << endl;

    cpu_startTime = clock();
  }
  if(tar_matrix==true)
  {
    for(unsigned int i=0;i<dim;i++)
    {

      memset(state,0,2*chip.state_size*sizeof(double));
      state[2*chip.idx_state[i]+0]=1;
      chip.pfa_time_evolution(state);
      for(unsigned int j=0;j<dim;j++)
      {
        input_matrix[2*i+0+offset*j]=state[2*chip.idx_state[j]+0];
        input_matrix[2*i+1+offset*j]=state[2*chip.idx_state[j]+1];
      }
    }
  }
  else if (tar_state==true)
  {

    memset(state,0,2*chip.state_size*sizeof(double));
    state[2*chip.idx_state[tar_idx]+0]=1;

    cerr << setw(width) << left << "#Drive Frequency/2pi [GHz]";
    cerr << setw(width) << left << "Time [ns]";
    for(unsigned int i=0;i<end;i++){cerr << setw(width) << left << "P("+to_string(TAR_IDX[i])+")";}
    cerr << "\n";

    chip.pfa_time_evolution(state);

    cerr << setw(width) << left << chip.para_ctl0[0][1]/(2*pi);
    cerr << setw(width) << left << (chip.T_end*chip.tau);
    for(unsigned int i=0;i<end;i++)
    {
      double hold=(state[2*TAR_IDX[i]+0]*state[2*TAR_IDX[i]+0]+state[2*TAR_IDX[i]+1]*state[2*TAR_IDX[i]+1]) ;
      cerr << setw(width) << left << hold;
    }
    cerr << "\n";

  }
  else if (tar_spec==true)
  {

    double hold=chip.para_ctl0[0][1];
    cerr << setw(width) << left << "#Drive Frequency/2pi [GHz]";
    cerr << setw(width) << left << "Time [ns]";
    for(unsigned int i=0;i<end;i++)
    {
      //TAR_IDX[i]=chip.idx_state[i];
      cerr << setw(width) << left << "P("+to_string(TAR_IDX[i])+")";
    }
    cerr << "\n";

    for(int i=-10;i<=10;i++)
    {
      chip.para_ctl0[0][1]=hold+2*pi*0.001*i;
      memset(state,0,2*chip.state_size*sizeof(double));
      state[2*tar_idx+0]=1;
      chip.pfa_time_evolution(state,TAR_IDX);
      cout <<  "\n";
    }
  }
  if(info==1)
  {
    cpu_endTime = clock();

    cpu_ElapseTime = ((cpu_endTime - cpu_startTime)/CLOCKS_PER_SEC);
    cout << "#Runtime [s] = " << cpu_ElapseTime << endl;
    cout << "#Simulation time [ns] = " << (chip.T_end*chip.tau) << endl;
  }
  delete[] state;
  delete[] input_matrix;

  return 0;
}
