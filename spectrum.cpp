#include "JUSQUACE.hpp"
#include "basis_subroutines.hpp"
#include "chip_library.hpp"

void compute_energie_labels(JUSQUACE_FD &chip,string *eigenvalues_label,const string gnuplot_lab,const unsigned int num_levels)
{

  //This sub. computes the labels for the eigenenergies of the non-interacting Hamiltoniian.
  //First, we extract the ket-vector subindices q_0,q_1, ... form the kets |...,q_1,q_0>.
  //Next, we compute the corresponding energies_z.
  //Then, we sort the energies together with the corresponding labels.

  string lab="";
  string lab_end=")'";
  string lab_min="";

  unsigned int size=chip.state_size;
  unsigned int hold_ui=0;
  unsigned int offset=0;
  unsigned int idx_min=0;
  unsigned int num_transmons=chip.Nt;
  unsigned int num_resonators=chip.Nk;
  unsigned int num_system_elements=num_transmons+num_resonators;
  unsigned int num_basis_states=4;
  u_int64_t *id=new u_int64_t[num_system_elements]();

  double ev=0;
  double min=0;
  double *eigenvalues=new double[size];




  //Compute energies and labels
  for(u_int64_t n=0;n<size;n++)
  {

    hold_ui=n;
    lab="";
    ev=0;



    for(u_int64_t i=0;i<num_system_elements;i++)
    {
      unsigned int fac_i=pow(num_basis_states,i);
      unsigned int fac_i_p_1=pow(num_basis_states,i+1);
      id[i]= hold_ui % fac_i_p_1;
      id[i]= id[i]/fac_i;
      hold_ui=hold_ui - id[i]*fac_i;
    }


    for(unsigned int i=0;i<num_resonators;i++)
    {
      unsigned int idx=num_resonators-i-1;
      lab+=to_string(id[idx+num_transmons])+",";
      ev+= chip.res_freq[idx]*id[idx+num_transmons];


    }

    for(unsigned int i=0;i<num_transmons;i++)
    {
      unsigned int idx=num_transmons-i-1;
      lab+=to_string(id[idx+offset])+",";
      //ev+= (EVS_data.eigenvalues_qubit[idx][id[idx+offset]]-EVS_data.eigenvalues_qubit[idx][0]);
      ev+=chip.freq_at_flux_zero[idx]*id[idx]+(chip.anharm_at_flux_zero[idx]*0.5)*id[idx]*(id[idx]-1);
      //ev+=chip.flux_freq(chip.freq_at_flux_zero[idx],0,chip.flux_offset[idx],chip.asymmetry_factor_qubit[idx])*id[idx]+(chip.anharm_at_flux_zero[idx]*0.5)*id[idx]*(id[idx]-1);

    }


    string hold =lab;
    eigenvalues_label[n]=gnuplot_lab+hold.substr(0, hold.size()-1)+lab_end;
    eigenvalues[n]=ev;
    //cerr << n << "\t" << (ev/(2*pi)) << "\t" << eigenvalues_label[n] << endl;
  }

  //Sort energies and labels
  for(unsigned int i=0;i<size;i++)
  {
    min=10000000;
    for(unsigned int j=i;j<(size-i);j++)
    {
      if(eigenvalues[j]<min)
      {
        min=eigenvalues[j];
        idx_min=j;
      }
    }
    //switch i and idx_min
    double hold_ev=eigenvalues[i];
    string hold_str=eigenvalues_label[i];
    eigenvalues[i]=eigenvalues[idx_min];
    eigenvalues_label[i]=eigenvalues_label[idx_min];
    eigenvalues[idx_min]=hold_ev;
    eigenvalues_label[idx_min]=hold_str;
  }

  delete [] eigenvalues;
  delete [] id;
}

int main(int argc, char *argv[])
{
  const unsigned int num_iterations=10;
  const unsigned int num_energies=10;
  const unsigned int width=25;
  const unsigned int tar_idx=1;
  const double delta=0.05;
  double flux=0;
  JUSQUACE_FD fd_data;
  ArchitectureII(fd_data);
  fd_data.N_ctl0=0;
  fd_data.N_ctl1=0;
  string *eigenvalues_label =new string[fd_data.state_size]();
  string filename="spec.dat";
  ofstream ofile;

  ofile.open(filename, ios_base::out);
  compute_energie_labels(fd_data,eigenvalues_label,"'~z{.4-} = (",num_energies);

  // Header line with variable name (here FLUX) and energie labels
  ofile << setw(width) << left << "#In this file we store the discrete energies_{z}(FLUX) in units of GHz\n";
  ofile << setw(width) << left << "#'FLUX/2pi'";
  for(unsigned int i=0;i<num_energies;i++){ofile << setw(width) << left << eigenvalues_label[i];}
  ofile << "\n";

  //Computaiton and storage of energies E_i(x) as a function of the variable x (here FLUX).
  for(unsigned int i=0;i<num_iterations;i++)
  {
    //Computaiton of energies E_i(x) as a function of the variable x (here FLUX).
    flux=pi*i*delta;
    fd_data.flux_offset[tar_idx]=flux;
    fd_data.update_flux_dependent_data();
    fd_data.initial_fd_data_structures();
    fd_data.update_time_dependent_data(0);
    fd_data.diagonalization(0,'N');

    //Storage of energies E_i(x) as a function of the variable x (here FLUX).
    ofile << setw(width) << left << flux/(pi);
    for(unsigned int j=0;j<num_energies;j++){ofile << setw(width) << left << (fd_data.eigenvalues_matrix[j]-fd_data.eigenvalues_matrix[0])/(2*pi);}
    ofile << "\n";
  }

  // Free dynamic memory and close data file
  delete [] eigenvalues_label;
  ofile.close();
  return 0;
}
