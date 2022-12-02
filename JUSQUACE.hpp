#include "pulse_library.hpp"

#ifndef DEV_DATA
#define DEV_DATA

class DEVICE_DATA
{
  protected:
    bool device_memory=false;//false=free and true=memory needs to be allocated
    static const bool bool_eigenvalues_with_high_accuracy=1;
    static const unsigned int Nts=4;
    static const unsigned int Nks=4;
    static const unsigned int max_num_flux_drive=1;
    static const unsigned int max_num_charge_drive=1;

  public:
    // Device system and interaction numbers (sizes)
    unsigned int Nt=0;
    unsigned int Nk=0;
    unsigned int Ni=0;
    unsigned int Ni_td=0;
    unsigned int N_ctl0=0;
    unsigned int N_ctl1=0;
    unsigned int T_start=0;
    unsigned int T_end=0;
    u_int64_t state_size=0;
    double tau=0;
    double time_offset_ctl0=0;
    double time_offset_ctl1=0;
    double time_total=0;

    // The following dynamic arrays (heap memory) containe the device parameters and some auxiliary parametrs,see amp_factor.
    double *EC=nullptr; // Capacitive energies
    double *EJ=nullptr; // Josephson energies
    double *freq_at_flux_zero=nullptr; //  device parameters
    double *anharm_at_flux_zero=nullptr; // device parameters
    double *asymmetry_factor_qubit=nullptr; // device parameters
    double *flux_offset=nullptr; // Flux at time zero
    double *para_Ni=nullptr;// Time-independent interacting strength parameters
    double *para_td_Ni=nullptr;// Time-dependent interacting strength parameters
    double *amp_factor=nullptr; // Precomputed constant factors which are device parameter dependent
    double *res_freq=nullptr; //  device parameters

  protected:
    // The following dynamic arrays (heap memory) contain indices
    // to control the time-dependent interacting strength g_{i,j}(t) of the model.
    unsigned int *idx_Ni=nullptr;
    unsigned int *idx_td_Ni=nullptr;
    unsigned int *idx_td_Ni_ctl0=nullptr;

  public:
    vector<unsigned int> idx_ctl0;
    vector<vector<double>> para_ctl0;

    vector<unsigned int> idx_ctl1;
    vector<vector<double>> para_ctl1;

    // The following vector data structures are used to store the pulse parameters of a device (system)
    // The idx_ctl data structures contain the target indices (index of the system element |n_{1}->1,n_{0}->0>)
    // The para_ctl data structures contain the actual pulse parameters, see charge_ctl/1 and charge_flux/0.
    vector<double> lamda_trans_shifted;
    vector<vector<double>> phase_CZ;

    string chip_name;
    unsigned int num_qubits=0;
    u_int64_t *idx_state=nullptr;

  protected:
    //***A set of subroutines which compute physical real-valued quantities as functions of the time and/or the dimensionless flux***//
    //The subroutine computes dimensionless charge values as a function of time dt for the controle pulse parameters ctl_parameters
    inline double charge_ctl(const vector<double> &ctl_parameters,const double dt)
    {
      const double offset=0;
      return single_qubit_charge_ctl(ctl_parameters,dt,offset);
    }
    //The subroutine computes dimensionless flux values as a function of time dt for the controle pulse parameters ctl_parameters
    inline double flux_ctl(const vector<double> &ctl_parameters,const double dt)
    {

      return two_qubit_flux_ctl(ctl_parameters,dt);
    }
    //The subroutine computes the dimensionless of the dimensionless flux values as a function of time dt for the controle pulse parameters ctl_parameters
    inline double flux_ctl_der(const vector<double> &ctl_parameters,const double dt)
    {
      return two_qubit_flux_ctl_der(ctl_parameters,dt);
    }
    //The subroutine is used to compute the low accuracy eigenvalues E_n(flux_t+flux_t0)=flux_freq(flux_t+flux_t0)*n+(anharmonicity/2)*n*(n-1)
    inline double flux_freq(double qubit_freq_flux_t0,double flux_t,double flux_t0,double asymmetry_factor)
    {

      double hold=0;
      hold+=(cos((flux_t+flux_t0))*cos((flux_t+flux_t0)));
      hold+=(asymmetry_factor*asymmetry_factor*sin((flux_t+flux_t0))*sin((flux_t+flux_t0)));
      hold=sqrt(sqrt(hold));
      hold=qubit_freq_flux_t0*(hold);
      return hold;
    }
    //The subroutine is used to compute the high accuracy eigenvalues E_n(flux_t+flux_t0)=flux_freq_with_high_accuracy(flux_t+flux_t0)*n+(flux_freq_anharmonicity_with_high_accuracy/2)*n*(n-1)
    inline double flux_freq_with_high_accuracy(double EC,double EJ1,double EJ2,double flux_t,double flux_t0)
    {
      double EJ=sqrt(EJ1*EJ1+EJ2*EJ2+2*EJ1*EJ2*cos(2*(flux_t+flux_t0)));
      double epsi=sqrt(EC/(2*EJ));
      double hold_esp=0;
      double hold_0=sqrt(2*EC*EJ);
      double hold_1=0;


      hold_esp=1.0;
      hold_1+=1.0000000000000000*hold_esp;//0 th

      hold_esp*=epsi;
      hold_1+=0.2500000000000000*hold_esp;//1 th

      hold_esp*=epsi;
      hold_1+=0.1640625000000000*hold_esp;//2 th

      hold_esp*=epsi;
      hold_1+=0.1484375000000000*hold_esp;//3 th

      hold_esp*=epsi;
      hold_1+=0.1623229980468750*hold_esp;//4 th

      hold_esp*=epsi;
      hold_1+=0.2029113769531250*hold_esp;//5 th

      hold_esp*=epsi;
      hold_1+=0.2814724445343018*hold_esp;//6 th

      hold_esp*=epsi;
      hold_1+=0.4256124496459961*hold_esp;//7 th

      hold_esp*=epsi;
      hold_1+=0.6934342137537897*hold_esp;//8 th

      hold_esp*=epsi;
      hold_1+=1.207704475149512*hold_esp;//9 th

      hold_esp*=epsi;
      hold_1+=2.235748448943923*hold_esp;//10 th

      hold_esp*=epsi;
      hold_1+=4.381048442388419*hold_esp;//11 th

      hold_esp*=epsi;
      hold_1+=9.058167325550627*hold_esp;//12 th

      hold_esp*=epsi;
      hold_1+=19.71115159176793*hold_esp;//13 th

      hold_esp*=epsi;
      hold_1+=45.04936094435116*hold_esp;//14 th

      hold_esp*=epsi;
      hold_1+=107.9430630773788*hold_esp;//15 th

      hold_esp*=epsi;
      hold_1+=270.7321342244121*hold_esp;//16 th

      hold_esp*=epsi;
      hold_1+=709.7234169531544*hold_esp;//17 th

      hold_esp*=epsi;
      hold_1+=1941.981675726763*hold_esp;//18 th

      hold_esp*=epsi;
      hold_1+=5539.060582830617*hold_esp;//19 th

      hold_esp*=epsi;
      hold_1+=16447.89461688524*hold_esp;//20 th

      hold_esp*=epsi;
      hold_1+=50784.73015187223*hold_esp;//21 th,

      hold_esp*=epsi;
      hold_1+=162849.5516964187*hold_esp;//22 th

      hold_esp*=epsi;
      hold_1+=541713.1767775880*hold_esp;//23 th

      hold_esp*=epsi;
      hold_1+=1.867231218078072e+6*hold_esp;//24 th

      hold_1=-0.25*EC*hold_1;

      hold_0=hold_0+hold_1;
      return hold_0;
    }
    //The subroutine is used to compute the high accuracy eigenvalues E_n(flux_t+flux_t0)=flux_freq_with_high_accuracy(flux_t+flux_t0)*n+(flux_freq_anharmonicity_with_high_accuracy/2)*n*(n-1)
    inline double flux_freq_anharmonicity_with_high_accuracy(double EC,double EJ1,double EJ2,double flux_t,double flux_t0)
    {

      double EJ=sqrt(EJ1*EJ1+EJ2*EJ2+2*EJ1*EJ2*cos(2*(flux_t+flux_t0)));
      double epsi=sqrt(EC/(2*EJ));
      double hold_esp=0;
      double hold_1=0;

      hold_esp=1.0;
      hold_1+=1.0000000000000000000000000*hold_esp;//0 th

      hold_esp*=epsi;
      hold_1+=0.5625000000000000*hold_esp;//1 th

      hold_esp*=epsi;
      hold_1+=0.6328125000000000*hold_esp;//2 th

      hold_esp*=epsi;
      hold_1+=0.8898925781250000*hold_esp;//3 th

      hold_esp*=epsi;
      hold_1+=1.431243896484375*hold_esp;//4 th

      hold_esp*=epsi;
      hold_1+=2.535112380981445*hold_esp;//5 th

      hold_esp*=epsi;
      hold_1+=4.844990015029907*hold_esp;//6 th

      hold_esp*=epsi;
      hold_1+=9.865587104111910*hold_esp;//7 th

      hold_esp*=epsi;
      hold_1+=21.22477681143209*hold_esp;//8 th

      hold_esp*=epsi;
      hold_1+=47.96277078389539*hold_esp;//9 th

      hold_esp*=epsi;
      hold_1+=113.3606107500927*hold_esp;//10 th

      hold_esp*=epsi;
      hold_1+=279.3520460809879*hold_esp;//11 th

      hold_esp*=epsi;
      hold_1+=716.0551841267438*hold_esp;//12 th

      hold_esp*=epsi;
      hold_1+=1905.741328339856*hold_esp;//13 th

      hold_esp*=epsi;
      hold_1+=5258.978031955065*hold_esp;//14 th

      hold_esp*=epsi;
      hold_1+=15030.98000512316*hold_esp;//15 th

      hold_esp*=epsi;
      hold_1+=44457.57052814981*hold_esp;//16 th

      hold_esp*=epsi;
      hold_1+=135977.8189939008*hold_esp;//17 th

      hold_esp*=epsi;
      hold_1+=429824.8237530430*hold_esp;//18 th

      hold_esp*=epsi;
      hold_1+=1.403401772681732e+6*hold_esp;//19 th

      hold_esp*=epsi;
      hold_1+=4.730659678149439e+6*hold_esp;//20 th

      hold_esp*=epsi;
      hold_1+=1.645505044661283e+7*hold_esp;//21 th

      hold_esp*=epsi;
      hold_1+=5.903415332047901e+7*hold_esp;//22 th

      hold_esp*=epsi;
      hold_1+=2.183327609822893e+8*hold_esp;//23 th

      hold_esp*=epsi;
      hold_1+=8.320026932974098e+8*hold_esp;//24 th


      hold_1=-0.25*EC*hold_1;

      return hold_1;
    }
    //The subroutine is used to the real-valued interaction strength as a function of flux(time).
    inline double time_dependent_interaction_strength(const double G,const double flux,const unsigned int i)
    {

      //flux should be flux plus offset!
      double hold=EJ[2*i+0]*EJ[2*i+0]+EJ[2*i+1]*EJ[2*i+1]+2*EJ[2*i+0]*EJ[2*i+1]*cos(2*flux);
      hold=sqrt(hold);
      hold=G*sqrt(sqrt(hold/(8*EC[i])));
      return hold;
    }

  public:

    inline void update_flux_dependent_data(void)
    {
      for(unsigned int i=0;i<Nt;i++)
      {
        if(bool_eigenvalues_with_high_accuracy==1)
        {
          freq_at_flux_zero[i]=flux_freq_with_high_accuracy(EC[i],EJ[2*i+0],EJ[2*i+1],0,flux_offset[i]);
          anharm_at_flux_zero[i]=flux_freq_anharmonicity_with_high_accuracy(EC[i],EJ[2*i+0],EJ[2*i+1],0,flux_offset[i]);
        }
        else if (bool_eigenvalues_with_high_accuracy==0)
        {
          freq_at_flux_zero[i]=flux_freq(freq_at_flux_zero[i],0,flux_offset[i],asymmetry_factor_qubit[i]);
          anharm_at_flux_zero[i]=anharm_at_flux_zero[i];
        }
      }
    }

    DEVICE_DATA()=default;
    DEVICE_DATA(vector<double> qubit_freq, vector<double> anharm,vector<double> flux_offset_in,vector<double> asymmetry_factor,vector<double> resonator_freq ,unsigned int num_transmons,unsigned int num_resonators,
      unsigned int num_flux_drives,unsigned int num_charge_drives,
      unsigned int num_interactions,vector<unsigned int> interaction_idx,vector<double> interaction_para,
      double time_grid_parameter,unsigned int num_iter_steps,
      unsigned int num_td_interactions=0,vector<unsigned int> interaction_td_idx={},vector<unsigned int> interaction_td_idx_ctl0={},vector<double> interaction_td_para={},
      vector<double> EC_input={},vector<double> EJ_input={},unsigned int num_qubit=0)
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not allocated" << endl;}};
      if(visible==true){cerr << "DD constructor" << endl;}

      Nt=num_transmons;
      Nk=num_resonators;
      Ni=num_interactions;
      num_qubits=num_qubit;
      Ni_td=num_td_interactions;
      N_ctl0=num_flux_drives;
      N_ctl1=num_charge_drives;
      T_start=0;
      T_end=num_iter_steps;
      tau=time_grid_parameter;
      state_size=1 << (2*(Nt+Nk));
      num_qubits=num_qubit;


      const string option="allocate";
      if(option=="allocate")
      {
        if(N_ctl0!=0)
        {

          idx_ctl0.resize(N_ctl0,0);
          para_ctl0.resize(N_ctl0,vector<double>(7,1));
          phase_CZ.resize(Nt,vector<double>(num_qubits,0));
        }
        else
        {
          warning("charge control pulses");
        }

        if(N_ctl1!=0)
        {

          idx_ctl1.resize(N_ctl1,0);
          para_ctl1.resize(N_ctl1,vector<double>(7,1));
        }
        else
        {
          warning("flux control pulses");
        }

        if(num_qubits!=0)
        {
          idx_state= new u_int64_t[(1 << num_qubits)]();
          lamda_trans_shifted.resize(num_qubits,0);
        }
        else
        {
          warning("computatoinal states indices");
        }

        if(Nt!=0)
        {

          EC = new double[Nt]();
          EJ = new double[2*Nt]();
          freq_at_flux_zero = new double[Nt]();
          anharm_at_flux_zero = new double[Nt]();
          asymmetry_factor_qubit = new double[Nt]();
          flux_offset = new double[Nt]();
          amp_factor= new double[2*Nt]();
        }
        else
        {
          warning("anharmonic osc. device parameters");
        }

        if(Nk!=0)
        {
          res_freq=new double[Nk]();
        }
        else
        {
          warning("harmonic osc. device parameters");
        }


        if(Ni!=0)
        {
          idx_Ni = new unsigned int[2*Ni]();
          para_Ni = new double[Ni]();
        }
        else
        {
          warning("time-independent interactions");
        }

        if(Ni_td!=0)
        {
          idx_td_Ni = new unsigned int[2*Ni_td]();
          idx_td_Ni_ctl0 = new unsigned int[Ni_td]();
          para_td_Ni = new double[Ni_td]();
        }
        else
        {
          warning("time-dependent interactions");
        }
      }

      for(unsigned int i=0;i<Nt;i++)
      {
        if((bool_eigenvalues_with_high_accuracy==1) &&(EC_input.size()*EJ_input.size())!=0)
        {
          EC[i]=EC_input[i];
          EJ[2*i+0]=EJ_input[2*i+0];
          EJ[2*i+1]=EJ_input[2*i+1];
          double Esum=(EJ[2*i+0]+EJ[2*i+1]);
          asymmetry_factor_qubit[i]=(EJ[2*i+1]-EJ[2*i+0])/Esum;
          amp_factor[2*i+0]=-(asymmetry_factor_qubit[i]/(2.0*sqrt(2.0)))*sqrt(sqrt(Esum/(2*EC[i])));
          amp_factor[2*i+1]=(asymmetry_factor_qubit[i]*asymmetry_factor_qubit[i]-1.0)/32.0;
          flux_offset[i]=flux_offset_in[i];
          freq_at_flux_zero[i]=flux_freq_with_high_accuracy(EC[i],EJ[2*i+0],EJ[2*i+1],0,flux_offset[i]);
          anharm_at_flux_zero[i]=flux_freq_anharmonicity_with_high_accuracy(EC[i],EJ[2*i+0],EJ[2*i+1],0,flux_offset[i]);
        }
        else if (bool_eigenvalues_with_high_accuracy==0)
        {
          flux_offset[i]=flux_offset_in[i];
          asymmetry_factor_qubit[i]=asymmetry_factor[i];
          freq_at_flux_zero[i]=qubit_freq[i];
          anharm_at_flux_zero[i]=anharm[i];
        }
      }
      for(unsigned int i=0;i<Nk;i++){res_freq[i]=resonator_freq[i];}
      for(unsigned int i=0;i<(2*Ni);i++){idx_Ni[i]=interaction_idx[i];}
      for(unsigned int i=0;i<Ni;i++){para_Ni[i]=interaction_para[i];}
      for(unsigned int i=0;i<(2*Ni_td);i++){idx_td_Ni[i]=interaction_td_idx[i];}
      for(unsigned int i=0;i<Ni_td;i++){para_td_Ni[i]=interaction_td_para[i];}
      for(unsigned int i=0;i<Ni_td;i++){idx_td_Ni_ctl0[i]=interaction_td_idx_ctl0[i];}
    }

    ~DEVICE_DATA()
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not freed" << endl;}};
      if(visible==true){cerr << "DD destructor" << endl;}

      const string option="free";
      if(option=="free")
      {

        if(N_ctl0!=0)
        {

          vector<unsigned int>().swap(idx_ctl0);
          vector<vector<double>>().swap(para_ctl0);
          vector<vector<double>>().swap(phase_CZ);
        }
        else
        {
          warning("charge control pulses");
        }


        if( N_ctl1!=0)
        {

          vector<unsigned int>().swap(idx_ctl1);
          vector<vector<double>>().swap(para_ctl1);
        }
        else
        {
          warning("flux control pulses");
        }

        if(num_qubits!=0)
        {
          delete [] idx_state;
          idx_state=nullptr;
          vector<double>().swap(lamda_trans_shifted);
        }
        else
        {
          warning("computatoinal states indices");
        }

        if (Nt!=0 && option=="free")
        {

          delete [] EC;
          EC=nullptr;
          delete [] EJ;
          EJ=nullptr;
          delete [] freq_at_flux_zero;
          freq_at_flux_zero=nullptr;
          delete [] anharm_at_flux_zero;
          anharm_at_flux_zero=nullptr;
          delete [] asymmetry_factor_qubit;
          anharm_at_flux_zero=nullptr;
          delete [] flux_offset;
          flux_offset=nullptr;
          delete [] amp_factor;
          amp_factor=nullptr;
        }
        else
        {
          warning("anharmonic osc. device parameters");
        }

        if(Nk!=0)
        {

          delete [] res_freq;
          res_freq=nullptr;
        }
        else
        {
          warning("harmonic osc. device parameters");
        }

        if(Ni !=0)
        {
          delete [] idx_Ni;
          idx_Ni=nullptr;
          delete [] para_Ni;
          para_Ni=nullptr;
        }
        else
        {
          warning("time-independent interactions");
        }

        if(Ni_td!=0)
        {
          delete [] idx_td_Ni;
          idx_td_Ni=nullptr;
          delete [] idx_td_Ni_ctl0;
          idx_td_Ni_ctl0=nullptr;
          delete [] para_td_Ni;
          para_td_Ni=nullptr;
        }
        else
        {
          warning("time-dependent interactions");
        }
      }
    }

    DEVICE_DATA(const DEVICE_DATA & other)
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not copied" << endl;}};
      if(visible==true){cerr << "DD copy constructor" << endl;}

      Nt=other.Nt;
      Nk=other.Nk;
      Ni=other.Ni;
      num_qubits=other.num_qubits;
      Ni_td=other.Ni_td;
      N_ctl0=other.N_ctl0;
      N_ctl1=other.N_ctl1;
      T_start=other.T_start;
      T_end=other.T_end;
      tau=other.tau;
      state_size=other.state_size;
      num_qubits=other.num_qubits;

      const string option="copy";
      if(option=="copy")
      {
        if(N_ctl0!=0)
        {

          idx_ctl0.resize(N_ctl0,0);
          copy(other.idx_ctl0.begin(), other.idx_ctl0.end(),back_inserter(idx_ctl0));
          para_ctl0.resize(N_ctl0,vector<double>(7,1));
          copy(other.para_ctl0.begin(), other.para_ctl0.end(),back_inserter(para_ctl0));
          phase_CZ.resize(Nt,vector<double>(num_qubits,0));
          copy(other.phase_CZ.begin(), other.phase_CZ.end(),back_inserter(phase_CZ));
        }
        else
        {
          warning("charge control pulses");
        }

        if(N_ctl1!=0)
        {

          idx_ctl1.resize(N_ctl1,0);
          copy(other.idx_ctl1.begin(), other.idx_ctl1.end(),back_inserter(idx_ctl1));
          para_ctl1.resize(N_ctl1,vector<double>(7,1));
          copy(other.para_ctl1.begin(), other.para_ctl1.end(),back_inserter(para_ctl1));
        }
        else
        {
          warning("flux control pulses");
        }

        if(num_qubits!=0)
        {
          idx_state= new u_int64_t[(1 << num_qubits)]();
          memcpy(idx_state,other.idx_state,sizeof(u_int64_t)*(1 << num_qubits));
          lamda_trans_shifted.resize(num_qubits,0);
          copy(other.lamda_trans_shifted.begin(), other.lamda_trans_shifted.end(),back_inserter(lamda_trans_shifted));
        }
        else
        {
          warning("computatoinal states indices");
        }

        if(Nt!=0)
        {

          EC = new double[Nt]();
          memcpy(EC,other.EC,sizeof(double)*1*Nt);
          EJ = new double[2*Nt]();
          memcpy(EJ,other.EJ,sizeof(double)*2*Nt);
          freq_at_flux_zero = new double[Nt]();
          memcpy(freq_at_flux_zero,other.freq_at_flux_zero,sizeof(double)*1*Nt);
          anharm_at_flux_zero = new double[Nt]();
          memcpy(anharm_at_flux_zero,other.anharm_at_flux_zero,sizeof(double)*1*Nt);
          asymmetry_factor_qubit = new double[Nt]();
          memcpy(asymmetry_factor_qubit,other.asymmetry_factor_qubit,sizeof(double)*1*Nt);
          flux_offset = new double[Nt]();
          memcpy(flux_offset,other.flux_offset,sizeof(double)*1*Nt);
          amp_factor= new double[2*Nt]();
          memcpy(amp_factor,other.amp_factor,sizeof(double)*2*Nt);
        }
        else
        {
          warning("anharmonic osc. device parameters");
        }

        if(Nk!=0)
        {
          res_freq=new double[Nk]();
          memcpy(res_freq,other.res_freq,sizeof(double)*1*Nk);
        }
        else
        {
          warning("harmonic osc. device parameters");
        }


        if(Ni!=0)
        {
          idx_Ni = new unsigned int[2*Ni]();
          memcpy(idx_Ni,other.idx_Ni,sizeof(unsigned int)*2*Ni);
          para_Ni = new double[Ni]();
          memcpy(para_Ni,other.para_Ni,sizeof(double)*1*Ni);
        }
        else
        {
          warning("time-independent interactions");
        }

        if(Ni_td!=0)
        {
          idx_td_Ni = new unsigned int[2*Ni_td]();
          memcpy(idx_td_Ni, other.idx_td_Ni, sizeof(unsigned int)*2*Ni_td);
          idx_td_Ni_ctl0 = new unsigned int[Ni_td]();
          memcpy(idx_td_Ni_ctl0,other.idx_td_Ni_ctl0,sizeof(unsigned int)*1*Ni_td);
          para_td_Ni = new double[Ni_td]();
          memcpy(para_td_Ni,other.para_td_Ni,sizeof(double)*Ni_td);
        }
        else
        {
          warning("time-dependent interactions");
        }
      }
    }

    DEVICE_DATA &operator=(const DEVICE_DATA & other)
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not copied/assignment" << endl;}};
      if(visible==true){cerr << "DD copy assignment constructor" << endl;}

      Nt=other.Nt;
      Nk=other.Nk;
      Ni=other.Ni;
      num_qubits=other.num_qubits;
      Ni_td=other.Ni_td;
      N_ctl0=other.N_ctl0;
      N_ctl1=other.N_ctl1;
      T_start=other.T_start;
      T_end=other.T_end;
      tau=other.tau;
      state_size=other.state_size;
      num_qubits=other.num_qubits;

      if(this!=&other)
      {
        const string option0="free";
        if(option0=="free")
        {

          if(N_ctl0!=0)
          {

            vector<unsigned int>().swap(idx_ctl0);
            vector<vector<double>>().swap(para_ctl0);
            vector<vector<double>>().swap(phase_CZ);
          }
          else
          {
            warning("charge control pulses");
          }


          if( N_ctl1!=0)
          {

            vector<unsigned int>().swap(idx_ctl1);
            vector<vector<double>>().swap(para_ctl1);
          }
          else
          {
            warning("flux control pulses");
          }

          if(num_qubits!=0)
          {
            delete [] idx_state;
            idx_state=nullptr;
            vector<double>().swap(lamda_trans_shifted);
          }
          else
          {
            warning("computatoinal states indices");
          }

          if (Nt!=0)
          {

            delete [] EC;
            EC=nullptr;
            delete [] EJ;
            EJ=nullptr;
            delete [] freq_at_flux_zero;
            freq_at_flux_zero=nullptr;
            delete [] anharm_at_flux_zero;
            anharm_at_flux_zero=nullptr;
            delete [] asymmetry_factor_qubit;
            anharm_at_flux_zero=nullptr;
            delete [] flux_offset;
            flux_offset=nullptr;
            delete [] amp_factor;
            amp_factor=nullptr;
          }
          else
          {
            warning("anharmonic osc. device parameters");
          }

          if(Nk!=0)
          {

            delete [] res_freq;
            res_freq=nullptr;
          }
          else
          {
            warning("harmonic osc. device parameters");
          }

          if(Ni !=0)
          {
            delete [] idx_Ni;
            idx_Ni=nullptr;
            delete [] para_Ni;
            para_Ni=nullptr;
          }
          else
          {
            warning("time-independent interactions");
          }

          if(Ni_td!=0)
          {
            delete [] idx_td_Ni;
            idx_td_Ni=nullptr;
            delete [] idx_td_Ni_ctl0;
            idx_td_Ni_ctl0=nullptr;
            delete [] para_td_Ni;
            para_td_Ni=nullptr;
          }
          else
          {
            warning("time-dependent interactions");
          }
        }

        const string option1="copy";
        if(option1=="copy")
        {
          if(N_ctl0!=0)
          {

            idx_ctl0.resize(N_ctl0,0);
            copy(other.idx_ctl0.begin(), other.idx_ctl0.end(),back_inserter(idx_ctl0));
            para_ctl0.resize(N_ctl0,vector<double>(7,1));
            copy(other.para_ctl0.begin(), other.para_ctl0.end(),back_inserter(para_ctl0));
            phase_CZ.resize(Nt,vector<double>(num_qubits,0));
            copy(other.phase_CZ.begin(), other.phase_CZ.end(),back_inserter(phase_CZ));
          }
          else
          {
            warning("charge control pulses");
          }

          if(N_ctl1!=0)
          {

            idx_ctl1.resize(N_ctl1,0);
            copy(other.idx_ctl1.begin(), other.idx_ctl1.end(),back_inserter(idx_ctl1));
            para_ctl1.resize(N_ctl1,vector<double>(7,1));
            copy(other.para_ctl1.begin(), other.para_ctl1.end(),back_inserter(para_ctl1));
          }
          else
          {
            warning("flux control pulses");
          }

          if(num_qubits!=0)
          {
            idx_state= new u_int64_t[(1 << num_qubits)]();
            memcpy(idx_state,other.idx_state,sizeof(u_int64_t)*(1 << num_qubits));
            lamda_trans_shifted.resize(num_qubits,0);
            copy(other.lamda_trans_shifted.begin(), other.lamda_trans_shifted.end(),back_inserter(lamda_trans_shifted));
          }
          else
          {
            warning("computatoinal states indices");
          }

          if(Nt!=0)
          {

            EC = new double[Nt]();
            memcpy(EC,other.EC,sizeof(double)*1*Nt);
            EJ = new double[2*Nt]();
            memcpy(EJ,other.EJ,sizeof(double)*2*Nt);
            freq_at_flux_zero = new double[Nt]();
            memcpy(freq_at_flux_zero,other.freq_at_flux_zero,sizeof(double)*1*Nt);
            anharm_at_flux_zero = new double[Nt]();
            memcpy(anharm_at_flux_zero,other.anharm_at_flux_zero,sizeof(double)*1*Nt);
            asymmetry_factor_qubit = new double[Nt]();
            memcpy(asymmetry_factor_qubit,other.asymmetry_factor_qubit,sizeof(double)*1*Nt);
            flux_offset = new double[Nt]();
            memcpy(flux_offset,other.flux_offset,sizeof(double)*1*Nt);
            amp_factor= new double[2*Nt]();
            memcpy(amp_factor,other.amp_factor,sizeof(double)*2*Nt);
          }
          else
          {
            warning("anharmonic osc. device parameters");
          }

          if(Nk!=0)
          {
            res_freq=new double[Nk]();
            memcpy(res_freq,other.res_freq,sizeof(double)*1*Nk);
          }
          else
          {
            warning("harmonic osc. device parameters");
          }


          if(Ni!=0)
          {
            idx_Ni = new unsigned int[2*Ni]();
            memcpy(idx_Ni,other.idx_Ni,sizeof(unsigned int)*2*Ni);
            para_Ni = new double[Ni]();
            memcpy(para_Ni,other.para_Ni,sizeof(double)*1*Ni);
          }
          else
          {
            warning("time-independent interactions");
          }

          if(Ni_td!=0)
          {
            idx_td_Ni = new unsigned int[2*Ni_td]();
            memcpy(idx_td_Ni, other.idx_td_Ni, sizeof(unsigned int)*2*Ni_td);
            idx_td_Ni_ctl0 = new unsigned int[Ni_td]();
            memcpy(idx_td_Ni_ctl0,other.idx_td_Ni_ctl0,sizeof(unsigned int)*1*Ni_td);
            para_td_Ni = new double[Ni_td]();
            memcpy(para_td_Ni,other.para_td_Ni,sizeof(double)*Ni_td);
          }
          else
          {
            warning("time-dependent interactions");
          }
        }
      }
      return *this;
    }

    DEVICE_DATA(DEVICE_DATA && other)
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not moved" << endl;}};
      if(visible==true){cerr << "DD move constructor" << endl;}

      Nt=other.Nt;
      Nk=other.Nk;
      Ni=other.Ni;
      num_qubits=other.num_qubits;
      Ni_td=other.Ni_td;
      N_ctl0=other.N_ctl0;
      N_ctl1=other.N_ctl1;
      T_start=other.T_start;
      T_end=other.T_end;
      tau=other.tau;
      state_size=other.state_size;
      num_qubits=other.num_qubits;

      const string option="move";
      if(option=="move")
      {
        if(N_ctl0!=0)
        {
          idx_ctl0=move(other.idx_ctl0);
          para_ctl0=move(other.para_ctl0);
          phase_CZ=move(other.phase_CZ);
        }
        else
        {
          warning("charge control pulses");
        }

        if(N_ctl1!=0)
        {

          idx_ctl1=move(other.idx_ctl1);
          para_ctl1=move(other.para_ctl1);
        }
        else
        {
          warning("flux control pulses");
        }

        if(num_qubits!=0)
        {
          idx_state=other.idx_state;
          other.idx_state=nullptr;
          lamda_trans_shifted=move(other.lamda_trans_shifted);
        }
        else
        {
          warning("computatoinal states indices");
        }

        if(Nt!=0)
        {

          EC = other.EC;
          other.EC=nullptr;
          EJ = other.EJ;
          other.EJ=nullptr;
          freq_at_flux_zero = other.freq_at_flux_zero;
          other.freq_at_flux_zero=nullptr;
          anharm_at_flux_zero = other.anharm_at_flux_zero;
          anharm_at_flux_zero=nullptr;
          asymmetry_factor_qubit = other.asymmetry_factor_qubit;
          other.asymmetry_factor_qubit=nullptr;
          flux_offset = other.flux_offset;
          other.flux_offset=nullptr;
          amp_factor= other.amp_factor;
          other.amp_factor=nullptr;
        }
        else
        {
          warning("anharmonic osc. device parameters");
        }

        if(Nk!=0)
        {
          res_freq=other.res_freq;
          other.res_freq=nullptr;
        }
        else
        {
          warning("harmonic osc. device parameters");
        }


        if(Ni!=0)
        {
          idx_Ni = other.idx_Ni;
          other.idx_Ni=nullptr;
          para_Ni = other.para_Ni;
          other.para_Ni=nullptr;
        }
        else
        {
          warning("time-independent interactions");
        }

        if(Ni_td!=0)
        {
          idx_td_Ni = other.idx_td_Ni;
          other.idx_td_Ni=nullptr;
          idx_td_Ni_ctl0 = other.idx_td_Ni_ctl0;
          other.idx_td_Ni_ctl0=nullptr;
          para_td_Ni = other.para_td_Ni;
          other.para_td_Ni=nullptr;
        }
        else
        {
          warning("time-dependent interactions");
        }
      }

      other.Nt=0;
      other.Nk=0;
      other.Ni=0;
      other.num_qubits=0;
      other.Ni_td=0;
      other.N_ctl0=0;
      other.N_ctl1=0;
      other.T_start=0;
      other.T_end=0;
      other.tau=0;
      other.state_size=0;
      other.num_qubits=0;


    }

    DEVICE_DATA &operator=(DEVICE_DATA && other)
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not moved/assignment" << endl;}};
      if(visible==true){cerr << "DD move assignment constructor" << endl;}



      if(this!=&other)
      {
        const string option0="free";
        if(option0=="free")
        {

          if(N_ctl0!=0)
          {

            vector<unsigned int>().swap(idx_ctl0);
            vector<vector<double>>().swap(para_ctl0);
            vector<vector<double>>().swap(phase_CZ);
          }
          else
          {
            warning("charge control pulses");
          }


          if( N_ctl1!=0)
          {

            vector<unsigned int>().swap(idx_ctl1);
            vector<vector<double>>().swap(para_ctl1);
          }
          else
          {
            warning("flux control pulses");
          }

          if(num_qubits!=0)
          {
            delete [] idx_state;
            idx_state=nullptr;
            vector<double>().swap(lamda_trans_shifted);
          }
          else
          {
            warning("computatoinal states indices");
          }

          if (Nt!=0)
          {

            delete [] EC;
            EC=nullptr;
            delete [] EJ;
            EJ=nullptr;
            delete [] freq_at_flux_zero;
            freq_at_flux_zero=nullptr;
            delete [] anharm_at_flux_zero;
            anharm_at_flux_zero=nullptr;
            delete [] asymmetry_factor_qubit;
            anharm_at_flux_zero=nullptr;
            delete [] flux_offset;
            flux_offset=nullptr;
            delete [] amp_factor;
            amp_factor=nullptr;
          }
          else
          {
            warning("anharmonic osc. device parameters");
          }

          if(Nk!=0)
          {

            delete [] res_freq;
            res_freq=nullptr;
          }
          else
          {
            warning("harmonic osc. device parameters");
          }

          if(Ni !=0)
          {
            delete [] idx_Ni;
            idx_Ni=nullptr;
            delete [] para_Ni;
            para_Ni=nullptr;
          }
          else
          {
            warning("time-independent interactions");
          }

          if(Ni_td!=0)
          {
            delete [] idx_td_Ni;
            idx_td_Ni=nullptr;
            delete [] idx_td_Ni_ctl0;
            idx_td_Ni_ctl0=nullptr;
            delete [] para_td_Ni;
            para_td_Ni=nullptr;
          }
          else
          {
            warning("time-dependent interactions");
          }
        }

        Nt=other.Nt;
        Nk=other.Nk;
        Ni=other.Ni;
        num_qubits=other.num_qubits;
        Ni_td=other.Ni_td;
        N_ctl0=other.N_ctl0;
        N_ctl1=other.N_ctl1;
        T_start=other.T_start;
        T_end=other.T_end;
        tau=other.tau;
        state_size=other.state_size;
        num_qubits=other.num_qubits;

        const string option1="move";
        if(option1=="move")
        {
          if(N_ctl0!=0)
          {
            idx_ctl0=move(other.idx_ctl0);
            para_ctl0=move(other.para_ctl0);
            phase_CZ=move(other.phase_CZ);
          }
          else
          {
            warning("charge control pulses");
          }

          if(N_ctl1!=0)
          {

            idx_ctl1=move(other.idx_ctl1);
            para_ctl1=move(other.para_ctl1);
          }
          else
          {
            warning("flux control pulses");
          }

          if(num_qubits!=0)
          {
            idx_state=other.idx_state;
            other.idx_state=nullptr;
            lamda_trans_shifted=move(other.lamda_trans_shifted);
          }
          else
          {
            warning("computatoinal states indices");
          }

          if(Nt!=0)
          {

            EC = other.EC;
            other.EC=nullptr;
            EJ = other.EJ;
            other.EJ=nullptr;
            freq_at_flux_zero = other.freq_at_flux_zero;
            other.freq_at_flux_zero=nullptr;
            anharm_at_flux_zero = other.anharm_at_flux_zero;
            other.anharm_at_flux_zero=nullptr;
            asymmetry_factor_qubit = other.asymmetry_factor_qubit;
            other.asymmetry_factor_qubit=nullptr;
            flux_offset = other.flux_offset;
            other.flux_offset=nullptr;
            amp_factor= other.amp_factor;
            other.amp_factor=nullptr;
          }
          else
          {
            warning("anharmonic osc. device parameters");
          }

          if(Nk!=0)
          {
            res_freq=other.res_freq;
            other.res_freq=nullptr;
          }
          else
          {
            warning("harmonic osc. device parameters");
          }


          if(Ni!=0)
          {
            idx_Ni = other.idx_Ni;
            other.idx_Ni=nullptr;
            para_Ni = other.para_Ni;
            other.para_Ni=nullptr;
          }
          else
          {
            warning("time-independent interactions");
          }

          if(Ni_td!=0)
          {
            idx_td_Ni = other.idx_td_Ni;
            other.idx_td_Ni=nullptr;
            idx_td_Ni_ctl0 = other.idx_td_Ni_ctl0;
            other.idx_td_Ni_ctl0=nullptr;
            para_td_Ni = other.para_td_Ni;
            other.para_td_Ni=nullptr;
          }
          else
          {
            warning("time-dependent interactions");
          }
        }

        other.Nt=0;
        other.Nk=0;
        other.Ni=0;
        other.num_qubits=0;
        other.Ni_td=0;
        other.N_ctl0=0;
        other.N_ctl1=0;
        other.T_start=0;
        other.T_end=0;
        other.tau=0;
        other.state_size=0;
        other.num_qubits=0;
      }
      return *this;
    }
};

#endif

#ifndef PFA_CPU
#define PFA_CPU

class JUSQUACE_PFA: public DEVICE_DATA
{
  private:
    // The following dynamic array (heap memory) contains the eigenvalues of the
    // anharmonic and harmonic osc Hamiltonians in the model Hamiltonian.s
    double *eigenvalues_device=nullptr;
    double *eigenvalues_device_exp=nullptr;
    // The following dynamic arrays (heap memory) contain the eigenvalues and matrices
    // which diagonalise the matrices (b^{\dagger}+b), (b^{\dagger}-b) and (b^{\dagger}b^{\dagger}-bb)
    // of the charge and flux driving term
    double *eigenvalues_lin_trafo_mat=nullptr;
    double *lin_trafo_mat=nullptr;
    double *lin_trafo_mat_adj=nullptr;

    // The following dynamic arrays (heap memory) contain the exponentials exp(-i const. eigenvalues[i])
    // computed with eigenvalues_device and eigenvalues_lin_trafo_mat.
    double *eigenvalues_inter_exp=nullptr;
    double *eigenvalues_inter_td_exp=nullptr;
    double *eigenvalues_ctl_exp=nullptr;//

    inline void diagonal_phase_transformation(double *state,const vector<double> &vec,const double Time)
    {
      for(unsigned int i=0;i<num_qubits;i++)
      {
        unsigned int size=state_size;
        double hold_real=0,hold_img=0;
        unsigned int hold_ui=0;
        vector<double> eigenvalues_exp(2*Nts,0);
        vector<u_int64_t> fac(Nt+1,0);
        vector<u_int64_t> id(Nt+1,0);

        for(unsigned int j=0;j<(Nt+1);j++)
        {
          fac[j]=pow(Nts,j);
        }

        for(unsigned int j=0;j<Nts;j++)
        {
          eigenvalues_exp[2*j+0]  =cos(j*Time*vec[i]);
          eigenvalues_exp[2*j+1]=sin(j*Time*vec[i]);
        }

        for(unsigned int j=0;j<size;j++)
        {
          hold_ui=j;
          for(unsigned int k=0;k<Nt;k++)
          {
            id[k]= hold_ui % fac[k+1];
            id[k]=id[k]/fac[k];
            hold_ui=hold_ui - id[k]*fac[k];
          }
          hold_real =state[2*j+0];
          hold_img  =state[2*j+1];
          state[2*j+0]  =hold_real*eigenvalues_exp[2*id[i]+0]   -eigenvalues_exp[2*id[i]+1]*hold_img;
          state[2*j+1]  =hold_real*eigenvalues_exp[2*id[i]+1]   +eigenvalues_exp[2*id[i]+0]*hold_img;
        }
      }
    }

    inline void expoential_i_of_eigenvalues_n_times_eigenvalues_m(const double *evals_n,const double *evals_m,double *evals_exp,const double lamda,const unsigned int N,const unsigned int M)
    {

      double hold=0;
      unsigned int idx=0;
      for(unsigned int n=0;n<N;n++)
      {
        for(unsigned int m=0;m<M;m++)
        {
          hold=-lamda*evals_n[n]*evals_m[m];
          evals_exp[2*idx+0]=cos(hold);
          evals_exp[2*idx+1]=sin(hold);
          idx++;
        }
      }
    }

    inline void update_time_dependent_data(const double dt)
    {
      double freq=0;
      double anharm=0;
      double hold=0;
      double tri_factor=0;
      unsigned int idx=0;
      unsigned int offset=0;
      //update all flux-time-dependent data structures
      for(unsigned int i=0;i<N_ctl0;i++)
      {
        idx=idx_ctl0[i];
        const double flux=flux_ctl(para_ctl0[i],dt);
        const double flux_derivative=flux_ctl_der(para_ctl0[i],dt);

        //cerr << " Time/AMP/FLUX/FLUX DER = " << dt << "/" << para_ctl0[i][2]/pi << "/" << flux  << "/" <<  flux_derivative << endl;

        // Compute the osc. freq and anharmonicity with a low (bool_eigenvalues_with_high_accuracy==0) or high (bool_eigenvalues_with_high_accuracy==1) accuracy approximation.
        if(bool_eigenvalues_with_high_accuracy==0)
        {
          //update_eigenvalues(flux,idx);
          freq=flux_freq(freq_at_flux_zero[idx],flux,flux_offset[idx],asymmetry_factor_qubit[idx]);
          anharm=anharm_at_flux_zero[idx];
        }
        else
        {
          //update_eigenvalues_with_high_accuracy(flux,idx);
          freq=flux_freq_with_high_accuracy(EC[idx],EJ[2*idx+0],EJ[2*idx+1],flux,flux_offset[idx]);
          anharm=flux_freq_anharmonicity_with_high_accuracy(EC[idx],EJ[2*idx+0],EJ[2*idx+1],flux,flux_offset[idx]);
        }
        // Update the data structure eigenvalues_device for all anharmonic osc.
        for(unsigned int j=0;j<Nts;j++){eigenvalues_device[j+idx*Nts]=freq*j+(anharm*0.5)*j*(j-1);}
        // Update the data structure exp(-i tau eigenvalues_device) for all anharmonic osc.
        hold=-tau*0.5;
        for(unsigned int j=0;j<Nts;j++)
        {

          eigenvalues_device_exp[2*j+0+idx*2*Nts]=cos(eigenvalues_device[j+idx*Nts]*hold);
          eigenvalues_device_exp[2*j+1+idx*2*Nts]=sin(eigenvalues_device[j+idx*Nts]*hold);
        }
        // Update the data structure exp(-i tau eigenvalues times eigenvalues) for the eigenvalues of the interaction terms
        for(unsigned int j=0;j<Ni_td;j++)
        {
          if(idx_ctl0[i]==idx_td_Ni_ctl0[j])
          {
            hold=time_dependent_interaction_strength(para_td_Ni[j],flux+flux_offset[idx_td_Ni_ctl0[j]],idx_td_Ni_ctl0[j]);
            // expoential_i_of_eigenvalues_n_times_eigenvalues_m(&eigenvalues_lin_trafo_mat[0*Nts],&eigenvalues_lin_trafo_mat[0*Nts],&eigenvalues_inter_td_exp[j*2*Nks*Nts], tau*hold,Nts,Nts);
            unsigned int counter=0;
            offset=j*2*Nks*Nts;
            for(unsigned int k=0;k<Nks;k++)
            {
              for(unsigned int l=0;l<Nts;l++)
              {
                double arg=-hold*tau*eigenvalues_lin_trafo_mat[k]*eigenvalues_lin_trafo_mat[l];
                eigenvalues_inter_td_exp[2*counter+0+offset]=cos(arg);
                eigenvalues_inter_td_exp[2*counter+1+offset]=sin(arg);
                counter++;
              }
            }
          }
        }

        tri_factor=0;
        hold=cos(flux+flux_offset[idx]);//in units of 1/pi not 1/2pi
        tri_factor+=hold*hold;
        hold=sin(flux+flux_offset[idx]);//in units of 1/pi not 1/2pi
        tri_factor+=(asymmetry_factor_qubit[idx]*asymmetry_factor_qubit[idx]*hold*hold);

        // Update the data structure exp(-i tau eigenvalues ) for the first part (b^{\dagger} - b) of the flux driving terms
        hold=amp_factor[2*idx+0]*pow(tri_factor, -0.875)*2*flux_derivative;
        hold=-tau*hold*0.5;
        offset=i*2*Nts+max_num_charge_drive*2*Nts;
        for(unsigned int j=0;j<Nts;j++)
        {
          eigenvalues_ctl_exp[2*j+0+offset]=cos(eigenvalues_lin_trafo_mat[j+1*Nts]*hold);
          eigenvalues_ctl_exp[2*j+1+offset]=sin(eigenvalues_lin_trafo_mat[j+1*Nts]*hold);
        }

        // Update the data structure exp(-i tau eigenvalues ) for the second part (b^{\dagger}b^{\dagger} - bb) of the flux driving terms
        hold=(amp_factor[2*idx+1]*sin(2*(flux+flux_offset[idx]))*2*flux_derivative)/tri_factor;//in units of 1/pi not 1/2pi
        hold=-tau*hold*0.5;
        offset=i*2*Nts+max_num_flux_drive*2*Nts+max_num_charge_drive*2*Nts;
        for(unsigned int j=0;j<Nts;j++)
        {
          eigenvalues_ctl_exp[2*j+0+offset]=cos(eigenvalues_lin_trafo_mat[j+2*Nts]*hold);
          eigenvalues_ctl_exp[2*j+1+offset]=sin(eigenvalues_lin_trafo_mat[j+2*Nts]*hold);
        }

      }

      // Update the data structure exp(-i tau eigenvalues ) for the charge driving terms (b^{\dagger} + b).
      for(unsigned int i=0;i<N_ctl1;i++)
      {
        hold=charge_ctl(para_ctl1[i],dt);
        //expoential_i_of_eigenvalues(&eigenvalues_lin_trafo_mat[0*Nts], &eigenvalues_ctl_exp[i*2*Nts+0+0],-tau*hold,Nts);
        hold=-tau*hold*1.0;
        offset=0;
        for(unsigned int j=0;j<Nts;j++)
        {
          eigenvalues_ctl_exp[2*j+0+offset]=cos(eigenvalues_lin_trafo_mat[j+0*Nts]*hold);
          eigenvalues_ctl_exp[2*j+1+offset]=sin(eigenvalues_lin_trafo_mat[j+0*Nts]*hold);
        }
      }
    }

    inline void initial_pfa_data_structures(void)
    {

      if(Nts!=0 || Nks!=0)
      {
        // Define offset value for matrix data structure (float row major format)
        const unsigned int offset=2*Nts;
        const unsigned int offset_Nts=2*Nts*Nts;
        // Lamda function defined to compute the adjoint of a complex matrix in row major format.
        auto adjoint = [](const double *matrix,double *matrix_adjoint,const unsigned int N,const unsigned int offset)
        {
            unsigned int idx,idy;
            for(unsigned int n=0;n<N;n++)
            {
              for(unsigned int m=0;m<N;m++)
              {
                idx=2*m+offset*n;
                idy=2*n+offset*m;
                matrix_adjoint[idx]=matrix[idy];
                matrix_adjoint[idx+1]=-matrix[idy+1];
              }
            }
          };


        // Rest all matrix and eigenvalue data structures
        memset(lin_trafo_mat, 0, 3*offset_Nts*sizeof(double));
        memset(lin_trafo_mat_adj, 0, 3*offset_Nts*sizeof(double));
        memset(eigenvalues_lin_trafo_mat, 0, 3*Nts*sizeof(double));


        // Initialisation and diagonalisation of the matrix (a^{\dagger}+a) for Nts states
        for(unsigned int i=0;i<Nts;i++)
        {
          if(i<(Nts-1))
          {
            lin_trafo_mat[2*(i+1)+0+offset*i+0*offset_Nts]+=sqrt(i+1);
          }
          if(i>0)
          {
            lin_trafo_mat[2*(i-1)+0+offset*i+0*offset_Nts]+=sqrt(i-0);
          }
        }

        LAPACKE_zheev( LAPACK_ROW_MAJOR, 'V', 'U', Nts, &lin_trafo_mat[0*offset_Nts], Nts, &eigenvalues_lin_trafo_mat[0*Nts]);
        adjoint(&lin_trafo_mat[0*offset_Nts],&lin_trafo_mat_adj[0*offset_Nts],Nts,2*Nts);

        // Initialisation and diagonalisation of the matrix (a^{\dagger}-a) for Nts states
        for(unsigned int i=0;i<Nts;i++)
        {
          if(i<(Nts-1)){lin_trafo_mat[2*(i+1)+1+offset*i+1*offset_Nts]+=-sqrt(i+1);}
          if(i>0){lin_trafo_mat[2*(i-1)+1+offset*i+1*offset_Nts]+=sqrt(i-0);}
        }

        LAPACKE_zheev( LAPACK_ROW_MAJOR, 'V', 'U', Nts, &lin_trafo_mat[1*offset_Nts], Nts, &eigenvalues_lin_trafo_mat[1*Nts]);
        adjoint(&lin_trafo_mat[1*offset_Nts],&lin_trafo_mat_adj[1*offset_Nts],Nts,2*Nts);

        // Initialisation and diagonalisation of the matrix (a^{\dagger}a^{\dagger}-aa) for Nts states
        for(unsigned int i=0;i<Nts;i++)
        {
          if(i<(Nts-2)){lin_trafo_mat[2*(i+2)+1+offset*i+2*offset_Nts]+=-sqrt((i+2)*(i+1));}
          if(i>1){lin_trafo_mat[2*(i-2)+1+offset*i+2*offset_Nts]+=sqrt((i-1)*(i-0));}
        }

        LAPACKE_zheev( LAPACK_ROW_MAJOR, 'V', 'U', Nts, &lin_trafo_mat[2*offset_Nts], Nts, &eigenvalues_lin_trafo_mat[2*Nts]);
        adjoint(&lin_trafo_mat[2*offset_Nts],&lin_trafo_mat_adj[2*offset_Nts],Nts,2*Nts);
      }

      if(num_qubits!=0)
      {
        u_int64_t dim=1 << num_qubits;
        for(uint64_t i=0;i<dim;i++)
        {
          idx_state[i]=0;
          for(unsigned int j=0;j<Nt;j++)
          {
            uint64_t mask = 1 << j;
            uint64_t hold = (i & mask) >> j;
            idx_state[i]+=pow(Nts,j)*hold;
          }
        }
      }

      if(Nt!=0)
      {
        for(unsigned int i=0;i<Nt;i++)
        {
          for(unsigned int j=0;j<Nts;j++)
          {
            eigenvalues_device[j+Nts*i]=freq_at_flux_zero[i]*j+(anharm_at_flux_zero[i]*0.5)*j*(j-1);
          }
          double hold=-tau*0.5;
          unsigned int offset_left=i*2*Nts;
          unsigned int offset_right=i*Nts;
          for(unsigned int j=0;j<Nts;j++)
          {
            eigenvalues_device_exp[2*j+0+offset_left]=cos(eigenvalues_device[j+offset_right]*hold);
            eigenvalues_device_exp[2*j+1+offset_left]=sin(eigenvalues_device[j+offset_right]*hold);
          }
        }
      }

      if(Nk!=0)
      {
        for(unsigned int i=0;i<Nk;i++)
        {

          for(unsigned int j=0;j<Nks;j++)
          {
            eigenvalues_device[j+i*Nks+Nt*Nts]=res_freq[i]*j;
          }
          double hold=-tau*0.5;
          unsigned int offset_left=i*2*Nks+Nt*2*Nts;
          unsigned int offset_right=i*Nks+Nt*Nts;
          for(unsigned int j=0;j<Nks;j++)
          {
            eigenvalues_device_exp[2*j+0+offset_left]=cos(eigenvalues_device[j+offset_right]*hold);
            eigenvalues_device_exp[2*j+1+offset_left]=sin(eigenvalues_device[j+offset_right]*hold);
          }
        }
      }

      if(Ni!=0)
      {
        for(unsigned int i=0;i<Ni;i++)
        {
          if(idx_Ni[2*i+1]>idx_Ni[2*i])
          {

            unsigned int counter=0;
            unsigned int offset=i*2*Nks*Nts;
            for(unsigned int j=0;j<Nks;j++)
            {
              for(unsigned int k=0;k<Nts;k++)
              {
                double arg=-para_Ni[i]*tau*eigenvalues_lin_trafo_mat[j]*eigenvalues_lin_trafo_mat[k];
                eigenvalues_inter_exp[2*counter+0+offset]=cos(arg);
                eigenvalues_inter_exp[2*counter+1+offset]=sin(arg);
                counter++;
              }
            }
          }
          else
          {
            cout << "#index error interaction"<<endl;
          }
        }
      }

      if(Ni_td!=0)
      {
        for(unsigned int i=0;i<Ni_td;i++)
        {
          if(idx_td_Ni[2*i+1]>idx_td_Ni[2*i])
          {
            //expoential_i_of_eigenvalues_n_times_eigenvalues_m(&eigenvalues_lin_trafo_mat[0*Nts],&eigenvalues_lin_trafo_mat[0*Nts],&eigenvalues_inter_td_exp[i*2*Nks*Nts], tau*para_td_Ni[i],4,4);
          }
          else
          {
            cout << "#index error interaction"<<endl;
          }
        }
      }

    }

    inline void diagonal_operator_2_bit(double *state,const double *transformation,const uint64_t target_idx,const uint64_t bit_offset,const uint64_t size);

    inline void diagonal_operator_2_bit_X_2_bit(double *state,const double *transformation,const uint64_t target_idx_1,const uint64_t bit_offset_1,const uint64_t target_idx_0,const uint64_t bit_offset_0,const uint64_t size);

    inline void basis_transformation_2_bit_real(double *state,const double *transformation,const uint64_t target_idx,const uint64_t bit_offset,const uint64_t size);

    inline void basis_transformation_2_bit_complex(double *state,const double *transformation,const uint64_t target_idx,const uint64_t bit_offset,const uint64_t size);

    inline void sample_time_evolution(const double *state,const double dt,const vector<uint64_t> tar_idx);
  public:
    inline void pfa_time_evolution(double *state, const vector<uint64_t> tar_idx=vector<uint64_t>());

    JUSQUACE_PFA() = default ;

    JUSQUACE_PFA(vector<double> qubit_freq, vector<double> anharm,vector<double> flux_offset_in,vector<double> asymmetry_factor,vector<double>resonator_freq ,unsigned int num_transmons,unsigned int num_resonators,
      unsigned int num_flux_drives,unsigned int num_charge_drives,
      unsigned int num_interactions,vector<unsigned int> interaction_idx,vector<double> interaction_para,
      double time_grid_parameter,unsigned int num_iter_steps,
      unsigned int num_td_interactions=0,vector<unsigned int> interaction_td_idx={},vector<unsigned int> interaction_td_idx_ctl0={},vector<double> interaction_td_para={},
      vector<double> EC_input={},vector<double> EJ_input={},unsigned int num_qubit=0):DEVICE_DATA(qubit_freq,anharm,flux_offset_in,asymmetry_factor,resonator_freq,num_transmons,num_resonators,num_flux_drives,num_charge_drives,num_interactions,interaction_idx,interaction_para,time_grid_parameter,num_iter_steps,num_td_interactions,interaction_td_idx,interaction_td_idx_ctl0,interaction_td_para,EC_input,EJ_input,num_qubit)
    {


      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not allocated and/or freed" << endl;}};
      if(visible==true){cerr << "PFA constructor" << endl;}

      const string option="allocate";
      if(option=="allocate")
      {
        if((Nt!=0 || Nk!=0))
        {
          eigenvalues_device = new double[Nt*Nts+Nk*Nks]();
          eigenvalues_device_exp = new double[Nt*2*Nts+Nk*2*Nks]();
        }
        else
        {
          warning("device data structures");
        }

        if((Nts+Nks)!=0)
        {

          eigenvalues_lin_trafo_mat= new double[3*Nts]();
          lin_trafo_mat=new double[3*2*Nts*Nts];
          lin_trafo_mat_adj=new double[3*2*Nts*Nts];
        }
        else
        {
          warning("pfa data structures");
        }

        if(Ni!=0)
        {

          eigenvalues_inter_exp = new double[Ni*2*Nks*Nts]();
        }
        else
        {

          warning("time-independent interactions");
        }

        if(Ni_td!=0)
        {

          eigenvalues_inter_td_exp = new double[Ni_td*2*Nks*Nts]();
        }
        else
        {

          warning("time-dependent interactions");
        }

        if((max_num_flux_drive !=0 || max_num_charge_drive!=0))
        {

          eigenvalues_ctl_exp = new double[max_num_charge_drive*2*Nts + max_num_flux_drive*2*Nts + max_num_flux_drive*2*Nts]();
        }
        else
        {
          warning("diagonal time-dependent data structures");
        }
      }
      initial_pfa_data_structures();
    }

    ~JUSQUACE_PFA()
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not freed" << endl;}};
      if(visible==true){cerr << "PFA destructor" << endl;}

      const string option="free";
      if(option=="free")
      {
        if((Nt!=0 || Nk!=0))
        {
          delete [] eigenvalues_device;
          eigenvalues_device=nullptr;
          delete [] eigenvalues_device_exp;
          eigenvalues_device_exp=nullptr;
        }
        else
        {
          warning("device data structures");
        }

        if((Nts+Nks)!=0)
        {

          delete [] eigenvalues_lin_trafo_mat;
          eigenvalues_lin_trafo_mat=nullptr;
          delete [] lin_trafo_mat;
          lin_trafo_mat=nullptr;
          delete [] lin_trafo_mat_adj;
          lin_trafo_mat_adj=nullptr;
        }
        else
        {
          warning("pfa data structures");
        }

        if (Ni!=0)
        {

          delete [] eigenvalues_inter_exp;
          eigenvalues_inter_exp=nullptr;
        }
        else
        {

          warning("time-independent interactions");
        }

        if (Ni_td!=0)
        {

          delete [] eigenvalues_inter_td_exp;
          eigenvalues_inter_td_exp=nullptr;
        }
        else
        {

          warning("time-dependent interactions");
        }

        if((max_num_flux_drive !=0 || max_num_charge_drive!=0))
        {

          delete [] eigenvalues_ctl_exp;
          eigenvalues_ctl_exp=nullptr;
        }
        else
        {
          warning("diagonal time-dependent data structures");
        }
      }
    }

    JUSQUACE_PFA(const JUSQUACE_PFA & other):DEVICE_DATA(other)
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not copied" << endl;}};
      if(visible==true){cerr << "PFA copy constructor" << endl;}

      const string option="copy";
      if(option=="copy")
      {
        if((Nt!=0 || Nk!=0))
        {
          eigenvalues_device = new double[Nt*Nts+Nk*Nks]();
          memcpy(eigenvalues_device,other.eigenvalues_device,sizeof(double)*(Nt*Nts+Nk*Nks));
          eigenvalues_device_exp = new double[Nt*2*Nts+Nk*2*Nks]();
          memcpy(eigenvalues_device_exp,other.eigenvalues_device_exp,sizeof(double)*(Nt*2*Nts+Nk*2*Nks));
        }
        else
        {
          warning("device data structures");
        }


        if((Nts+Nks)!=0)
        {

          eigenvalues_lin_trafo_mat= new double[3*Nts]();
          memcpy(eigenvalues_lin_trafo_mat,other.eigenvalues_lin_trafo_mat,sizeof(double)*3*Nts);
          lin_trafo_mat=new double[3*2*Nts*Nts];
          memcpy(lin_trafo_mat,other.lin_trafo_mat,sizeof(double)*(3*2*Nts*Nts));
          lin_trafo_mat_adj=new double[3*2*Nts*Nts];
          memcpy(lin_trafo_mat_adj,other.lin_trafo_mat_adj,sizeof(double)*(3*2*Nts*Nts));
        }
        else
        {
          warning("pfa data structures");
        }

        if(Ni!=0)
        {

          eigenvalues_inter_td_exp = new double[Ni_td*2*Nks*Nts]();
          memcpy(eigenvalues_inter_td_exp,other.eigenvalues_inter_td_exp,sizeof(double)*(Ni_td*2*Nks*Nts));
        }
        else
        {

          warning("time-independent interactions");
        }

        if(Ni_td!=0)
        {

          eigenvalues_inter_td_exp = new double[Ni_td*2*Nks*Nts]();
          memcpy(eigenvalues_inter_td_exp,other.eigenvalues_inter_td_exp,sizeof(double)*(Ni_td*2*Nks*Nts));
        }
        else
        {

          warning("time-dependent interactions");
        }

        if((max_num_flux_drive !=0 || max_num_charge_drive!=0))
        {

          eigenvalues_ctl_exp = new double[max_num_charge_drive*2*Nts + max_num_flux_drive*2*Nts + max_num_flux_drive*2*Nts]();
          memcpy(eigenvalues_ctl_exp,other.eigenvalues_ctl_exp,sizeof(double)*(max_num_charge_drive*2*Nts + max_num_flux_drive*2*Nts + max_num_flux_drive*2*Nts));
        }
        else
        {
          warning("diagonal time-dependent data structures");
        }
      }

    }

    JUSQUACE_PFA &operator=(const JUSQUACE_PFA & other)
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not copied/assignment" << endl;}};
      if(visible==true){cerr << "PFA copy assignment constructor" << endl;}

      if(this!=&other)
      {
        DEVICE_DATA::operator=(other);

        const string option0="free";
        if(option0=="free")
        {
          if((Nt!=0 || Nk!=0))
          {
            delete [] eigenvalues_device;
            eigenvalues_device=nullptr;
            delete [] eigenvalues_device_exp;
            eigenvalues_device_exp=nullptr;
          }
          else
          {
            warning("device data structures");
          }

          if((Nts+Nks)!=0)
          {

            delete [] eigenvalues_lin_trafo_mat;
            eigenvalues_lin_trafo_mat=nullptr;
            delete [] lin_trafo_mat;
            lin_trafo_mat=nullptr;
            delete [] lin_trafo_mat_adj;
            lin_trafo_mat_adj=nullptr;
          }
          else
          {
            warning("pfa data structures");
          }

          if (Ni!=0)
          {

            delete [] eigenvalues_inter_exp;
            eigenvalues_inter_exp=nullptr;
          }
          else
          {

            warning("time-independent interactions");
          }

          if (Ni_td!=0)
          {

            delete [] eigenvalues_inter_td_exp;
            eigenvalues_inter_td_exp=nullptr;
          }
          else
          {

            warning("time-dependent interactions");
          }

          if((max_num_flux_drive !=0 || max_num_charge_drive!=0))
          {

            delete [] eigenvalues_ctl_exp;
            eigenvalues_ctl_exp=nullptr;
          }
          else
          {
            warning("diagonal time-dependent data structures");
          }
        }

        const string option1="copy";
        if(option1=="copy")
        {
          if((Nt!=0 || Nk!=0))
          {
            eigenvalues_device = new double[Nt*Nts+Nk*Nks]();
            memcpy(eigenvalues_device,other.eigenvalues_device,sizeof(double)*(Nt*Nts+Nk*Nks));
            eigenvalues_device_exp = new double[Nt*2*Nts+Nk*2*Nks]();
            memcpy(eigenvalues_device_exp,other.eigenvalues_device_exp,sizeof(double)*(Nt*2*Nts+Nk*2*Nks));
          }
          else
          {
            warning("device data structures");
          }


          if((Nts+Nks)!=0)
          {

            eigenvalues_lin_trafo_mat= new double[3*Nts]();
            memcpy(eigenvalues_lin_trafo_mat,other.eigenvalues_lin_trafo_mat,sizeof(double)*3*Nts);
            lin_trafo_mat=new double[3*2*Nts*Nts];
            memcpy(lin_trafo_mat,other.lin_trafo_mat,sizeof(double)*(3*2*Nts*Nts));
            lin_trafo_mat_adj=new double[3*2*Nts*Nts];
            memcpy(lin_trafo_mat_adj,other.lin_trafo_mat_adj,sizeof(double)*(3*2*Nts*Nts));
          }
          else
          {
            warning("pfa data structures");
          }

          if(Ni!=0)
          {

            eigenvalues_inter_td_exp = new double[Ni_td*2*Nks*Nts]();
            memcpy(eigenvalues_inter_td_exp,other.eigenvalues_inter_td_exp,sizeof(double)*(Ni_td*2*Nks*Nts));
          }
          else
          {

            warning("time-independent interactions");
          }

          if(Ni_td!=0)
          {

            eigenvalues_inter_td_exp = new double[Ni_td*2*Nks*Nts]();
            memcpy(eigenvalues_inter_td_exp,other.eigenvalues_inter_td_exp,sizeof(double)*(Ni_td*2*Nks*Nts));
          }
          else
          {

            warning("time-dependent interactions");
          }

          if((max_num_flux_drive !=0 || max_num_charge_drive!=0))
          {

            eigenvalues_ctl_exp = new double[max_num_charge_drive*2*Nts + max_num_flux_drive*2*Nts + max_num_flux_drive*2*Nts]();
            memcpy(eigenvalues_ctl_exp,other.eigenvalues_ctl_exp,sizeof(double)*(max_num_charge_drive*2*Nts + max_num_flux_drive*2*Nts + max_num_flux_drive*2*Nts));
          }
          else
          {
            warning("diagonal time-dependent data structures");
          }
        }
      }
      return *this;
    }

    JUSQUACE_PFA(JUSQUACE_PFA && other):DEVICE_DATA(move(other))
    {


      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not moved" << endl;}};
      if(visible==true){cerr << "PFA move constructors" << endl;}

      const string option="move";
      if(option=="move")
      {
        if((Nt!=0 || Nk!=0))
        {
          eigenvalues_device = other.eigenvalues_device;
          other.eigenvalues_device=nullptr;
          eigenvalues_device_exp = other.eigenvalues_device_exp;
          other.eigenvalues_device_exp=nullptr;
        }
        else
        {
          warning("device data structures");
        }


        if((Nts+Nks)!=0)
        {

          eigenvalues_lin_trafo_mat= other.eigenvalues_lin_trafo_mat;
          other.eigenvalues_lin_trafo_mat=nullptr;
          lin_trafo_mat=other.lin_trafo_mat;
          other.lin_trafo_mat=nullptr;
          lin_trafo_mat_adj=other.lin_trafo_mat_adj;
          lin_trafo_mat_adj=nullptr;
        }
        else
        {
          warning("pfa data structures");
        }

        if(Ni!=0)
        {

          eigenvalues_inter_td_exp = other.eigenvalues_inter_td_exp;
          other.eigenvalues_inter_td_exp=nullptr;
        }
        else
        {

          warning("time-independent interactions");
        }

        if(Ni_td!=0)
        {

          eigenvalues_inter_td_exp = other.eigenvalues_inter_td_exp;
          other.eigenvalues_inter_td_exp=nullptr;
        }
        else
        {

          warning("time-dependent interactions");
        }

        if((max_num_flux_drive !=0 || max_num_charge_drive!=0))
        {

          eigenvalues_ctl_exp = other.eigenvalues_ctl_exp;
          other.eigenvalues_ctl_exp=nullptr;
        }
        else
        {
          warning("diagonal time-dependent data structures");
        }
      }
    }

    JUSQUACE_PFA &operator=(JUSQUACE_PFA && other)
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not moved/assignment" << endl;}};
      if(visible==true){cerr << "PFA move assignment constructor" << endl;}

      if(this!=&other)
      {
        DEVICE_DATA::operator=(move(other));

        const string option0="free";
        if(option0=="free")
        {
          if((Nt!=0 || Nk!=0))
          {
            delete [] eigenvalues_device;
            eigenvalues_device=nullptr;
            delete [] eigenvalues_device_exp;
            eigenvalues_device_exp=nullptr;
          }
          else
          {
            warning("device data structures");
          }

          if((Nts+Nks)!=0)
          {

            delete [] eigenvalues_lin_trafo_mat;
            eigenvalues_lin_trafo_mat=nullptr;
            delete [] lin_trafo_mat;
            lin_trafo_mat=nullptr;
            delete [] lin_trafo_mat_adj;
            lin_trafo_mat_adj=nullptr;
          }
          else
          {
            warning("pfa data structures");
          }

          if (Ni!=0)
          {

            delete [] eigenvalues_inter_exp;
            eigenvalues_inter_exp=nullptr;
          }
          else
          {

            warning("time-independent interactions");
          }

          if (Ni_td!=0)
          {

            delete [] eigenvalues_inter_td_exp;
            eigenvalues_inter_td_exp=nullptr;
          }
          else
          {

            warning("time-dependent interactions");
          }

          if((max_num_flux_drive !=0 || max_num_charge_drive!=0))
          {

            delete [] eigenvalues_ctl_exp;
            eigenvalues_ctl_exp=nullptr;
          }
          else
          {
            warning("diagonal time-dependent data structures");
          }
        }

        const string option1="move";
        if(option1=="move")
        {
          if((Nt!=0 || Nk!=0))
          {
            eigenvalues_device = other.eigenvalues_device;
            other.eigenvalues_device=nullptr;
            eigenvalues_device_exp = other.eigenvalues_device_exp;
            other.eigenvalues_device_exp=nullptr;
          }
          else
          {
            warning("device data structures");
          }


          if((Nts+Nks)!=0)
          {

            eigenvalues_lin_trafo_mat= other.eigenvalues_lin_trafo_mat;
            other.eigenvalues_lin_trafo_mat=nullptr;
            lin_trafo_mat=other.lin_trafo_mat;
            other.lin_trafo_mat=nullptr;
            lin_trafo_mat_adj=other.lin_trafo_mat_adj;
            other.lin_trafo_mat_adj=nullptr;
          }
          else
          {
            warning("pfa data structures");
          }

          if(Ni!=0)
          {

            eigenvalues_inter_exp = other.eigenvalues_inter_exp;
            other.eigenvalues_inter_exp=nullptr;
          }
          else
          {

            warning("time-independent interactions");
          }

          if(Ni_td!=0)
          {
            eigenvalues_inter_td_exp = other.eigenvalues_inter_td_exp;
            other.eigenvalues_inter_td_exp=nullptr;
          }
          else
          {

            warning("time-dependent interactions");
          }

          if((max_num_flux_drive !=0 || max_num_charge_drive!=0))
          {

            eigenvalues_ctl_exp = other.eigenvalues_ctl_exp;
            other.eigenvalues_ctl_exp=nullptr;
          }
          else
          {
            warning("diagonal time-dependent data structures");
          }
        }
      }
      return *this;
    }
};

inline void JUSQUACE_PFA::diagonal_operator_2_bit(double *state,const double *transformation,const uint64_t target_idx,const uint64_t bit_offset,const uint64_t size)
{
  double hold_real=0;
  double hold_img=0;
  uint64_t mask = (1 << ((2*target_idx+0)+ bit_offset)) + (1 << ((2*target_idx+1)+bit_offset));
  uint64_t idx=0;
  #pragma omp for
  for(unsigned int i=0;i<size;i++)
  {
    idx=(i & mask);
    idx=(i & mask) >> (2*target_idx + bit_offset);
    hold_real =state[2*i+0]*transformation[2*idx+0]- transformation[2*idx+1]*state[2*i+1];
    hold_img  =state[2*i+0]*transformation[2*idx+1]+ transformation[2*idx+0]*state[2*i+1];
    state[2*i+0]=hold_real;
    state[2*i+1]=hold_img;
  }
}

inline void JUSQUACE_PFA::diagonal_operator_2_bit_X_2_bit(double *state,const double *transformation,const uint64_t target_idx_1,const uint64_t bit_offset_1,const uint64_t target_idx_0,const uint64_t bit_offset_0,const uint64_t size)
{
  const uint64_t mask_0 = (1 << ((2*target_idx_0+0) + bit_offset_0)) + (1 << ((2*target_idx_0+1)+bit_offset_0));
  const uint64_t mask_1 = (1 << ((2*target_idx_1+0) + bit_offset_1)) + (1 << ((2*target_idx_1+1)+bit_offset_1));
  double hold_real=0;
  double hold_img=0;
  uint64_t idx_0=0;
  uint64_t idx_1=0;
  #pragma omp for
  for(unsigned int i=0;i<size;i++)
  {
    idx_0=(i & mask_0) >> (2*target_idx_0+bit_offset_0);
    idx_1=(i & mask_1) >> (2*target_idx_1+bit_offset_1);
    idx_0=idx_1*4+idx_0;
    hold_real =state[2*i+0]*transformation[2*idx_0+0]- transformation[2*idx_0+1]*state[2*i+1];
    hold_img  =state[2*i+0]*transformation[2*idx_0+1]+ transformation[2*idx_0+0]*state[2*i+1];
    state[2*i+0]=hold_real;
    state[2*i+1]=hold_img;
  }
}

inline void JUSQUACE_PFA::basis_transformation_2_bit_real(double *state,const double *transformation,const uint64_t target_idx,const uint64_t bit_offset,const uint64_t size)
{

  const uint64_t inc  = 1 << (2*(target_idx+0)+bit_offset);
  const uint64_t incl = 1 << (2*(target_idx+1)+bit_offset);
  uint64_t idx  = 0;
  double hold_state[2*4]={};

  #pragma omp for collapse(2)
  for(uint64_t K = 0; K < size; K += incl)
  {
    for(uint64_t M = 0; M < inc; M++)
    {
      for(unsigned int i=0;i<4;i++)
      {
        for(unsigned int j=0;j<4;j++)
        {
          idx = (K|M)+ j*inc;
          hold_state[2*i+0]+= state[2*idx+0]*transformation[2*j+0+8*i];
          hold_state[2*i+1]+= state[2*idx+1]*transformation[2*j+0+8*i];
        }
      }
      for(unsigned int i=0;i<4;i++)
      {
        idx = (K|M)+ i*inc;
        state[2*idx+0]=hold_state[2*i+0];
        state[2*idx+1]=hold_state[2*i+1];
        hold_state[2*i+0]=0;
        hold_state[2*i+1]=0;
      }
    }
  }
}

inline void JUSQUACE_PFA::basis_transformation_2_bit_complex(double *state,const double *transformation,const uint64_t target_idx,const uint64_t bit_offset,const uint64_t size)
{
  const uint64_t inc  = 1 << (2*(target_idx+0)+bit_offset);
  const uint64_t incl = 1 << (2*(target_idx+1)+bit_offset);
  uint64_t idx  = 0;
  double hold_state[2*4]={};

  #pragma omp for collapse(2)
  for(uint64_t K = 0; K < size; K += incl)
  {
    for(uint64_t M = 0; M < inc; M++)
    {
      for(unsigned int i=0;i<4;i++)
      {
        for(unsigned int j=0;j<4;j++)
        {
          idx = (K|M)+ j*inc;
          hold_state[2*i+0]+= state[2*idx+0]*transformation[2*j+0+8*i] - transformation[2*j+1+8*i]*state[2*idx+1];
          hold_state[2*i+1]+= state[2*idx+0]*transformation[2*j+1+8*i] + transformation[2*j+0+8*i]*state[2*idx+1];
        }
      }
      for(unsigned int i=0;i<4;i++)
      {
        idx = (K|M)+ i*inc;
        state[2*idx+0]=hold_state[2*i+0];
        state[2*idx+1]=hold_state[2*i+1];
        hold_state[2*i+0]=0;
        hold_state[2*i+1]=0;
      }
    }
  }
}

inline void JUSQUACE_PFA::sample_time_evolution(const double *state,const double dt, const vector<uint64_t> tar_idx)
{
  const unsigned int width=30;
  const unsigned int ID=0;
  double hold=0;
  cerr << setw(width) << left << para_ctl0[ID][1]/(2*pi);
  cerr << setw(width) << left <<  dt;
  for(unsigned int i=0;i<tar_idx.size();i++)
  {
    hold=state[2*tar_idx[i]+0]*state[2*tar_idx[i]+0]+state[2*tar_idx[i]+1]*state[2*tar_idx[i]+1];
    cerr << setw(width) << left << hold;
  }
  cerr << "\n";
}



void JUSQUACE_PFA::pfa_time_evolution(double *state, const vector<uint64_t> tar_idx)
{

  double dt=0;
  const unsigned int offset_Nts=2*Nts*Nts;
  const bool sample=(tar_idx.size()!=0);
  //Initial the time-dependent dipole-dipole interaction strengths;
  for(unsigned int i=0;i<Ni_td;i++)
  {

    unsigned int idx=idx_td_Ni_ctl0[i];
    double hold= time_dependent_interaction_strength(para_td_Ni[i],0+flux_offset[idx],idx);
    expoential_i_of_eigenvalues_n_times_eigenvalues_m(&eigenvalues_lin_trafo_mat[0*Nts],&eigenvalues_lin_trafo_mat[0*Nts],&eigenvalues_inter_td_exp[i*2*Nks*Nts], tau*hold,Nts,Nts);
  }

  // Time-evolution loop
  #pragma omp parallel
  for(unsigned int t=T_start;t<T_end;t++)
  {
    //**************** update qubit eigenvalues ****************//
    #pragma omp single
    {
      dt=(t+0.5)*tau;
      update_time_dependent_data(dt);
    }

    //**************** exp(-i tau/2 H0) ****************//
    for(unsigned int i=0;i<Nt;i++)
    {
      diagonal_operator_2_bit(state,&eigenvalues_device_exp[i*2*Nts],i,0,state_size);
    }
    for(unsigned int i=0;i<Nk;i++)
    {

      diagonal_operator_2_bit(state,&eigenvalues_device_exp[i*2*Nks+Nt*2*Nts],i,2*Nt,state_size);
    }

    #if (ADIABATIC==0)
    // //**************** exp(-i tau/2 DrivetermFlux) ****************//
    for(unsigned int i=0;i<N_ctl0;i++)// Vi^d
    {

      basis_transformation_2_bit_complex(state,&lin_trafo_mat_adj[1*offset_Nts],idx_ctl0[i],0,state_size);
      diagonal_operator_2_bit(state,&eigenvalues_ctl_exp[i*2*Nts+max_num_charge_drive*2*Nts],idx_ctl0[i],0,state_size);
      basis_transformation_2_bit_complex(state,&lin_trafo_mat[1*offset_Nts],idx_ctl0[i],0,state_size);
      basis_transformation_2_bit_complex(state,&lin_trafo_mat_adj[2*offset_Nts],idx_ctl0[i],0,state_size);
      diagonal_operator_2_bit(state,&eigenvalues_ctl_exp[i*2*Nts+max_num_flux_drive*2*Nts+max_num_charge_drive*2*Nts],idx_ctl0[i],0,state_size);
      basis_transformation_2_bit_complex(state,&lin_trafo_mat[2*offset_Nts],idx_ctl0[i],0,state_size);
    }
    #endif


    //**************** V^d such that V * V^d =I  ****************//
    for(unsigned int i=0;i<Nt;i++)// Vi^d
    {
      basis_transformation_2_bit_real(state,&lin_trafo_mat_adj[0],i,0,state_size);

    }
    for(unsigned int i=0;i<Nk;i++)// Vi^d
    {
      basis_transformation_2_bit_real(state,&lin_trafo_mat_adj[0],i,2*Nt,state_size);
    }



    //**************** // g_i_j (b+b^d)j (b+b^d)i  ****************//
    for(unsigned int i=0;i<Ni;i++)// g_i_j (ai+ai^d) (aj+aj^d)
    {
      diagonal_operator_2_bit_X_2_bit(state,&eigenvalues_inter_exp[i*2*Nks*Nts],idx_Ni[2*i+1],0,idx_Ni[2*i+0],0,state_size);
    }



    for(unsigned int i=0;i<Ni_td;i++)// g_i_j(t) (ai+ai^d) (aj+aj^d)
    {

      diagonal_operator_2_bit_X_2_bit(state,&eigenvalues_inter_td_exp[i*2*Nks*Nts],idx_td_Ni[2*i+1],0,idx_td_Ni[2*i+0],0,state_size);
    }




    //**************** // omega(t) (ai+ai^d)  ****************//
    for(unsigned int i=0;i<N_ctl1;i++)
    {
      diagonal_operator_2_bit(state,&eigenvalues_ctl_exp[i*2*Nts],idx_ctl1[i],0,state_size);
    }



    //**************** V such that V * V^d =I  ****************//
    for(unsigned int i=0;i<Nt;i++)// Vi^d
    {
      basis_transformation_2_bit_real(state,&lin_trafo_mat[0],i,0,state_size);
    }
    for(unsigned int i=0;i<Nk;i++)// Vi^d
    {
      basis_transformation_2_bit_real(state,&lin_trafo_mat[0],i,2*Nt,state_size);
    }

    #if (ADIABATIC==0)
    // //**************** exp(-i tau/2 DrivetermFlux) ****************//

    for(unsigned int i=0;i<N_ctl0;i++)// Vi^d
    {

      basis_transformation_2_bit_complex(state,&lin_trafo_mat_adj[1*offset_Nts],idx_ctl0[i],0,state_size);
      diagonal_operator_2_bit(state,&eigenvalues_ctl_exp[i*2*Nts+max_num_charge_drive*2*Nts],idx_ctl0[i],0,state_size);
      basis_transformation_2_bit_complex(state,&lin_trafo_mat[1*offset_Nts],idx_ctl0[i],0,state_size);

      basis_transformation_2_bit_complex(state,&lin_trafo_mat_adj[2*offset_Nts],idx_ctl0[i],0,state_size);
      diagonal_operator_2_bit(state,&eigenvalues_ctl_exp[i*2*Nts+max_num_flux_drive*2*Nts+max_num_charge_drive*2*Nts],idx_ctl0[i],0,state_size);
      basis_transformation_2_bit_complex(state,&lin_trafo_mat[2*offset_Nts],idx_ctl0[i],0,state_size);
    }
    #endif

    // //**************** exp(-i tau/2 H0) ****************//
    for(unsigned int i=0;i<Nt;i++)
    {

      diagonal_operator_2_bit(state,&eigenvalues_device_exp[i*2*Nts],i,0,state_size);
    }
    for(unsigned int i=0;i<Nk;i++)
    {
      diagonal_operator_2_bit(state,&eigenvalues_device_exp[i*2*Nks+Nt*2*Nts],i,2*Nt,state_size);
    }

    if(((t+0)%1000==0) && (t!=(T_end-1)) && sample){sample_time_evolution(state,dt,tar_idx);}

  }
  if(N_ctl0==1)
  {
    time_offset_ctl0+=T_end*tau;
  }
  if(N_ctl1==1)
  {
    time_offset_ctl1+=T_end*tau;
  }
  time_total+=T_end*tau;
}
#endif

#ifndef FD_CPU
#define FD_CPU

class JUSQUACE_FD: public DEVICE_DATA
{
  private:
    double flux_amp0=0;
    double flux_amp1=0;
    double charge_amp=0;

    double *para_td_Ni_hold=nullptr;

    double *statecopy=nullptr;
    double *matrix=nullptr;
    double *matrix_adjoint=nullptr;
    double *eigenvalues_matrix_exp=nullptr;

  public:
    double *eigenvalues_device=nullptr;
    double *eigenvalues_matrix=nullptr;

    //The subroutine which initialises the data structure for the diagonal terms of the Hamiltonian.
    inline void initial_fd_data_structures(void)
    {

      for(unsigned int i=0;i<Nt;i++)
      {
        for(unsigned int j=0;j<Nts;j++)
        {
          eigenvalues_device[j+i*Nts]=freq_at_flux_zero[i]*j+(anharm_at_flux_zero[i]*0.5)*j*(j-1);
        }
      }

      for(unsigned int i=0;i<Nk;i++)
      {
        for(unsigned int j=0;j<Nks;j++)
        {

          eigenvalues_device[j+i*Nks+Nt*Nts]=res_freq[i]*j;
        }
      }
    }
  private:
    //***A set of subroutines which implement algebraic operations to simulate the time evolution of a quantum system***//
    //The subroutine computes the adjoint of a complex matrix.
    inline void adjoint(const double *matrix,double *matrix_adjoint,const unsigned int N,const unsigned int offset);
    //The subroutine multiples the diagonal matrix exp(-i matrixelements[i] tau) with the state vector elements state[2*i+0] and state[2*i+1]
    inline void exp_diagonal_matrix_vector_multiplication(double *state,const double *matrixelements,const unsigned int N,const double tau);
    //The subroutine multiples a non-diagonal matrix with the vector state
    inline void full_matrix_vector_multiplication(double *vector,const double *vectorcopy,const double *matrix,const unsigned int N);

    //***A set of subroutines which update a Hermitian matrix (Hamiltonian) to simulate the time evolution of a quantum system***//
    //The subroutine updates all diagonal terms of the Hamiltonian matrix
    inline void update_diagonal_terms(double *matrix,double *eigenvalues,const unsigned int n_start,const unsigned int state_size);
    //The subroutine updates all interaction terms of the Hamiltonian matrix
    inline void update_interaction_terms(double *matrix,const double lamda,const unsigned int idx_N,const unsigned int idx_K,const unsigned int n_start,const unsigned int k_start,const unsigned int state_size);
    //The subroutine updates all charge dricing terms of the Hamiltonian matrix
    inline void update_charge_driving_terms(double *matrix,const double lamda,const unsigned int idx_N,const unsigned int n_start,const unsigned int state_size);
    //The subroutine updates the first part of the flux driving terms of the Hamiltonian matrix
    inline void update_flux_driving_terms_part_one(double *matrix,const double lamda,const unsigned int idx_N,const unsigned int n_start,const unsigned int state_size);
    //The subroutine updates the second part of the flux driving terms of the Hamiltonian matrix
    inline void update_flux_driving_terms_part_two(double *matrix,const double lamda,const unsigned int idx_N,const unsigned int n_start,const unsigned int state_size);
    // The subroutine initial_matrix(const double dt) initialises the Hamiltonian at time dt.
    inline void initial_matrix(const double dt);

  public:
    //***A set of subroutines which allocate and free memory and update time-dependent data strucutres to simulate the time evolution of a quantum system***//
    //The subroutine updates all time-dependent matrix elements and parameters for a fixed time step dt.
    inline void update_time_dependent_data(const double dt)
    {
      double flux=0;
      double flux_derivative=0;
      double freq_time_dep=0;
      double alpha=0;
      double tri_factor=0;
      unsigned int idx=0;
      // flux control time dependencies
      for(unsigned int i=0;i<N_ctl0;i++)
      {
        flux=flux_ctl(para_ctl0[i],dt);
        idx=idx_ctl0[i];
        // Update the eigenvalues of the anharmonic oscillators
        if(bool_eigenvalues_with_high_accuracy==0)
        {
          freq_time_dep=flux_freq(freq_at_flux_zero[idx],flux,flux_offset[idx],asymmetry_factor_qubit[idx]);
          alpha=anharm_at_flux_zero[idx];
        }
        else
        {
          freq_time_dep=flux_freq_with_high_accuracy(EC[idx],EJ[2*idx+0],EJ[2*idx+1],flux,flux_offset[idx]);
          alpha=flux_freq_anharmonicity_with_high_accuracy(EC[idx],EJ[2*idx+0],EJ[2*idx+1],flux,flux_offset[idx]);
        }
        for(unsigned int j=0;j<Nts;j++)
        {
          // eigenvalues_qubit[idx][j]=freq_time_dep*j+(alpha*0.5)*j*(j-1);
          eigenvalues_device[j+idx*Nts]=freq_time_dep*j+(alpha*0.5)*j*(j-1);
        }

        // Update the time-dependent interaction strength parameters
        for(unsigned int j=0;j<Ni_td;j++)
        {
          if(idx_ctl0[i]==idx_td_Ni_ctl0[j])
          {
            para_td_Ni[j]=time_dependent_interaction_strength(para_td_Ni_hold[j],flux+flux_offset[idx_td_Ni_ctl0[j]],idx_td_Ni_ctl0[j]);
          }
        }

        // Update the time-dependent amplitudes for the flux driving term
        flux=flux+flux_offset[idx];
        flux_derivative=flux_ctl_der(para_ctl0[i],dt);
        double hold=0;
        tri_factor=0;
        hold=cos(flux);//in units of 1/pi not 1/2pi
        tri_factor+=hold*hold;
        hold=sin(flux);//in units of 1/pi not 1/2pi
        tri_factor+=(asymmetry_factor_qubit[idx]*asymmetry_factor_qubit[idx]*hold*hold);

        flux_amp0=amp_factor[2*idx+0]*pow(tri_factor, -0.875)*2*flux_derivative;
        flux_amp1=(amp_factor[2*idx+1]*sin(2*flux)*2*flux_derivative)/tri_factor;


      }
      // charge control time dependencies
      for(unsigned int i=0;i<N_ctl1;i++)
      {
        charge_amp=charge_ctl(para_ctl1[i],dt);
      }
    }
    // The subroutine diagonalization(const double dt) calls initial_matrix(const double dt) and LAPACK for the full diagonalisation.
    inline void diagonalization(const double dt,const char opt='V');
    //***The main subroutine which simulates the time evolution of a quantum system***//
    // The subroutine fd_transmon_basis(double *state) implements a loop over the discreate time steps (t+0.5)*tau, where tau is the time-grid parameters
    // and implements the matrix-vector products state_new=matrix((t+0.5)*tau) exp(-i\tau eigenvalues((t+0.5)*tau)) matrix_adjoint((t+0.5)*tau) state_old.
    void fd_time_evolution(double *state);

    JUSQUACE_FD() = default ;

    JUSQUACE_FD(vector<double> qubit_freq, vector<double> anharm,vector<double> flux_offset_in,vector<double> asymmetry_factor,vector<double> resonator_freq ,unsigned int num_transmons,unsigned int num_resonators,
      unsigned int num_flux_drives,unsigned int num_charge_drives,
      unsigned int num_interactions,vector<unsigned int> interaction_idx,vector<double> interaction_para,
      double time_grid_parameter,unsigned int num_iter_steps,
      unsigned int num_td_interactions=0,vector<unsigned int> interaction_td_idx={},vector<unsigned int> interaction_td_idx_ctl0={},vector<double> interaction_td_para={},
      vector<double> EC_input={},vector<double> EJ_input={},unsigned int num_qubit=0):DEVICE_DATA(qubit_freq,anharm,flux_offset_in,asymmetry_factor,resonator_freq,num_transmons,num_resonators,num_flux_drives,num_charge_drives,num_interactions,interaction_idx,interaction_para,time_grid_parameter,num_iter_steps,num_td_interactions,interaction_td_idx,interaction_td_idx_ctl0,interaction_td_para,EC_input,EJ_input,num_qubit)
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not allocated" << endl;}};
      if(visible==true){cerr << "FD constructor" << endl;}


      const string option="allocate";
      if(option=="allocate")
      {
        if(state_size!=0)
        {
          matrix = new double[2*state_size*state_size]();
          matrix_adjoint = new double[2*state_size*state_size]();
          eigenvalues_matrix = new double[1*state_size]();
          eigenvalues_matrix_exp = new double[2*state_size]();
          statecopy = new double[2*state_size]();
        }
        else
        {
          warning("matrix data structure");
        }

        if((Nt!=0 || Nk!=0))
        {
          eigenvalues_device=new double[Nt*Nts+Nk*Nks]();
        }
        else
        {
          warning("eigenvalues device data structure");
        }

        if((Ni_td!=0))
        {
          para_td_Ni_hold=new double[Ni_td]();
        }
        else
        {
          warning("aux. time-dependent interaction strength data structure");
        }
      }

      for(unsigned int i=0;i<Ni_td;i++){para_td_Ni_hold[i]=para_td_Ni[i];}

      initial_fd_data_structures();
    }

    ~JUSQUACE_FD()
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not freed" << endl;}};
      if(visible==true){cerr << "FD destructor" << endl;}

      const string option="free";
      if(option=="free")
      {
        if(state_size!=0)
        {

          delete [] matrix;
          matrix=nullptr;
          delete [] matrix_adjoint;
          matrix_adjoint=nullptr;
          delete [] eigenvalues_matrix;
          matrix_adjoint=nullptr;
          delete [] eigenvalues_matrix_exp;
          eigenvalues_matrix_exp=nullptr;
          delete [] statecopy;
          statecopy=nullptr;
        }
        else
        {
          warning("matrix data structure");
        }

        if( (Nt!=0 || Nk!=0))
        {
          delete [] eigenvalues_device;
          eigenvalues_device=nullptr;
        }
        else
        {
          warning("eigenvalues device data structure");
        }

        if((Ni_td!=0))
        {
          delete [] para_td_Ni_hold;
          para_td_Ni_hold=nullptr;
        }
        else
        {
          warning("aux. time-dependent interaction strength data structure");
        }
      }
    }

    JUSQUACE_FD(const JUSQUACE_FD & other):DEVICE_DATA(other)
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not copied" << endl;}};
      if(visible==true){cerr << "FD copy constructor" << endl;}

      const string option="copy";
      if(option=="copy")
      {
        if(state_size!=0)
        {
          matrix = new double[2*state_size*state_size]();
          memcpy(matrix,other.matrix,sizeof(double)*2*state_size*state_size);
          matrix_adjoint = new double[2*state_size*state_size]();
          memcpy(matrix_adjoint,other.matrix_adjoint,sizeof(double)*2*state_size*state_size);
          eigenvalues_matrix = new double[1*state_size]();
          memcpy(eigenvalues_matrix,other.eigenvalues_matrix,sizeof(double)*1*state_size);
          eigenvalues_matrix_exp = new double[2*state_size]();
          memcpy(eigenvalues_matrix_exp,other.eigenvalues_matrix_exp,sizeof(double)*2*state_size);
          statecopy = new double[2*state_size]();
          memcpy(statecopy,other.statecopy,sizeof(double)*2*state_size);
        }
        else
        {
          warning("matrix data structure");
        }

        if((Nt!=0 || Nk!=0))
        {
          eigenvalues_device=new double[Nt*Nts+Nk*Nks]();
          memcpy(eigenvalues_device,other.eigenvalues_device,sizeof(double)*(Nt*Nts+Nk*Nks));
        }
        else
        {
          warning("eigenvalues device data structure");
        }

        if((Ni_td!=0))
        {
          para_td_Ni_hold=new double[Ni_td]();
          memcpy(para_td_Ni_hold,other.para_td_Ni_hold,sizeof(double)*1*Ni_td);
        }
        else
        {
          warning("aux. time-dependent interaction strength data structure");
        }
      }

    }

    JUSQUACE_FD &operator=(const JUSQUACE_FD & other)
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not copied/assignment" << endl;}};
      if(visible==true){cerr << "FD copy assignment constructor" << endl;}

      if(this!=&other)
      {
        DEVICE_DATA::operator=(other);

        const string option0="free";
        if(option0=="free")
        {
          if(state_size!=0)
          {

            delete [] matrix;
            matrix=nullptr;
            delete [] matrix_adjoint;
            matrix_adjoint=nullptr;
            delete [] eigenvalues_matrix;
            matrix_adjoint=nullptr;
            delete [] eigenvalues_matrix_exp;
            eigenvalues_matrix_exp=nullptr;
            delete [] statecopy;
            statecopy=nullptr;
          }
          else
          {
            warning("matrix data structure");
          }

          if( (Nt!=0 || Nk!=0))
          {
            delete [] eigenvalues_device;
            eigenvalues_device=nullptr;
          }
          else
          {
            warning("eigenvalues device data structure");
          }

          if((Ni_td!=0))
          {
            delete [] para_td_Ni_hold;
            para_td_Ni_hold=nullptr;
          }
          else
          {
            warning("aux. time-dependent interaction strength data structure");
          }
        }

        const string option1="copy";
        if(option1=="copy")
        {
          if(state_size!=0)
          {
            matrix = new double[2*state_size*state_size]();
            memcpy(matrix,other.matrix,sizeof(double)*2*state_size*state_size);
            matrix_adjoint = new double[2*state_size*state_size]();
            memcpy(matrix_adjoint,other.matrix_adjoint,sizeof(double)*2*state_size*state_size);
            eigenvalues_matrix = new double[1*state_size]();
            memcpy(eigenvalues_matrix,other.eigenvalues_matrix,sizeof(double)*1*state_size);
            eigenvalues_matrix_exp = new double[2*state_size]();
            memcpy(eigenvalues_matrix_exp,other.eigenvalues_matrix_exp,sizeof(double)*2*state_size);
            statecopy = new double[2*state_size]();
            memcpy(statecopy,other.statecopy,sizeof(double)*2*state_size);
          }
          else
          {
            warning("matrix data structure");
          }

          if((Nt!=0 || Nk!=0))
          {
            eigenvalues_device=new double[Nt*Nts+Nk*Nks]();
            memcpy(eigenvalues_device,other.eigenvalues_device,sizeof(double)*(Nt*Nts+Nk*Nks));
          }
          else
          {
            warning("eigenvalues device data structure");
          }

          if((Ni_td!=0))
          {
            para_td_Ni_hold=new double[Ni_td]();
            memcpy(para_td_Ni_hold,other.para_td_Ni_hold,sizeof(double)*1*Ni_td);
          }
          else
          {
            warning("aux. time-dependent interaction strength data structure");
          }
        }
      }
      return *this;
    }

    JUSQUACE_FD(JUSQUACE_FD && other):DEVICE_DATA(move(other))
    {
      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not moved" << endl;}};
      if(visible==true){cerr << "FD move constructor" << endl;}

      const string option="move";
      if(option=="move")
      {
        if(state_size!=0)
        {
          matrix = other.matrix;
          other.matrix=nullptr;
          matrix_adjoint = other.matrix_adjoint;
          other.matrix_adjoint =nullptr;
          eigenvalues_matrix = other.eigenvalues_matrix;
          other.eigenvalues_matrix = nullptr;
          eigenvalues_matrix_exp = other.eigenvalues_matrix_exp;
          other.eigenvalues_matrix_exp = nullptr;
          statecopy = other.statecopy;
          other.statecopy=nullptr;
        }
        else
        {
          warning("matrix data structure");
        }

        if((Nt!=0 || Nk!=0))
        {
          eigenvalues_device=other.eigenvalues_device;
          eigenvalues_device=nullptr;
        }
        else
        {
          warning("eigenvalues device data structure");
        }

        if((Ni_td!=0))
        {
          para_td_Ni_hold=other.para_td_Ni_hold;
          other.para_td_Ni_hold=nullptr;
        }
        else
        {
          warning("aux. time-dependent interaction strength data structure");
        }
      }

    }

    JUSQUACE_FD &operator=(JUSQUACE_FD && other)
    {

      bool visible=false;
      auto warning = [visible](const string str) {if(visible==true){ cerr << "Dynamic memory for "+str+" not moved/assignment" << endl;}};
      if(visible==true){cerr << "FD move assignment constructor" << endl;}

      if(this!=&other)
      {
        DEVICE_DATA::operator=(move(other));

        const string option0="free";
        if(option0=="free")
        {
          if(state_size!=0)
          {

            delete [] matrix;
            matrix=nullptr;
            delete [] matrix_adjoint;
            matrix_adjoint=nullptr;
            delete [] eigenvalues_matrix;
            matrix_adjoint=nullptr;
            delete [] eigenvalues_matrix_exp;
            eigenvalues_matrix_exp=nullptr;
            delete [] statecopy;
            statecopy=nullptr;
          }
          else
          {
            warning("matrix data structure");
          }

          if( (Nt!=0 || Nk!=0))
          {
            delete [] eigenvalues_device;
            eigenvalues_device=nullptr;
          }
          else
          {
            warning("eigenvalues device data structure");
          }

          if((Ni_td!=0))
          {
            delete [] para_td_Ni_hold;
            para_td_Ni_hold=nullptr;
          }
          else
          {
            warning("aux. time-dependent interaction strength data structure");
          }
        }

        const string option1="move";
        if(option1=="move")
        {
          if(state_size!=0)
          {
            matrix = other.matrix;
            other.matrix=nullptr;
            matrix_adjoint = other.matrix_adjoint;
            other.matrix_adjoint =nullptr;
            eigenvalues_matrix = other.eigenvalues_matrix;
            other.eigenvalues_matrix = nullptr;
            eigenvalues_matrix_exp = other.eigenvalues_matrix_exp;
            other.eigenvalues_matrix_exp = nullptr;
            statecopy = other.statecopy;
            other.statecopy=nullptr;
          }
          else
          {
            warning("matrix data structure");
          }

          if((Nt!=0 || Nk!=0))
          {
            eigenvalues_device=other.eigenvalues_device;
            other.eigenvalues_device=nullptr;
          }
          else
          {
            warning("eigenvalues device data structure");
          }

          if((Ni_td!=0))
          {
            para_td_Ni_hold=other.para_td_Ni_hold;
            other.para_td_Ni_hold=nullptr;
          }
          else
          {
            warning("aux. time-dependent interaction strength data structure");
          }
        }
      }
      return *this;
    }
};

inline void JUSQUACE_FD::adjoint(const double *matrix,double *matrix_adjoint,const unsigned int N,const unsigned int offset)
{
  unsigned int idx,idy;
  for(unsigned int n=0;n<N;n++)
  {
    for(unsigned int m=0;m<N;m++)
    {
      idx=2*m+offset*n;
      idy=2*n+offset*m;
      matrix_adjoint[idx]=matrix[idy];
      matrix_adjoint[idx+1]=-matrix[idy+1];
    }
  }
}

inline void JUSQUACE_FD::exp_diagonal_matrix_vector_multiplication(double *state,const double *matrixelements,const unsigned int N,const double tau)
{
  double ex,ey,sx,sy;
  for(unsigned int n=0;n<N;n++)
  {
    ex=cos(-tau*matrixelements[n]);
    ey=sin(-tau*matrixelements[n]);
    sx=state[2*n];
    sy=state[2*n+1];
    state[2*n]  =sx*ex - ey*sy;
    state[2*n+1]=sx*ey + ex*sy;
  }
}

inline void JUSQUACE_FD::full_matrix_vector_multiplication(double *vector,const double *vectorcopy,const double *matrix,const unsigned int N)
{
  unsigned int idx_0,idx_1,offset=2*N;
  double hold_real,hold_img;
  for(unsigned int i=0;i<N;i++)
  {
    hold_real=0;
    hold_img=0;
    for(unsigned int j=0;j<N;j++)
    {
      idx_0=2*j+offset*i;
      idx_1=2*j;
      hold_real += matrix[idx_0]*vectorcopy[idx_1]   - vectorcopy[idx_1+1]*matrix[idx_0+1];
      hold_img  += matrix[idx_0]*vectorcopy[idx_1+1] + vectorcopy[idx_1]*matrix[idx_0+1];
    }
    idx_1=2*i;
    vector[idx_1]=hold_real;
    vector[idx_1+1]=hold_img;
  }
}

inline void JUSQUACE_FD::update_diagonal_terms(double *matrix,double *eigenvalues,const unsigned int n_start,const unsigned int state_size)
{
  unsigned int shift=n_start;
  unsigned int offset = 2*state_size;
  unsigned int idx_matrix=0;
  unsigned int idx_n=0;

  for(unsigned int i=0;i<state_size;i++)
  {
    idx_matrix=2*i+offset*i;
    //**********************************************//
    idx_n= (i >> shift) & 3 ;
    //**********************************************//
    matrix[idx_matrix]+=eigenvalues[idx_n];
  }
}

inline void JUSQUACE_FD::update_interaction_terms(double *matrix,const double lamda,const unsigned int idx_N,const unsigned int idx_K,const unsigned int n_start,const unsigned int k_start,const unsigned int state_size)
{
  unsigned int shift_n=n_start;
  unsigned int shift_k=k_start;
  unsigned int idx_matrix=0;
  unsigned int idx_n=0;
  unsigned int idx_k=0;
  unsigned int offset = 2*state_size;
  unsigned int offset_n= 1 << (2*idx_N);
  unsigned int offset_k= 1 << (2*idx_K);

  for(unsigned int i=0;i<state_size;i++)
  {
    //**********************************************//
    idx_n= (i >> shift_n) & 3 ;
    idx_k= (i >> shift_k) & 3 ;
    //**********************************************//
    if(idx_n>0 && idx_k>0)
    {
      idx_matrix=2*i+offset*(i-offset_n-offset_k);
      matrix[idx_matrix]+=lamda*sqrt(idx_n)*sqrt(idx_k);
    }
    if(idx_n>0 && idx_k<3)
    {
      idx_matrix=2*i+offset*(i-offset_n+offset_k);
      matrix[idx_matrix]+=lamda*sqrt(idx_n)*sqrt(idx_k+1);
    }
    if(idx_n<3 && idx_k>0)
    {
      idx_matrix=2*i+offset*(i+offset_n-offset_k);
      matrix[idx_matrix]+=lamda*sqrt(idx_n+1)*sqrt(idx_k);
    }
    if(idx_n<3 && idx_k<3)
    {
      idx_matrix=2*i+offset*(i+offset_n+offset_k);
      matrix[idx_matrix]+=lamda*sqrt(idx_n+1)*sqrt(idx_k+1);
    }
  }
}

inline void JUSQUACE_FD::update_charge_driving_terms(double *matrix,const double lamda,const unsigned int idx_N,const unsigned int n_start,const unsigned int state_size)
{
  unsigned int shift_n=n_start;
  unsigned int idx_matrix=0;
  unsigned int idx_n=0;
  unsigned int offset = 2*state_size;
  unsigned int offset_n= 1 << (2*idx_N);


  for(unsigned int i=0;i<state_size;i++)
  {
    //**********************************************//
    idx_n= (i >> shift_n) & 3 ;
    //**********************************************//
    if(idx_n>0)
    {
      idx_matrix=2*i+offset*(i-offset_n);
      matrix[idx_matrix]+=lamda*sqrt(idx_n);
    }
    if(idx_n<3)
    {
      idx_matrix=2*i+offset*(i+offset_n);
      matrix[idx_matrix]+=lamda*sqrt(idx_n+1);
    }
  }
}

inline void JUSQUACE_FD::update_flux_driving_terms_part_one(double *matrix,const double lamda,const unsigned int idx_N,const unsigned int n_start,const unsigned int state_size)
{
  unsigned int shift_n=n_start;
  unsigned int idx_matrix=0;
  unsigned int idx_n=0;
  unsigned int offset = 2*state_size;
  unsigned int offset_n= 1 << (2*idx_N);


  for(unsigned int i=0;i<state_size;i++)
  {
    //**********************************************//
    idx_n= (i >> shift_n) & 3 ;
    //**********************************************//
    if(idx_n>0)
    {
      idx_matrix=2*i+1+offset*(i-offset_n);
      matrix[idx_matrix]+=-lamda*sqrt(idx_n-0);
    }
    if(idx_n<3)
    {
      idx_matrix=2*i+1+offset*(i+offset_n);
      matrix[idx_matrix]+=lamda*sqrt(idx_n+1);
    }
  }
}

inline void JUSQUACE_FD::update_flux_driving_terms_part_two(double *matrix,const double lamda,const unsigned int idx_N,const unsigned int n_start,const unsigned int state_size)
{
  unsigned int shift_n=n_start;
  unsigned int idx_matrix=0;
  unsigned int idx_n=0;
  unsigned int offset = 2*state_size;
  unsigned int offset_n= 1 << (2*idx_N);


  for(unsigned int i=0;i<state_size;i++)
  {
    //**********************************************//
    idx_n= (i >> shift_n) & 3 ;
    //**********************************************//
    if(idx_n>1)
    {
      idx_matrix=2*i+1+offset*(i-2*offset_n);
      matrix[idx_matrix]+=-lamda*sqrt((idx_n-1)*(idx_n-0));
    }
    if(idx_n<2)
    {
      idx_matrix=2*i+1+offset*(i+2*offset_n);
      matrix[idx_matrix]+=lamda*sqrt((idx_n+2)*(idx_n+1));
    }
  }
}

inline void JUSQUACE_FD::initial_matrix(const double dt)
{
  for(unsigned int i=0;i<N_ctl0;i++)
  {

    update_flux_driving_terms_part_one(matrix,flux_amp0,idx_ctl0[i],2*idx_ctl0[i],state_size);
    update_flux_driving_terms_part_two(matrix,flux_amp1,idx_ctl0[i],2*idx_ctl0[i],state_size);
  }
  for(unsigned int i=0;i<N_ctl1;i++)
  {

    update_charge_driving_terms(matrix,charge_amp,idx_ctl1[i],2*idx_ctl1[i],state_size);
  }
  for(unsigned int i=0;i<Nt;i++)
  {
    update_diagonal_terms(matrix,&eigenvalues_device[i*Nts],2*i,state_size);
  }
  for(unsigned int i=0;i<Nk;i++)
  {
    update_diagonal_terms(matrix,&eigenvalues_device[i*Nks+Nt*Nts],2*(i+Nt),state_size);
  }
  for(unsigned int i=0;i<Ni;i++)
  {

    update_interaction_terms(matrix,para_Ni[i],idx_Ni[2*i], idx_Ni[2*i+1],2*idx_Ni[2*i],2*idx_Ni[2*i+1],state_size);
  }
  for(unsigned int i=0;i<Ni_td;i++)
  {

    update_interaction_terms(matrix,para_td_Ni[i],idx_td_Ni[2*i], idx_td_Ni[2*i+1],2*idx_td_Ni[2*i],2*idx_td_Ni[2*i+1],state_size);
  }
}

inline void JUSQUACE_FD::diagonalization(const double dt,const char opt)
{
  memset(matrix, 0, sizeof(double)*2*state_size*state_size);
  memset(matrix_adjoint, 0, sizeof(double)*2*state_size*state_size);
  memset(eigenvalues_matrix, 0, sizeof(double)*1*state_size);
  memset(eigenvalues_matrix_exp, 0, sizeof(double)*2*state_size);
  initial_matrix(dt);
  LAPACKE_zheev( LAPACK_ROW_MAJOR, opt, 'U', state_size, &matrix[0], state_size, &eigenvalues_matrix[0]);
  adjoint(matrix,matrix_adjoint,state_size,2*state_size);
}

void JUSQUACE_FD::fd_time_evolution(double *state)
{

  double dt=0;
  for(unsigned int i=0;i<Ni_td;i++){para_td_Ni[i]=time_dependent_interaction_strength(para_td_Ni_hold[i],0+flux_offset[idx_td_Ni_ctl0[i]],idx_td_Ni_ctl0[i]);}

  for(unsigned int t=T_start;t<T_end;t++)
  {

    dt=(t+0.5)*tau;

    update_time_dependent_data(dt);

    diagonalization(dt);

    memset(statecopy, 0, sizeof(double)*2*state_size);

    memcpy(statecopy, state, sizeof(double)*2*state_size);

    full_matrix_vector_multiplication(state,statecopy,matrix_adjoint,state_size);

    exp_diagonal_matrix_vector_multiplication(state,eigenvalues_matrix,state_size,tau);

    memcpy(statecopy, state, sizeof(double)*2*state_size);

    full_matrix_vector_multiplication(state,statecopy,matrix,state_size);
  }
}

#endif
