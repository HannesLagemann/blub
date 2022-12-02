
template<class T>
void ArchitectureII(T &Data,const double time_grid_parameter=0.001)
{

  // //************* System size ***************//

  unsigned int num_transmons=2;
  unsigned int num_resonators=1;
  unsigned int num_qubits=2;

  //************* System Parameter ***************//

  vector<double> EC={6.71217898201047,6.5125182190863};
  vector<double> EJ={19.7280631289123,19.7280631289123*2.0,30.2646640224655,30.2646640224655*2.0};
  vector<double> qubit_freq={2*pi*4.1999396,2*pi*5.2000022};
  vector<double> anharm={-0.32005948*2*pi,-0.29500447*2*pi};
  vector<double> asymmetry_factor={0.5,0.5};
  vector<double> flux_offset={0,0};
  vector<double> phase_shifts={0,0};// you can remove the dynamic phase exp(-i T omega) with these factors (example: gate optimization)
  vector<double> resonator_freq={2*pi*45};

  //************* Set the number of control pulses (flux and charge) ***************//

  unsigned int num_flux_drives=1;
  unsigned int num_charge_drives=1;

  //************* Couple qubit 0 to qubit 2 and qubit 1 to qubit 2 ***************//

  unsigned int num_interactions=1;
  vector<unsigned int> interaction_idx(2*num_interactions,0);
  vector<double> interaction_para(num_interactions,0);
  const double G=0.300;//GHz‚

  interaction_para[0]=2*pi*G*effective_interaction_strength_factor_tunable_transmon(EC[0],EJ[2*0+0],EJ[2*0+1],2*flux_offset[0]);
  interaction_idx[0]=0;
  interaction_idx[1]=2;

  unsigned int num_interactions_time_dependent=1;
  vector<unsigned int> interaction_idx_time_dependent_ctl0(num_interactions_time_dependent,0);
  vector<unsigned int> interaction_idx_time_dependent(2*num_interactions_time_dependent,0);
  vector<double> interaction_para_time_dependent(num_interactions_time_dependent,0);

  interaction_para_time_dependent[0]=2*pi*G;
  interaction_idx_time_dependent[0]=1;
  interaction_idx_time_dependent[1]=2;
  interaction_idx_time_dependent_ctl0[0]=1;

  //************* Set the time grid parameter tau and the number of iteration steps ***************//

  //double time_grid_parameter = 0.001;
  unsigned int num_iter_steps =125000;//

  //************* Feed the class with your input data ***************//

  Data= T(qubit_freq, anharm,flux_offset,asymmetry_factor,resonator_freq,num_transmons,num_resonators,
  num_flux_drives,num_charge_drives,
  num_interactions,interaction_idx,interaction_para,
  time_grid_parameter,num_iter_steps,
  num_interactions_time_dependent,interaction_idx_time_dependent,interaction_idx_time_dependent_ctl0,interaction_para_time_dependent,
  EC,EJ,num_qubits);


  //************* Set aux. parameters ***************//

  Data.chip_name="ARCII";
  Data.lamda_trans_shifted[0]=2*pi*4.19560457393345;
  Data.lamda_trans_shifted[1]=2*pi*5.1946805628581;


  //************* Define the flux control pulse ***************//
  const bool gate_set[2+1]={false,false,true};

  // A sine envelope e(t) (SEP) with the rise and fall time Data.para_ctl0[idx_CZ01][1]
  // and the cosine function cos(Data.para_ctl0[idx_CZ01][1] t) form the pulse
  // flux(t)= Data.para_ctl0[idx_CZ01][2] e(t) cos(Data.para_ctl0[idx_CZ01][1]*t-Data.para_ctl0[idx_CZ01][3])
  // with the pulse duration Data.para_ctl0[idx_CZ01][6], i.e. e(Data.para_ctl0[idx_CZ01][1])=0.
  if(num_charge_drives!=0)
  {
    // Controll charge X_pi/2 rotation qubit 0
    if(gate_set[0]==true)
    {

      unsigned int idx_RX0=0;
      Data.T_end=52250;
      Data.idx_ctl1[idx_RX0]=0;
      //double fac_0=sqrt(sqrt((8*EC[Data.idx_ctl1[0]])/ (EJ[2*Data.idx_ctl1[0]]+EJ[2*Data.idx_ctl1[0]+1]) ))/(2*EC[Data.idx_ctl1[0]]);
      Data.para_ctl1[idx_RX0][0]=Data.T_end*Data.tau;
      Data.para_ctl1[idx_RX0][1]=2*pi*4.19560457393345;//4.19560457393345;//1.0;// freq*2*pi
      Data.para_ctl1[idx_RX0][2]=0.0584;//amp0
      Data.para_ctl1[idx_RX0][3]=0;//0.0;//  phi1
      Data.para_ctl1[idx_RX0][4]=0.0722505531310642; //amp2
      Data.para_ctl1[idx_RX0][5]=0.0; //phi2
      Data.para_ctl1[idx_RX0][6]=12.0823585295687; // T_fall/rise
    }

    //Controll charge X_pi/2 rotation qubit 1
    if(gate_set[1]==true)
    {
      unsigned int idx_RX1=1;
      Data.T_end=52950;
      Data.idx_ctl1[idx_RX1]=1;
      //double fac_1=sqrt(sqrt((8*EC[Data.idx_ctl1[0]])/ (EJ[2*Data.idx_ctl1[0]]+EJ[2*Data.idx_ctl1[0]+1]) ))/(2*EC[Data.idx_ctl1[0]]);
      Data.para_ctl1[idx_RX1][0]=Data.T_end*Data.tau;// T
      Data.para_ctl1[idx_RX1][1]=2*pi*5.1946805628581;//freq*2*pi
      Data.para_ctl1[idx_RX1][2]=0.0654;//amp0
      Data.para_ctl1[idx_RX1][3]=0.0;//0.0;//  phi1
      Data.para_ctl1[idx_RX1][4]=0.07 ; //amp2
      Data.para_ctl1[idx_RX1][5]=0.0; //phi2
      Data.para_ctl1[idx_RX1][6]=10.0; // T_fall/rise
    }
  }

  // A sine envelope e(t) (SEP) with the rise and fall time Data.para_ctl0[idx_CZ01][1]
  // and the cosine function cos(Data.para_ctl0[idx_CZ01][1] t) form the pulse
  // flux(t)= Data.para_ctl0[idx_CZ01][2] e(t) cos(Data.para_ctl0[idx_CZ01][1]*t-Data.para_ctl0[idx_CZ01][3])
  // with the pulse duration Data.para_ctl0[idx_CZ01][6], i.e. e(Data.para_ctl0[idx_CZ01][1])=0.
  string PULSE="SEP";
  if(num_flux_drives!=0)
  {

    if(gate_set[2]==true)
    {
      if(PULSE=="SEP")
      {
        unsigned int idx_CZ01=0;
        Data.T_end=300000;//
        Data.idx_ctl0[idx_CZ01]=0;
        Data.para_ctl0[idx_CZ01][0]=Data.T_end*Data.tau*0.5;//Time rise and fall
        Data.para_ctl0[idx_CZ01][1]=2*pi*45.00;//freq*2*pi
        Data.para_ctl0[idx_CZ01][2]=pi*0.005;//amp0*pi
        Data.para_ctl0[idx_CZ01][3]=0;//phase 0
        Data.para_ctl0[idx_CZ01][4]=0;//amp 1
        Data.para_ctl0[idx_CZ01][5]=5;//phase 1
        Data.para_ctl0[idx_CZ01][6]=Data.T_end*Data.tau;//Total time

        Data.phase_CZ[idx_CZ01][0]=0;
        Data.phase_CZ[idx_CZ01][1]=0;
      }
    }
  }
}
