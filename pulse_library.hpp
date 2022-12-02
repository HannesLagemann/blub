#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>
#include <string>
#include <vector>
#include <float.h>
#include <math.h>
#include <time.h>
#include <chrono>
#include <cstring>
#include <sys/time.h>
#include <random>
#define lapack_complex_float float
#define lapack_complex_double double
#ifdef __APPLE__
#include </opt/homebrew/Cellar/openblas/0.3.15_1/include/cblas.h>
#include </opt/homebrew/Cellar/openblas/0.3.15_1/include/lapacke.h>
#endif
#ifdef __linux__
#include <mkl.h>
#include <mkl_lapacke.h>
//#include <cblas.h>
//#include <lapacke.h>
#endif

//#pragma warning disable 186
#pragma GCC diagnostic ignored "-Wtype-limits"
#pragma GCC diagnostic ignored "-Wpedantic"

using namespace std;

#ifndef PI
#define PI
static const double pi=3.1415926535897932384626433832795028841971693993751058209749445923;
#endif

#ifndef SQ_PULSE_FUNC
#define SQ_PULSE_FUNC
#if (SQ_PULSE==0)
//GAUSS
inline double single_qubit_charge_ctl(const vector<double> &ctl_parameters,const double dt,const double offset=0)
{

  //ctl_parameters[3]=((-ctl_parameters[1]*offset)+ctl_parameters[3]);
  double phase=((-ctl_parameters[1]*offset)+ctl_parameters[3]);
  double ttss=2*ctl_parameters[6]*ctl_parameters[6]; // 2*sigma*sigma
  double hold=(-ctl_parameters[0]*ctl_parameters[0])/(ttss*4);
  double hold_Gauss=0;

  hold_Gauss=(dt-ctl_parameters[0]*0.5)*(dt-ctl_parameters[0]*0.5);
  hold_Gauss=-hold_Gauss/ttss;
  hold_Gauss=exp(hold_Gauss);

  double hold_der=(ctl_parameters[2]*hold_Gauss)/(1-exp(hold)); // save partial result for the derivertive

  hold_Gauss=(hold_Gauss-exp(hold))/(1-exp(hold));
  hold_Gauss=ctl_parameters[2]*hold_Gauss;
  hold_Gauss=hold_Gauss*cos(ctl_parameters[1]*dt-phase);

  double hold_Gauss_Dot=0;
  hold_Gauss_Dot=-2*((dt-ctl_parameters[0]*0.5)/ttss)*hold_der;
  hold_Gauss_Dot=ctl_parameters[4]*hold_Gauss_Dot;
  hold_Gauss_Dot=hold_Gauss_Dot*sin(ctl_parameters[1]*dt-phase);

  //return (hold_Gauss+hold_Gauss_Dot);
  return 100;
}

#endif
#endif

#ifndef TQ_PULSE_FUNC
#define TQ_PULSE_FUNC

#if (TQ_PULSE==0)
//UMP
inline double two_qubit_flux_ctl(const vector<double> &ctl_parameters,const double dt)
{

  //cerr << ctl_parameters[2]/(pi) << endl;
  double L0=ctl_parameters[0];
  double sigma=sqrt(2)*ctl_parameters[1];
  double hold=0;
  double t=dt-6*sigma;
  hold+=erf(t/sigma);
  hold-=erf((t-L0)/sigma);
  hold=0.5*ctl_parameters[2]*hold;
  return hold;
  //return 100;
}

#elif (TQ_PULSE==1)
//BMP
inline double two_qubit_flux_ctl(const vector<double> &ctl_parameters,const double dt)
{

  double alpha=1.0;
  double L0=ctl_parameters[0];
  double L1=(alpha/(alpha+1))*L0;
  double sigma=sqrt(2)*ctl_parameters[1];
  double hold=0;
  double t=dt-6*sigma;
  hold+=erf(t/sigma);
  hold-=erf((t-L1)/sigma) ;
  hold-=alpha*erf((t-L1)/sigma);
  hold+=alpha*erf((t-L0)/sigma);
  hold=0.5*ctl_parameters[2]*hold;
  return hold;
}



#elif (TQ_PULSE==2)
//SEP
inline double two_qubit_flux_ctl(const vector<double> &ctl_parameters,const double dt,const double offset=0)
{


  double phase=((-ctl_parameters[1]*offset)+ctl_parameters[3]);
  double T_flat=ctl_parameters[6]-2*ctl_parameters[0]-ctl_parameters[5];
  double T_inc=ctl_parameters[0];
  double T_dec=T_inc+T_flat;
  double T_buf=ctl_parameters[6]-ctl_parameters[5];
  double hold=pi/(2*T_inc);
  double drive_value=0;


  if(dt<=T_inc)
  {
    drive_value=(ctl_parameters[2]*sin(hold*dt)*(cos(ctl_parameters[1]*dt-phase)));
    return drive_value;
  }
  if(T_inc<dt && dt<T_dec)
  {
    drive_value=(ctl_parameters[2]*(cos(ctl_parameters[1]*dt-phase)));
    return drive_value;

  }
  if(dt>=T_dec && dt<=T_buf)
  {
    drive_value=(ctl_parameters[2]*sin(pi*0.5+hold*(dt-T_dec))*(cos(ctl_parameters[1]*dt-phase)));
    return drive_value;
  }
  return 0;
}




#endif



#endif

#ifndef TQ_PULSE_FUNC_DER
#define TQ_PULSE_FUNC_DER


#if (ADIABATIC==0)

inline double two_qubit_flux_ctl_der(const vector<double> &ctl_parameters,const double dt,const double offset=0)
{
  const double tau=0.001;
  double hold=0;
  hold+=two_qubit_flux_ctl(ctl_parameters,(dt+1*tau));
  hold-=two_qubit_flux_ctl(ctl_parameters,(dt-1*tau));
  hold=hold/(2*tau);
  return hold;
  //return 10;
}

#else

inline double two_qubit_flux_ctl_der(const vector<double> &ctl_parameters,const double dt,const double offset=0)
{
  return 0;
}

#endif
#endif
