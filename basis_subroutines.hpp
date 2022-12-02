#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <string.h>
#include <fstream>
#include <float.h>
#ifdef __linux__
#include <experimental/filesystem>
namespace filesystem = std::experimental::filesystem;
//#include <filesystem>
//namespace filesystem = std::filesystem;
#endif
#ifdef __APPLE__
//#include <experimental/filesystem>
#endif
#ifdef MPI_MARCO
#include <mpi.h>
#endif
#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace std;
#ifndef sqrt5
#define sqrt5 2.236067977499789696
#endif

#ifndef BASIS_SUB
#define BASIS_SUB


static const string hline="--------------------------------------------------\n";

inline void warning(string str)
{
  cout << "warning: " << str << endl;
}

template<typename T>
inline void set_state_zero(T &state,const unsigned int size)
{

  for(unsigned int i=0;i<size;i++)
  {
    state[2*i+0]=0;
    state[2*i+1]=0;
  }
}

template<typename T>
inline void print_prob_state(const T &state,const unsigned int size)
{
  double hold=0;
  for(unsigned int i=0;i<size;i++)
  {
    cout << setw(5) << left  << i;
    hold=state[2*i+0]*state[2*i+0]+state[2*i+1]*state[2*i+1];
    cout << setw(40) << left << hold;
    cout << endl;
  }
}

template<typename T>
inline void print_complex_matrix(const T &matrix,const unsigned int end_r,const unsigned int end_c,const unsigned int offset)
{
  cout << "-----------------------------------------------------------------------------------------------------" << endl;
  for(unsigned int c=0;c<end_c;c++)
  {
    string hold="#"+to_string(c)+":";
    cout << setw(5) << left  << hold;
    for(unsigned int r=0;r<end_r;r++)
    {
      hold=to_string(matrix[2*r+offset*c])+"+ i "+to_string(matrix[2*r+1+offset*c]);
      cout << setw(40) << left << hold;
    }
    cout << endl;
  }
  cout << "-----------------------------------------------------------------------------------------------------" << endl;
}

template<typename T>
inline void print_prob_matrix(const T &matrix,const unsigned int num_row,const unsigned int num_col,const unsigned int offset)
{
  unsigned int size_i=num_row;
  unsigned int size_j=num_col/2;
  unsigned int idx;
  double hold;
  cout << "-----------------------------------------------------------------------------------------------------" << endl;
  for(unsigned int i=0;i<size_i;i++)
  {
    string hold_str="#"+to_string(i)+":";
    cout << setw(5) << left << hold_str;
    for(unsigned int j=0;j<size_j;j++)
    {
      idx=2*j+offset*i;
      hold=matrix[idx]*matrix[idx]+matrix[idx+1]*matrix[idx+1];
      cout << setw(40) << left << hold;
    }
    cout << endl;

  }
  cout << "-----------------------------------------------------------------------------------------------------" << endl;
}

template <typename T>
inline void state_copy(const T &state,T &statecopy,const unsigned int size)
{
  unsigned int idx=0;
  for(unsigned int n=0;n<size;n++)
  {
    idx=2*n;
    statecopy[idx]=state[idx];
    statecopy[idx+1]=state[idx+1];
  }
}

template <typename T>
inline void copy_matrix0_to_matrix1(const T &matrix0,T &matrix1,unsigned int size)
{
  for(unsigned int i=0;i<size;i++)
  {
    matrix1[i]=matrix0[i];
  }
}

template <typename T>
inline double norm(const T &state,unsigned int size)
{
  double hold=0;
  for(unsigned int n=0;n<size;n++)
  {
    hold+=state[2*n+0]*state[2*n+0]+state[2*n+1]*state[2*n+1];
  }
  return hold;
}

inline void matXvec(double *matrix,double *vec, unsigned int size)
{
  unsigned int idx=0;
  double *hold=new double[2*size]();
  for(unsigned int i=0;i<size;i++)
  {
    hold[2*i+0]=0;
    hold[2*i+1]=0;
    for(unsigned int j=0;j<size;j++)
    {
      idx=2*j+0+2*size*i;
      hold[2*i+0]+=matrix[idx]*vec[2*i+0] - vec[2*j+1]*matrix[idx+1];
      hold[2*i+1]+=matrix[idx]*vec[2*i+1] + vec[2*j+0]*matrix[idx+1];
    }
  }
  for(unsigned int i=0;i<(2*size);i++)
  {
    vec[i]=hold[i];
  }
  delete[] hold;
}

inline void vector_matrix_multiplication_nxn(double *vector,double *matrix,unsigned int N)
{
  unsigned int idx_0=0,idx_1=0,offset=2*N;
  double hold_real=0,hold_img=0;
  double *vectorcopy=new double[2*N];
  for(unsigned int i=0;i<N;i++)
  {
    vectorcopy[2*i+0]=vector[2*i+0];
    vectorcopy[2*i+1]=vector[2*i+1];
  }
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
  delete[] vectorcopy;
}

inline void vector_exp_matrix_multiplication_nxn(double *vector,double *diagonal_matrix,double scalar,unsigned int size)
{

  double hold_real=0;
  double hold_img=0;
  double hold_cos=0;
  double hold_sin=0;
  double hold=0;
  for(unsigned int i=0;i<size;i++)
  {
    hold=scalar*diagonal_matrix[i];
    hold_cos=cos(hold);
    hold_sin=sin(hold);
    hold_real=vector[2*i+0];
    hold_img=vector[2*i+1];
    vector[2*i+0]=hold_real*hold_cos + hold_sin*hold_img;
    vector[2*i+1]=(-hold_real*hold_sin) + hold_cos*hold_img;
  }
}

inline double statistical_distace(const double *state_0,const double *state_1,const unsigned int dim)
{
  double hold_0=0;
  double hold_1=0;
  double result=0;
  for(unsigned int i=0;i<dim;i++)
  {
    hold_0=state_0[2*i+0]*state_0[2*i+0]+state_0[2*i+1]*state_0[2*i+1];
    hold_1=state_1[2*i+0]*state_1[2*i+0]+state_1[2*i+1]*state_1[2*i+1];
    result+=abs(hold_0-hold_1);
  }
  result=0.5*result;
  return result;
}

//*************************************************** Physical functions ***************************************************//

inline double flux_freq_with_high_accuracy(double EC,double EJ1,double EJ2,double flux_t,double flux_t0)
{
  double EJ=sqrt(EJ1*EJ1+EJ2*EJ2+2*EJ1*EJ2*cos(1*(flux_t+flux_t0)));
  //double d=(EJ2-EJ1)/(EJ1+EJ2);
  //double EJ=(EJ1+EJ2)*sqrt(cos((flux_t+flux_t0))*cos((flux_t+flux_t0))+d*d*sin((flux_t+flux_t0))*sin((flux_t+flux_t0)));
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

inline double flux_freq_anharmonicity_with_high_accuracy(double EC,double EJ1,double EJ2,double flux_t,double flux_t0)
{

  double EJ=sqrt(EJ1*EJ1+EJ2*EJ2+2*EJ1*EJ2*cos(1*(flux_t+flux_t0)));
  //double d=(EJ2-EJ1)/(EJ1+EJ2);
  //double EJ=(EJ1+EJ2)*sqrt(cos((flux_t+flux_t0))*cos((flux_t+flux_t0))+d*d*sin((flux_t+flux_t0))*sin((flux_t+flux_t0)));
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

inline double fixed_frequency_transmon_effective_interaction_strength(double EC, double EJ)
{
  return sqrt(sqrt((EJ/(8*EC))));
}

inline double flux_tunable_transmon_effective_interaction_strength(double EC, double EJ1,double EJ2,double offset)
{
  return sqrt(sqrt(sqrt(EJ1*EJ1+EJ2*EJ2+2*EJ1*EJ2*cos(offset))/(8*EC)));
}

inline double effective_interaction_strength_factor_tunable_transmon(double EC, double EJ1,double EJ2,double flux)
{
  double EJ_time=sqrt(EJ1*EJ1+EJ2*EJ2+2*EJ1*EJ2*cos(flux));
  double result=sqrt(sqrt(EJ_time/(8*EC)));
  return result;
}

inline double effective_interaction_strength_factor_fixed_transmon(double EC, double EJ)
{
  double result=sqrt(sqrt(EJ/(8*EC)));
  return result;
}

#endif
