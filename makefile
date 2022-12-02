CPPFLAGS	= -O2 -std=c++17 -pedantic -Wall -Wunknown-pragmas -Wextra -lstdc++fs # -fopenmp -lomp -Winline
INTELFLAGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread -diag-disable=remark
AUX				=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread

CPP=icpc
Adia=0#FALSE: Non-adiabatic
SQ_P=0 #Gauss envelope for harmonic charge drive
TQ_P=2 #Sine envelope for harmonic flux drive
MARCOS= -D SQ_PULSE=$(SQ_P) -D TQ_PULSE=$(TQ_P) -D ADIABATIC=$(Adia) #-D USE_OPENMP

time_evol: time_evol.cpp
	@$(CPP) $(MARCOS) time_evol.cpp -o main $(CPPFLAGS) $(INTELFLAGS);

trotter_error: trotter_error.cpp
	@$(CPP) $(MARCOS) trotter_error.cpp -o main $(CPPFLAGS) $(INTELFLAGS);

spectrum: spectrum.cpp
	@$(CPP) $(MARCOS) spectrum.cpp -o main $(CPPFLAGS) $(INTELFLAGS);
