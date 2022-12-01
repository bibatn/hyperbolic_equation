CXX=g++
CXXFLAGS= -O3 -Wall

TARGET = hyperbolic_equation

all: $(TARGET)

hyperbolic_equation:
	mpic++ -O3 -Wall -o hyperbolic_equation -c hyperbolic_equation.cpp

build: hyperbolic_equation
	mpic++ hyperbolic_equation -O3 -std=c++11 -fopenmp main.cpp -o main

run:
	mpirun -np 4 ./main 128 1 out.txt

submit-polus-parallel:
	g++ -O3 -std=c++11 -fopenmp main_mpi.cpp -o task3 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm
	for N in 128 256 512 ; do \
		for p in 1 4 8 16 32 ; do \
			for i in {1..5} ; do \
				bsub -n $$p -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./main $$N 1 out\_$$p\_$$N\_1.txt ; \
				bsub -n $$p -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./main $$N out\_$$p\_$$N\_pi.txt ; \
			done \
		done \
	done
