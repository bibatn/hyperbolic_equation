CXX=g++
CXXFLAGS= -O3 -Wall

TARGET = hyperbolic_equation

all: $(TARGET)

build: hyperbolic_equation.h
	mpic++ -O3 -std=c++11 -fopenmp hyperbolic_equation.cpp -o hyperbolic_equation

run:
	mpirun -np 4 ./hyperbolic_equation 128 1 out.txt

submit-polus-parallel:
	for N in 128 256 512 ; do \
		for p in 1 4 8 16 32 ; do \
			for i in {1..5} ; do \
				bsub -n $$p -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./hyperbolic_equation $$N 1 out\_$$p\_$$N\_1.txt ; \
				bsub -n $$p -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./hyperbolic_equation $$N out\_$$p\_$$N\_pi.txt ; \
			done \
		done \
	done
