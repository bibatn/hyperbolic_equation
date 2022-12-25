CXX=g++
CXXFLAGS = -O3 -Wall
OPT = -acc -ta=nvidia, -Minfo=accel
TARGET = hyperbolic_equation

all: $(TARGET)

build: hyperbolic_equation.h
	mpic++ -O3 -std=c++11 -fopenmp hyperbolic_equation.cpp -o hyperbolic_equation

build_pgc++: hyperbolic_equation.h
	pgc++ -std=c++11 $(OPT) hyperbolic_equation.cpp -o hyperbolic_equation

build_: hyperbolic_equation.h
	mpic++ -O3 -std=c++11 $(OPT) hyperbolic_equation.cpp -o hyperbolic_equation

build_p: hyperbolic_equation.h
	pgc++ -O3 -std=c++11 $(OPT) hyperbolic_equation.cpp -o hyperbolic_equation

run:
	mpirun -np 4 ./hyperbolic_equation 128 1 out.txt

run_polus:
	bsub -n 1 -gpu "num=2" -o o.txt -e err.txt mpiexec ./hyperbolic_equation 128 1 out.txt

submit-polus-parallel:
	for N in 512 ; do \
		for p in 1 4 16; do \
			for i in {1..2} ; do \
				bsub -n $$p -m polus-c3-ib -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./hyperbolic_equation $$N 1 out\_$$p\_$$N\_1.txt ; \
			done \
		done \
	done
