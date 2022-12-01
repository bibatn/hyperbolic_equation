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
