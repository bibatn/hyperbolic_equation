#include "hyperbolic_equation.h"
#include <chrono>
#include <iostream>
#include <cmath>
#include <mpi.h>
#include <fstream>
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

double Functions::AnalyticalSolution(double x, double y, double z, double t) const
{
    return sin(M_PI * x / g.L_x) * sin(M_PI * y / g.L_y) * sin(2 * M_PI * z / g.L_z) * cos(a * t);
}

double Functions::Phi(double x, double y, double z) const
{
    return AnalyticalSolution(x, y, z, 0);
}



int main(int argc, char** argv) {
    // MPI initialization
    int mp, np;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mp);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int timeSteps = 20;
    int K = 1000;
    double T = 1.0;

    // parameters parsing
    // N (x,y,z points), L (Lx=Ly=Lz=pi), out_file
    int N = atoi(argv[1]);
    double L = (argc == 4) ? strtod(argv[2], NULL) : M_PI;
    char* filename = (argc == 4) ? argv[3] : argv[2];

    double start_time = MPI_Wtime();

    // equation solving
    Grid grid = Grid(L, L, L, N, T, K);
    SolverMPI solver = SolverMPI(grid);
    double error = solver.Solve(timeSteps);

    // elapsed time computing
    double end_time = MPI_Wtime();
    double delta = end_time - start_time;
    double max_time;
    MPI_Reduce(&delta, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (mp == 0) {
        std::ofstream fout(filename, std::ios_base::app);
        fout << N << " " << np << " " << error << " " << max_time << std::endl;
        fout.close();
    }
    MPI_Finalize();
    return 0;
}




