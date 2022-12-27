#include <cmath>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <stdexcept>
#include "openacc.h"
struct Block
{
    Block(int x_min, int x_max, int y_min, int y_max, int z_min, int z_max) :
            x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), z_min(z_min), z_max(z_max)
    {
        x_size = x_max - x_min + 1;
        y_size = y_max - y_min + 1;
        z_size = z_max - z_min + 1;
        size = x_size * y_size * z_size;
    }

    void expand_Block()
    {
        x_min = x_min - 1;
        y_min = y_min - 1;
        z_min = z_min - 1;
        x_size = x_max - x_min + 1;
        y_size = y_max - y_min + 1;
        z_size = z_max - z_min + 1;
        size = x_size * y_size * z_size;
    }

    void narrow_Block()
    {
        x_min = x_min + 1;
        y_min = y_min + 1;
        z_min = z_min + 1;
        x_size = x_max - x_min + 1;
        y_size = y_max - y_min + 1;
        z_size = z_max - z_min + 1;
        size = x_size * y_size * z_size;
    }

    int x_min, y_min, z_min;
    int x_max, y_max, z_max;
    int x_size, y_size, z_size;
    int size;
};
struct Grid
{
    Grid(double L_x, double L_y, double L_z, int N, double T, int K)
    {
        this->L_x = L_x;
        this->L_y = L_y;
        this->L_z = L_z;
        this->N = N;

        this->T = T;
        this->K = K;

        this->h_x = L_x / N;
        this->h_y = L_y / N;
        this->h_z = L_z / N;

        this->tau = T / K;
    }

    double L_x, L_y, L_z;
    double h_x, h_y, h_z;
    int N;

    double T;
    double tau;
    int K;
};

#pragma acc routine
double analyticalSolution(double x, double y, double z, double t, double a, const Grid g) {
    return sin(M_PI * x / g.L_x) * sin(M_PI * y / g.L_y) * sin(2 * M_PI * z / g.L_z) * cos(a * t);
}

#pragma acc routine
double phi(double x, double y, double z, double a, const Grid g) {
    return analyticalSolution(x, y, z, 0, a, g);
}

#pragma acc routine
inline int index(int i, int j, int k, const Block b) {
    // get the linear index inside the array of the given grid block
    return (i - b.x_min) * b.y_size * b.z_size + (j - b.y_min) * b.z_size + (k - b.z_min);
}

enum Axis
{
    X, Y, Z
};


class Functions
{
    double a; // a_t in the analytical solution
    Grid g; // grid with the main data related to the solving
public:
    explicit Functions(Grid g): g(g)
    {
        a = M_PI * sqrt(1.0 / (g.L_x * g.L_x) + 1.0 / (g.L_y * g.L_y) + 4.0 / (g.L_z * g.L_z));
    }

    double AnalyticalSolution(double x, double y, double z, double t) const;

    double Phi(double x, double y, double z) const;

};


class Index
{
    Grid g;
public:
    explicit Index(Grid g) : g(g) {}

    inline int operator()(int i, int j, int k) const
    {
        // 3d-array of points is flattened, get the linear index
        return (i * (g.N + 1) + j) * (g.N + 1) + k;
    }

    inline int operator()(int i, int j, int k, Block b) const
    {
        // get the linear index inside the array of the given grid block
        return (i - b.x_min) * b.y_size * b.z_size + (j - b.y_min) * b.z_size + (k - b.z_min);
    }

};

void split(int x_min, int x_max, int y_min, int y_max, int z_min, int z_max, int size, Axis axis, std::vector<Block> &blocks)
{
    // split the grid recursively into blocks
    if (size == 1)
    {
        blocks.emplace_back(x_min, x_max, y_min, y_max, z_min, z_max);
        return;
    }

    if (size % 2 == 1)
    { // if size is odd we make it even
        if (axis == X)
        {
            int x = x_min + (x_max - x_min) / size;
            blocks.emplace_back(x_min, x, y_min, y_max, z_min, z_max);
            x_min = x + 1;
            axis = Y;
        }
        else if (axis == Y)
        {
            int y = y_min + (y_max - y_min) / size;
            blocks.emplace_back(x_min, x_max, y_min, y, z_min, z_max);
            y_min = y + 1;
            axis = Z;
        }
        else
        {
            int z = z_min + (z_max - z_min) / size;
            blocks.emplace_back(x_min, x_max, y_min, y_max, z_min, z);
            z_min = z + 1;
            axis = X;
        }

        size--;
    }

    // now the size is even
    if (axis == X)
    {
        int x = (x_min + x_max) / 2;
        split(x_min, x, y_min, y_max, z_min, z_max, size / 2, Y, blocks);
        split(x + 1, x_max, y_min, y_max, z_min, z_max, size / 2, Y, blocks);
    }
    else if (axis == Y)
    {
        int y = (y_min + y_max) / 2;
        split(x_min, x_max, y_min, y, z_min, z_max, size / 2, Z, blocks);
        split(x_min, x_max, y + 1, y_max, z_min, z_max, size / 2, Z, blocks);
    }
    else
    {
        int z = (z_min + z_max) / 2;
        split(x_min, x_max, y_min, y_max, z_min, z, size / 2, X, blocks);
        split(x_min, x_max, y_min, y_max, z + 1, z_max, size / 2, X, blocks);
    }
}

class SolverMPI
{
    Functions f;
    Grid g;
    Index ind; // for getting flattened indexes in the 3-d array
    double * u0;//add destructor
    double * u1;
    double * u2;
    std::vector< std::pair<int, Block> > blocksToSend;
    std::vector< std::pair<int, Block> > blocksToReceive;
    int * offset_vector;
    double * dataToReceive;
    int proc_rank, proc_size;

public:

    explicit SolverMPI(Grid g) : g(g), f(g), ind(g)
    {
        proc_rank = 0; proc_size = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
#pragma acc enter data copyin(this)
    }

    std::vector<double> GetSendData0(const Block block, const Block otherBlock) const
    {
        std::vector<double> dataToSend(otherBlock.size);

        for (int i = otherBlock.x_min; i <= otherBlock.x_max; i++)
            for (int j = otherBlock.y_min; j <= otherBlock.y_max; j++)
                for (int k = otherBlock.z_min; k <= otherBlock.z_max; k++)
                    dataToSend[ind(i, j, k, otherBlock)] = u0[ind(i, j, k, block)];

        return dataToSend;
    }
    
    std::vector<double> GetSendData1(const Block block, const Block otherBlock) const
    {
        std::vector<double> dataToSend(otherBlock.size);

        for (int i = otherBlock.x_min; i <= otherBlock.x_max; i++)
            for (int j = otherBlock.y_min; j <= otherBlock.y_max; j++)
                for (int k = otherBlock.z_min; k <= otherBlock.z_max; k++)
                    dataToSend[ind(i, j, k, otherBlock)] = u1[ind(i, j, k, block)];

        return dataToSend;
    }

#pragma acc routine
    void Exchange(int uInd, const Block block)
    {
    //        std::vector<double> dataToReceive(blocksToReceive.size());
    //        std::cout << "I'm here! " << std::endl;
        for (int i = 0; i < blocksToReceive.size(); i++) {
            std::vector<double> dataToSend;
            if (uInd == 0)
                dataToSend = GetSendData0(block, blocksToSend[i].second);
            else if (uInd == 1)
                dataToSend = GetSendData1(block, blocksToSend[i].second);
            MPI_Sendrecv(dataToSend.data(), blocksToSend[i].second.size, MPI_DOUBLE, blocksToSend[i].first, 777, &dataToReceive[offset_vector[i]], blocksToReceive[i].second.size, MPI_DOUBLE, blocksToReceive[i].first, 777, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        }
    }


#pragma acc routine
    double FindU(int i, int j, int k, const Block b) const {

        if (b.x_min <= i and i <= b.x_max and b.y_min <= j and j <= b.y_max and b.z_min <= k and k <= b.z_max) {
            return u1[ind(i, j, k, b)];
        }

        for (int r_i = 0; r_i < blocksToReceive.size(); r_i++) {
            Block otherB = blocksToReceive[r_i].second;

            if (i < otherB.x_min or i > otherB.x_max or
                j < otherB.y_min or j > otherB.y_max or
                k < otherB.z_min or k > otherB.z_max)
                continue;

            return dataToReceive[offset_vector[r_i] + ind(i, j, k, otherB)];
        }
        return 1;
    }
#pragma acc routine
    double FindU0(int i, int j, int k, const Block b) const {

        if (b.x_min <= i and i <= b.x_max and b.y_min <= j and j <= b.y_max and b.z_min <= k and k <= b.z_max) {
            return u0[ind(i, j, k, b)];
        }

        for (int r_i = 0; r_i < blocksToReceive.size(); r_i++) {
            Block otherB = blocksToReceive[r_i].second;

            if (i < otherB.x_min or i > otherB.x_max or
                j < otherB.y_min or j > otherB.y_max or
                k < otherB.z_min or k > otherB.z_max)
                continue;

            return dataToReceive[offset_vector[r_i] + ind(i, j, k, otherB)];
        }
        return 1;
    }

    double LaplaceOperator(int i, int j, int k, const Block b) const
    {
        double dx = (FindU(i, j - 1, k, b) - 2 * u1[ind(i, j, k, b)] + FindU(i, j + 1, k, b)) / (g.h_y * g.h_y);
        double dy = (FindU(i - 1, j, k, b) - 2 * u1[ind(i, j, k, b)] + FindU(i + 1, j, k, b)) / (g.h_x * g.h_x);
        double dz = (FindU(i, j, k - 1, b) - 2 * u1[ind(i, j, k, b)] + FindU(i, j, k + 1, b)) / (g.h_z * g.h_z);
        return dx + dy + dz;
    }

    double LaplaceOperator0(int i, int j, int k, const Block b) const
    {
        double dx = (FindU0(i, j - 1, k, b) - 2 * u0[ind(i, j, k, b)] + FindU0(i, j + 1, k, b)) / (g.h_y * g.h_y);
        double dy = (FindU0(i - 1, j, k, b) - 2 * u0[ind(i, j, k, b)] + FindU0(i + 1, j, k, b)) / (g.h_x * g.h_x);
        double dz = (FindU0(i, j, k - 1, b) - 2 * u0[ind(i, j, k, b)] + FindU0(i, j, k + 1, b)) / (g.h_z * g.h_z);
        return dx + dy + dz;
    }

    double ComputeLayerError(double t, const Block b) const
    {
        double errorLocal = 0;
        // maximum difference between values of u analytical and u computed
        for (int i = b.x_min; i <= b.x_max; i++)
            for (int j = b.y_min; j <= b.y_max; j++)
                for (int k = b.z_min; k <= b.z_max; k++)
                    errorLocal = std::max(errorLocal, fabs(u1[ind(i, j, k, b)] -
                                                           f.AnalyticalSolution(i * g.h_x, j * g.h_y, k * g.h_z, t)));
        double error = 0;
        MPI_Reduce(&errorLocal, &error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        return error;
    }

    void FillBoundaryValues0(double t, const Block b)
    {
        // Variant 3 -> first kind for x, periodic for y, first kind for z
        if (b.x_min == 0) {
            for (int i = b.y_min; i <= b.y_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u0[ind(b.x_min, i, j, b)] = 0;
        }

        if (b.x_max == g.N) {
            for (int i = b.y_min; i <= b.y_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u0[ind(b.x_max, i, j, b)] = 0;
        }

        if (b.y_min == 0) {
            for(int i = b.x_min; i<=b.x_max; i++)
                for(int j = b.z_min; j<=b.z_max; j++)
                    u0[ind(i, b.y_min, j, b)] = 0;
        }

        if(b.y_max == g.N){
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u0[ind(i, b.y_max, j, b)] = 0;
        }

        if (b.z_min == 0) {
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.y_min; j <= b.y_max; j++)
                    u0[ind(i, j, b.z_min, b)] = f.AnalyticalSolution(i * g.h_x, j * g.h_y, 0, t);
        }

        if (b.z_max == g.N) {
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.y_min; j <= b.y_max; j++)
                    u0[ind(i, j, b.z_max, b)] = f.AnalyticalSolution(i * g.h_x, j * g.h_y, g.L_z, t);
        }
    }

    void FillBoundaryValues1(double t, const Block b)
    {
        // Variant 3 -> first kind for x, periodic for y, first kind for z
        if (b.x_min == 0) {
            for (int i = b.y_min; i <= b.y_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u1[ind(b.x_min, i, j, b)] = 0;
        }

        if (b.x_max == g.N) {
            for (int i = b.y_min; i <= b.y_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u1[ind(b.x_max, i, j, b)] = 0;
        }

        if (b.y_min == 0) {
            for(int i = b.x_min; i<=b.x_max; i++)
                for(int j = b.z_min; j<=b.z_max; j++)
                    u1[ind(i, b.y_min, j, b)] = 0;
        }

        if(b.y_max == g.N){
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u1[ind(i, b.y_max, j, b)] = 0;
        }

        if (b.z_min == 0) {
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.y_min; j <= b.y_max; j++)
                    u1[ind(i, j, b.z_min, b)] = f.AnalyticalSolution(i * g.h_x, j * g.h_y, 0, t);
        }

        if (b.z_max == g.N) {
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.y_min; j <= b.y_max; j++)
                    u1[ind(i, j, b.z_max, b)] = f.AnalyticalSolution(i * g.h_x, j * g.h_y, g.L_z, t);
        }
    }

    void FillBoundaryValues2(double t, const Block b)
    {
        // Variant 3 -> first kind for x, periodic for y, first kind for z
        if (b.x_min == 0) {
            for (int i = b.y_min; i <= b.y_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u2[ind(b.x_min, i, j, b)] = 0;
        }

        if (b.x_max == g.N) {
            for (int i = b.y_min; i <= b.y_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u2[ind(b.x_max, i, j, b)] = 0;
        }

        if (b.y_min == 0) {
            for(int i = b.x_min; i<=b.x_max; i++)
                for(int j = b.z_min; j<=b.z_max; j++)
                    u2[ind(i, b.y_min, j, b)] = 0;
        }

        if(b.y_max == g.N){
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u2[ind(i, b.y_max, j, b)] = 0;
        }

        if (b.z_min == 0) {
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.y_min; j <= b.y_max; j++)
                    u2[ind(i, j, b.z_min, b)] = f.AnalyticalSolution(i * g.h_x, j * g.h_y, 0, t);
        }

        if (b.z_max == g.N) {
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.y_min; j <= b.y_max; j++)
                    u2[ind(i, j, b.z_max, b)] = f.AnalyticalSolution(i * g.h_x, j * g.h_y, g.L_z, t);
        }
    }

    void InitValues(const Block b)
    {
        // boundary (i = 0,N or j = 0,N or k = 0,N)
        FillBoundaryValues0(0, b);
        FillBoundaryValues1(g.tau, b);

        // compute the boundaries of the current block
        int x1 = std::max(b.x_min, 1); int x2 = std::min(b.x_max, g.N - 1);
        int y1 = std::max(b.y_min, 1); int y2 = std::min(b.y_max, g.N - 1);
        int z1 = std::max(b.z_min, 1); int z2 = std::min(b.z_max, g.N - 1);

        // initial values for inner points in u_0

#pragma acc update device(u0[0:b.size])
#pragma acc kernels loop independent
        for (int i = x1; i <= x2; i++)
            for (int j = y1; j <= y2; j++)
                for (int k = z1; k <= z2; k++)
                    u0[index(i, j, k, b)] = f.Phi(i * g.h_x, j * g.h_y, k * g.h_z);


        Exchange(0, b);
        // initial values for inner points in u_1
#pragma acc kernels loop independent
        for (int i = x1; i <= x2; i++)
            for (int j = y1; j <= y2; j++)
                for (int k = z1; k <= z2; k++)
                    u1[ind(i, j, k, b)] = u0[ind(i, j, k, b)] + g.tau * g.tau / 2 * LaplaceOperator0(i, j, k, b);
    }

    void GetNextU(int step, const Block b)
    {
        // compute the boundaries of the current block
        int x1 = std::max(b.x_min, 1); int x2 = std::min(b.x_max, g.N - 1);
        int y1 = std::max(b.y_min, 1); int y2 = std::min(b.y_max, g.N - 1);
        int z1 = std::max(b.z_min, 1); int z2 = std::min(b.z_max, g.N - 1);

        Exchange(1, b);
        // calculate u_n+1 inside the area
#pragma acc kernels loop independent
        for (int i = x1; i <= x2; i++)
            for (int j = y1; j <= y2; j++)
                for (int k = z1; k <= z2; k++)
                    u2[index(i, j, k, b)] = 2 * u1[index(i, j, k, b)] -
                                                   u0[index(i, j, k, b)] +
                                                   g.tau * g.tau * LaplaceOperator( i, j, k, b);
        FillBoundaryValues2(step * g.tau, b);
    }

    bool IsInside(int xmin1, int xmax1, int ymin1, int ymax1, int xmin2, int xmax2, int ymin2, int ymax2) const
    {
        return xmin2 <= xmin1 && xmax1 <= xmax2 && ymin2 <= ymin1 && ymax1 <= ymax2;
    }

    void GetNeighbours(const std::vector<Block> &blocks)
    {
        Block block = blocks[proc_rank];

        for (int i = 0; i < proc_size; i++)
        {
            if (i == proc_rank)
                continue;

            Block otherBlock = blocks[i];
            if (block.x_min == otherBlock.x_max + 1 or otherBlock.x_min == block.x_max + 1)
            {

                int xSend = block.x_min == otherBlock.x_max + 1 ? block.x_min : block.x_max;
                int xRecv = otherBlock.x_min == block.x_max + 1 ? otherBlock.x_min : otherBlock.x_max;
                int y_min, y_max, z_min, z_max;
                if (IsInside(block.y_min, block.y_max, block.z_min, block.z_max,
                             otherBlock.y_min, otherBlock.y_max, otherBlock.z_min, otherBlock.z_max)) {
                    y_min = block.y_min; y_max = block.y_max; z_min = block.z_min; z_max = block.z_max;
                } else if (IsInside(otherBlock.y_min, otherBlock.y_max, otherBlock.z_min, otherBlock.z_max,
                                    block.y_min, block.y_max, block.z_min, block.z_max)) {
                    y_min = otherBlock.y_min; y_max = otherBlock.y_max; z_min = otherBlock.z_min; z_max = otherBlock.z_max;
                } else
                    continue;
                // add block as a rectangle (it's a border between two processes)
                blocksToSend.emplace_back(i, Block(xSend, xSend, y_min, y_max, z_min, z_max));
                blocksToReceive.emplace_back(i, Block(xRecv, xRecv, y_min, y_max, z_min, z_max));
                continue;
            }
            if (block.y_min == otherBlock.y_max + 1 or otherBlock.y_min == block.y_max + 1)
            {
                int ySend = block.y_min == otherBlock.y_max + 1 ? block.y_min : block.y_max;
                int yRecv = otherBlock.y_min == block.y_max + 1 ? otherBlock.y_min : otherBlock.y_max;
                int x_min, x_max, z_min, z_max;
                if (IsInside(block.x_min, block.x_max, block.z_min, block.z_max,
                             otherBlock.x_min, otherBlock.x_max, otherBlock.z_min, otherBlock.z_max)) {
                    x_min = block.x_min; x_max = block.x_max; z_min = block.z_min; z_max = block.z_max;
                } else if (IsInside(otherBlock.x_min, otherBlock.x_max, otherBlock.z_min, otherBlock.z_max,
                                    block.x_min, block.x_max, block.z_min, block.z_max)) {
                    x_min = otherBlock.x_min; x_max = otherBlock.x_max; z_min = otherBlock.z_min; z_max = otherBlock.z_max;
                } else
                    continue;
                blocksToSend.emplace_back(i, Block(x_min, x_max, ySend, ySend, z_min, z_max));
                blocksToReceive.emplace_back(i, Block(x_min, x_max, yRecv, yRecv, z_min, z_max));
                continue;
            }
            if (block.z_min == otherBlock.z_max + 1 or otherBlock.z_min == block.z_max + 1)
            {
                int zSend = block.z_min == otherBlock.z_max + 1 ? block.z_min : block.z_max;
                int zRecv = otherBlock.z_min == block.z_max + 1 ? otherBlock.z_min : otherBlock.z_max;
                int x_min, x_max, y_min, y_max;
                if (IsInside(block.x_min, block.x_max, block.y_min, block.y_max,
                             otherBlock.x_min, otherBlock.x_max, otherBlock.y_min, otherBlock.y_max)) {
                    x_min = block.x_min; x_max = block.x_max; y_min = block.y_min; y_max = block.y_max;
                } else if (IsInside(otherBlock.x_min, otherBlock.x_max, otherBlock.y_min, otherBlock.y_max,
                                    block.x_min, block.x_max, block.y_min, block.y_max)) {
                    x_min = otherBlock.x_min; x_max = otherBlock.x_max; y_min = otherBlock.y_min; y_max = otherBlock.y_max;
                } else
                    continue;
                blocksToSend.emplace_back(i, Block(x_min, x_max, y_min, y_max, zSend, zSend));
                blocksToReceive.emplace_back(i, Block(x_min, x_max, y_min, y_max, zRecv, zRecv));
                continue;
            }
        }
    }

    double Solve(int steps)
    {
        // split grid between processes
        std::vector<Block> blocks;
        split(0, g.N, 0, g.N, 0, g.N, proc_size, X, blocks);
        Block block = blocks[proc_rank];
//        block.expand_Block();

        // allocate space for u
        u0 = new double[block.size];
        u1 = new double[block.size];
        u2 = new double[block.size];
#pragma acc enter data create(u0[0:block.size])
#pragma acc enter data create(u1[0:block.size])
#pragma acc enter data create(u2[0:block.size])
//        block.narrow_Block();
        // fill blocksToSend and blocksToReceive vectors
        GetNeighbours(blocks);

        offset_vector = new int [blocksToReceive.size()];
        offset_vector[0] = 0;
        int data_size = blocksToReceive[0].second.size;
        for (int i = 1; i < blocksToReceive.size(); i++)
        {
            offset_vector[i] = offset_vector[i-1] + blocksToReceive[i-1].second.size;
            data_size = data_size + blocksToReceive[i].second.size;
        }
        dataToReceive = new double [data_size];
#pragma acc enter data create(offset_vector[0:blocksToReceive.size()])
#pragma acc enter data create(dataToReceive[0:data_size])
        // init u_0 and u_1
        InitValues(block);

        // calculate the next time layers for u
        for (int step = 2; step <= steps; step++)
        {
            GetNextU(step, block);
            std::swap(u0,u2);
            std::swap(u0, u1);
        }
        double layerError = ComputeLayerError(steps * g.tau, block);
#pragma acc wait
#pragma acc exit data delete(u0)
#pragma acc exit data delete(u1)
#pragma acc exit data delete(u2)
#pragma acc exit data delete(offset_vector)
#pragma acc exit data delete(dataToReceive)
#pragma acc exit data delete(this)

        return layerError;
    }

};

