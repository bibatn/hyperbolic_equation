#include <cmath>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <stdexcept>

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
    std::vector< std::vector<double> > u;
    std::vector< std::pair<int, Block> > blocksToSend;
    std::vector< std::pair<int, Block> > blocksToReceive;
    int proc_rank, proc_size;

public:

    explicit SolverMPI(Grid g) : g(g), f(g), ind(g)
    {
        proc_rank = 0; proc_size = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
    }

    std::vector<double> GetSendData(int uInd, const Block block, const Block otherBlock) const
    {
        std::vector<double> dataToSend(otherBlock.size);

#pragma omp parallel for collapse(3)
        for (int i = otherBlock.x_min; i <= otherBlock.x_max; i++)
            for (int j = otherBlock.y_min; j <= otherBlock.y_max; j++)
                for (int k = otherBlock.z_min; k <= otherBlock.z_max; k++)
                    dataToSend[ind(i, j, k, otherBlock)] = u[uInd][ind(i, j, k, block)];

        return dataToSend;
    }

    std::vector< std::vector<double> > Exchange(int uInd, const Block block) const
    {
        std::vector< std::vector<double> > dataToReceive(blocksToReceive.size());
        std::vector<MPI_Request> requests(2);
        std::vector<MPI_Status> statuses(2);

        for (int i = 0; i < blocksToReceive.size(); i++) {
            std::vector<double> dataToSend = GetSendData(uInd, block, blocksToSend[i].second);
            dataToReceive[i] = std::vector<double>(blocksToReceive[i].second.size);
            MPI_Isend(dataToSend.data(), blocksToSend[i].second.size, MPI_DOUBLE, blocksToSend[i].first, 0, MPI_COMM_WORLD, &requests[0]);
            MPI_Irecv(dataToReceive[i].data(), blocksToReceive[i].second.size, MPI_DOUBLE, blocksToReceive[i].first, 0, MPI_COMM_WORLD, &requests[1]);
            MPI_Waitall(2, requests.data(), statuses.data());
        }
        return dataToReceive;
    }

    double FindU(int uInd, int i, int j, int k, const Block b, const std::vector< std::vector<double> > &recieved) const {

        if (b.x_min <= i and i <= b.x_max and b.y_min <= j and j <= b.y_max and b.z_min <= k and k <= b.z_max) {
            return u[uInd][ind(i, j, k, b)];
        }

        for (int r_i = 0; r_i < blocksToReceive.size(); r_i++) {
            Block otherB = blocksToReceive[r_i].second;

            if (i < otherB.x_min or i > otherB.x_max or
                j < otherB.y_min or j > otherB.y_max or
                k < otherB.z_min or k > otherB.z_max)
                continue;

            return recieved[r_i][ind(i, j, k, otherB)];
        }
        throw std::runtime_error("u value non found");
    }

    double LaplaceOperator(int uInd, int i, int j, int k, const Block b, const std::vector< std::vector<double> > &recieved) const
    {
        double dx = (FindU(uInd, i, j - 1, k, b, recieved) - 2 * u[uInd][ind(i, j, k, b)] + FindU(uInd, i, j + 1, k, b, recieved)) / (g.h_y * g.h_y);
        double dy = (FindU(uInd, i - 1, j, k, b, recieved) - 2 * u[uInd][ind(i, j, k, b)] + FindU(uInd, i + 1, j, k, b, recieved)) / (g.h_x * g.h_x);
        double dz = (FindU(uInd,i, j, k - 1, b, recieved) - 2 * u[uInd][ind(i, j, k, b)] + FindU(uInd, i, j, k + 1, b, recieved)) / (g.h_z * g.h_z);
        return dx + dy + dz;
    }

    double ComputeLayerError(int uInd, double t, const Block b) const
    {
        double errorLocal = 0;
        // maximum difference between values of u analytical and u computed
#pragma omp parallel for collapse(3) reduction(max: errorLocal)
        for (int i = b.x_min; i <= b.x_max; i++)
            for (int j = b.y_min; j <= b.y_max; j++)
                for (int k = b.z_min; k <= b.z_max; k++)
                    errorLocal = std::max(errorLocal, fabs(u[uInd][ind(i, j, k, b)] -
                                                           f.AnalyticalSolution(i * g.h_x, j * g.h_y, k * g.h_z, t)));
        double error = 0;
        MPI_Reduce(&errorLocal, &error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        return error;
    }

    void FillBoundaryValues(int uInd, double t, const Block b)
    {
        // Variant 3 -> first kind for x, periodic for y, first kind for z
        if (b.x_min == 0) {
#pragma omp parallel for collapse(2)
            for (int i = b.y_min; i <= b.y_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u[uInd][ind(b.x_min, i, j, b)] = 0;
        }

        if (b.x_max == g.N) {
#pragma omp parallel for collapse(2)
            for (int i = b.y_min; i <= b.y_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u[uInd][ind(b.x_max, i, j, b)] = 0;
        }

        if (b.y_min == 0) {
#pragma omp parallel for collapse(2)
            for(int i = b.x_min; i<=b.x_max; i++)
                for(int j = b.z_min; j<=b.z_max; j++)
                    u[uInd][ind(i, b.y_min, j, b)] = 0;
        }

        if(b.y_max == g.N){
#pragma omp parallel for collapse(2)
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u[uInd][ind(i, b.y_max, j, b)] = 0;
        }

        if (b.z_min == 0) {
#pragma omp parallel for collapse(2)
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.y_min; j <= b.y_max; j++)
                    u[uInd][ind(i, j, b.z_min, b)] = f.AnalyticalSolution(i * g.h_x, j * g.h_y, 0, t);
        }

        if (b.z_max == g.N) {
#pragma omp parallel for collapse(2)
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.y_min; j <= b.y_max; j++)
                    u[uInd][ind(i, j, b.z_max, b)] = f.AnalyticalSolution(i * g.h_x, j * g.h_y, g.L_z, t);
        }
    }

    void InitValues(const Block b)
    {
        // boundary (i = 0,N or j = 0,N or k = 0,N)
        FillBoundaryValues(0, 0, b);
        FillBoundaryValues(1, g.tau, b);

        // compute the boundaries of the current block
        int x1 = std::max(b.x_min, 1); int x2 = std::min(b.x_max, g.N - 1);
        int y1 = std::max(b.y_min, 1); int y2 = std::min(b.y_max, g.N - 1);
        int z1 = std::max(b.z_min, 1); int z2 = std::min(b.z_max, g.N - 1);

        // initial values for inner points in u_0
#pragma omp parallel for collapse(3)
        for (int i = x1; i <= x2; i++)
            for (int j = y1; j <= y2; j++)
                for (int k = z1; k <= z2; k++)
                    u[0][ind(i, j, k, b)] = f.Phi(i * g.h_x, j * g.h_y, k * g.h_z);

        std::vector< std::vector<double> > recieved = Exchange(0, b);
        // initial values for inner points in u_1
#pragma omp parallel for collapse(3)
        for (int i = x1; i <= x2; i++)
            for (int j = y1; j <= y2; j++)
                for (int k = z1; k <= z2; k++)
                    u[1][ind(i, j, k, b)] = u[0][ind(i, j, k, b)] + g.tau * g.tau / 2 * LaplaceOperator(0, i, j, k, b, recieved);
    }

    void GetNextU(int step, const Block b)
    {
        // compute the boundaries of the current block
        int x1 = std::max(b.x_min, 1); int x2 = std::min(b.x_max, g.N - 1);
        int y1 = std::max(b.y_min, 1); int y2 = std::min(b.y_max, g.N - 1);
        int z1 = std::max(b.z_min, 1); int z2 = std::min(b.z_max, g.N - 1);

        std::vector< std::vector<double> > received = Exchange((step + 2) % 3, b);
        // calculate u_n+1 inside the area
#pragma omp parallel for collapse(3)
        for (int i = x1; i <= x2; i++)
            for (int j = y1; j <= y2; j++)
                for (int k = z1; k <= z2; k++)
                    u[step % 3][ind(i, j, k, b)] = 2 * u[(step + 2) % 3][ind(i, j, k, b)] -
                                                   u[(step + 1) % 3][ind(i, j, k, b)] +
                                                   g.tau * g.tau * LaplaceOperator((step + 2) % 3, i, j, k, b, received);
        FillBoundaryValues(step % 3, step * g.tau, b);
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
        u.resize(3);
        for (int i = 0; i < 3; i++)
            u[i].resize(block.size);
//        block.narrow_Block();
        // fill blocksToSend and blocksToReceive vectors
        GetNeighbours(blocks);

        // init u_0 and u_1
        InitValues(block);

        // calculate the next time layers for u
        for (int step = 2; step <= steps; step++)
        {
            GetNextU(step, block);
        }

        return ComputeLayerError(steps % 3, steps * g.tau, block);
    }

};

