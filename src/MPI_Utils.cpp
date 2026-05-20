#include "MPI_Utils.h"

#include "Globals.h"
#include "Grid.h"

#include <algorithm>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace
{
int mpi_rank = 0;
int mpi_size = 1;
bool mpi_initialized_here = false;

int owned_begin()
{
    if (mpi_size <= 1)
        return 0;
    return (No_Physical_Cells * mpi_rank) / mpi_size;
}

int owned_end()
{
    if (mpi_size <= 1)
        return No_Physical_Cells;
    return (No_Physical_Cells * (mpi_rank + 1)) / mpi_size;
}

void flatten_owned_state(const vector<V_D> &source, const int components, vector<double> &buffer)
{
    const int begin = owned_begin();
    const int end = owned_end();
    buffer.assign(static_cast<size_t>(end - begin) * components, 0.0);

    for (int c = begin; c < end; ++c)
    {
        if (c >= static_cast<int>(source.size()))
            continue;
        for (int v = 0; v < components && v < static_cast<int>(source[c].size()); ++v)
            buffer[static_cast<size_t>(c - begin) * components + v] = source[c][v];
    }
}

void scatter_physical_state(vector<V_D> &target, const vector<double> &buffer, const int components)
{
    if (target.size() < static_cast<size_t>(No_Physical_Cells))
        return;

    for (int c = 0; c < No_Physical_Cells; ++c)
    {
        if (target[c].size() < static_cast<size_t>(components))
            target[c].resize(components, 0.0);
        for (int v = 0; v < components; ++v)
            target[c][v] = buffer[static_cast<size_t>(c) * components + v];
    }
}

#ifdef USE_MPI
void allgather_state(vector<V_D> &state, const int components)
{
    if (mpi_size <= 1 || No_Physical_Cells <= 0)
        return;

    vector<double> sendbuf;
    flatten_owned_state(state, components, sendbuf);

    vector<int> recvcounts(mpi_size, 0);
    vector<int> displs(mpi_size, 0);
    for (int r = 0; r < mpi_size; ++r)
    {
        const int begin = (No_Physical_Cells * r) / mpi_size;
        const int end = (No_Physical_Cells * (r + 1)) / mpi_size;
        recvcounts[r] = (end - begin) * components;
        displs[r] = begin * components;
    }

    vector<double> recvbuf(static_cast<size_t>(No_Physical_Cells) * components, 0.0);
    MPI_Allgatherv(sendbuf.data(), static_cast<int>(sendbuf.size()), MPI_DOUBLE,
                   recvbuf.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                   MPI_COMM_WORLD);

    scatter_physical_state(state, recvbuf, components);
}
#endif
} // namespace

void CFD_MPI_Initialize(int *argc, char ***argv)
{
#ifdef USE_MPI
    int already_initialized = 0;
    MPI_Initialized(&already_initialized);
    if (!already_initialized)
    {
        MPI_Init(argc, argv);
        mpi_initialized_here = true;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#else
    (void)argc;
    (void)argv;
#endif
}

void CFD_MPI_Finalize()
{
#ifdef USE_MPI
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (mpi_initialized_here && !finalized)
        MPI_Finalize();
#endif
}

bool CFD_MPI_Is_Enabled()
{
    return mpi_size > 1;
}

bool CFD_MPI_Is_Root()
{
    return mpi_rank == 0;
}

int CFD_MPI_Rank()
{
    return mpi_rank;
}

int CFD_MPI_Size()
{
    return mpi_size;
}

bool CFD_MPI_Owns_Cell(int cellIndex)
{
    if (cellIndex < 0 || cellIndex >= No_Physical_Cells)
        return false;
    return cellIndex >= owned_begin() && cellIndex < owned_end();
}

void CFD_MPI_Build_Local_Leaf_Cell_List(std::vector<int> &leafCells)
{
    Build_Leaf_Cell_List(leafCells);
    if (!CFD_MPI_Is_Enabled())
        return;

    leafCells.erase(std::remove_if(leafCells.begin(), leafCells.end(),
                                   [](const int cell) { return !CFD_MPI_Owns_Cell(cell); }),
                    leafCells.end());
}

double CFD_MPI_Global_Min(double value)
{
#ifdef USE_MPI
    if (mpi_size > 1)
    {
        double global_value = value;
        MPI_Allreduce(&value, &global_value, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        return global_value;
    }
#endif
    return value;
}

double CFD_MPI_Global_Max(double value)
{
#ifdef USE_MPI
    if (mpi_size > 1)
    {
        double global_value = value;
        MPI_Allreduce(&value, &global_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return global_value;
    }
#endif
    return value;
}

void CFD_MPI_Global_Sum4(double values[4])
{
#ifdef USE_MPI
    if (mpi_size > 1)
    {
        double global_values[4] = {0.0, 0.0, 0.0, 0.0};
        MPI_Allreduce(values, global_values, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (int i = 0; i < 4; ++i)
            values[i] = global_values[i];
    }
#else
    (void)values;
#endif
}

void CFD_MPI_Synchronize_Solution_State()
{
#ifdef USE_MPI
    if (mpi_size <= 1)
        return;

    allgather_state(U_Cells, NUM_FLUX_COMPONENTS);
    allgather_state(Primitive_Cells, 11);
#endif
}

void CFD_MPI_Barrier()
{
#ifdef USE_MPI
    if (mpi_size > 1)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
}
