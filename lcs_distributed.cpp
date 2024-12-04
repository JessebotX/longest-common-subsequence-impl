#include "core/utils.h"
#include <mpi.h>
#include <stdio.h>

#define LCS_VERSION       "distributed"
#define DEFAULT_THREADS   1
#define DEFAULT_STRING_X  "abcd"
#define DEFAULT_STRING_Y  "acbad"

struct ThreadData {
    int id;
    double timeTaken;

    ThreadData(int _id, double _timeTaken)
        : id(_id), timeTaken(_timeTaken) {}
};

using namespace std;

void getCellsToCalculate(vector<int>& vec, vector<vector<int>>& dp,
                         int rows, int cols,
                         int thread_id, uint num_threads,
                         int diag) {

    // Find start position of the diagonal
    // (all threads compute the same number here)
    int start_row = std::max(0, diag - (cols - 1));
    int start_col = std::min(diag, cols - 1);
    
    // number of elements in this diagonal 
    // (all threads calculate the same number here)
    int elements_in_diagonal = std::min(start_col + 1, rows - start_row);

    // number of threads that get an extra cell to calculate 
    // (all threads calculate same num here)
    int z = elements_in_diagonal % num_threads;

    // min number of cells each thread gets 
    // (all threads calculate same number here)
    int base_count = elements_in_diagonal / num_threads;

    // the starting index along the current diagonal that the current thread needs to start at
    // this is the only number that is uniquely computed per thread
    int start_idx = std::min(thread_id, z) * (base_count + 1) + std::max(0, thread_id - z) * (base_count);

    // put the results in the vector
    vec[0] = start_row + start_idx + 1;
    vec[1] = start_col - start_idx + 1;
    vec[2] = base_count + (thread_id < (z) ? 1 : 0);
}

void findLCS(const string& X, const string& Y, ThreadData& processData, int world_rank, int world_size) {
    timer t1;
    t1.start();
    // --------------------------------------------------------------
    int n = X.length();
    int m = Y.length();
    uint num_diagonals = n + m - 1;

    vector<vector<int>> dp(n + 1, vector<int>(m + 1));

    vector<int> cellsToCalculate(3);
    int x, y, z;
    for (int i = 0; i < num_diagonals; i++) {
        getCellsToCalculate(cellsToCalculate, dp, n, m, world_rank, world_size, i);
        x = cellsToCalculate[0];
        y = cellsToCalculate[1];
        z = cellsToCalculate[2];

        for (int j = 0; j < z; j++) {
            if (X[x-1] == Y[y - 1]) {
                dp[x][y] = dp[x - 1][y - 1] + 1;
            } else {
                dp[x][y] = max(dp[x - 1][y], dp[x][y - 1]);
            }
            x += 1;
            y -= 1;
            //printf("%d\n", dp[x][y]);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    unordered_set<string> lcsSet;

    // --------------------------------------------------------------
    processData.timeTaken = t1.stop();
    return;
}

int main(int argc, char** argv) {
    // === BEGIN cli parsing ===
    cxxopts::Options options("lcs", "Find longest common subsequence.");
    options.add_options(
        "",
        { 
            { "t", "Number of threads", cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS) },
            { "x", "1st sequence", cxxopts::value<string>()->default_value(DEFAULT_STRING_X) },
            { "y", "2nd sequence", cxxopts::value<string>()->default_value(DEFAULT_STRING_Y) },
        }
    );
    auto cli = options.parse(argc, argv);

    // -t ignored in distributed
    uint nThreads = cli["t"].as<uint>();
    if (LCS_VERSION == "serial" && nThreads != 1) { // error when serial version > 1 thread
        printf("Error : Number of threads for serial version must be 1.\n");
        return 1;
    }

    string X = cli["x"].as<string>(); // sequence 1
    string Y = cli["y"].as<string>(); // sequence 2
    // === END cli parsing ===

    // === BEGIN mpi ===
    MPI_Init(NULL, NULL);

    // Get number of processses and rank of process
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) { // Print preface (from root)
        printf("LCS Version : %s\n", LCS_VERSION);
        printf("Number of threads : %u\n", world_size);
        printf("Sequence X : %s\n", X.c_str());
        printf("Sequence Y : %s\n", Y.c_str());
        printf("Finding longest common subsequence...\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // ------------------------------------------------------------
    timer mainTimer;
    mainTimer.start();

    ThreadData processData(world_rank, 0.0);

    findLCS(X, Y, processData, world_rank, world_size);

    // ------------------------------------------------------------
    double timeTaken = mainTimer.stop();

    MPI_Barrier(MPI_COMM_WORLD);
    // ------ BEGIN print ------

    // Stats per process
    if (world_rank == 0) {
        printf("id, time_taken\n");        
    }
    MPI_Barrier(MPI_COMM_WORLD);

    printf("%d, %.4f\n", world_rank, processData.timeTaken);

    // Total time taken
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0) {
        printf("Time taken (in seconds) : %.4f\n", timeTaken);
    }
    // ------ END print ------

    // === END mpi ===
    MPI_Finalize();

    return 0;
}
