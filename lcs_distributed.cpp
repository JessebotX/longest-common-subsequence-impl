#include "core/utils.h"
#include <mpi.h>
#include <stdio.h>

#define LCS_VERSION       "distributed"
#define DEFAULT_PROCESSES   1
#define DEFAULT_STRING_X  "abcd"
#define DEFAULT_STRING_Y  "acbad"

struct ProcessData {
    int id;
    double timeTaken;

    ProcessData(int _id, double _timeTaken)
        : id(_id), timeTaken(_timeTaken) {}
};

struct CellsToCalculate {
    int startingRow;
    int startingColumn;
    int numElements;       // number of elements to process per process (minPerProcess || minPerProcess + 1)
    int minPerProcess;     // min # of elements each process gets
    int elementsInDiagonal;

    CellsToCalculate()
        : startingRow(0), startingColumn(0), 
        numElements(0), minPerProcess(0), elementsInDiagonal(0)
    {}
}

using namespace std;

void getCellsToCalculate(CellsToCalculate& cellsToCalculate,
                         int rows, int cols,
                         int process_id, uint num_processes,
                         int diag) {

    // Find start position of the diagonal
    // (all processes compute the same number here)
    int start_row = std::max(0, diag - (cols - 1));
    int start_col = std::min(diag, cols - 1);
    
    // number of elements in this diagonal 
    // (all processes calculate the same number here)
    int elements_in_diagonal = std::min(start_col + 1, rows - start_row);

    // number of processes that get an extra cell to calculate 
    // (all processes calculate same num here)
    int z = elements_in_diagonal % num_processes;

    // min number of cells each process gets 
    // (all processes calculate same number here)
    int base_count = elements_in_diagonal / num_processes;

    // the starting index along the current diagonal that the current process needs to start at
    // this is the only number that is uniquely computed per process
    int start_idx = std::min(process_id, z) * (base_count + 1) + std::max(0, process_id - z) * (base_count);

    // put the results in the struct
    cellsToCalculate.startingRow = start_row + start_idx + 1;
    cellsToCalculate.startingColumn = start_col - start_idx + 1;
    cellsToCalculate.numElements = base_count + (process_id < (z) ? 1 : 0);
    cellsToCalculate.minPerProcess = base_count;
    cellsToCalculate.elementsInDiagonal = elements_in_diagonal;
}

void findLCS(const string& X, const string& Y, ProcessData& processData, int world_rank, int world_size) {
    timer t1;
    t1.start();
    // --------------------------------------------------------------
    int n = X.length();
    int m = Y.length();
    uint num_diagonals = n + m - 1;

    //vector<vector<int>> dp(n + 1, vector<int>(m + 1));
    int** dp = new int*[n + 1];
    for (int i = 0; i < (n + 1); i++) {
        dp[i] = new int[m + 1];
    }

    // First: broadcast 1 sequence to all nodes
    if (X[0] == Y[0]) {
        dp[0][0] = 1;
    } else {
        dp[0][0] = 0;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //vector<int> cellsToCalculate(3);
    CellsToCalculate cellsToCalculate();
    int x, y, nElements;

    for (int i = 1; i < num_diagonals; i++) {   

        getCellsToCalculate(cellsToCalculate, n, m, world_rank, world_size, i);
        x = cellsToCalculate.startingColumn;
        y = cellsToCalculate.startingRow;
        int initialX = x;
        int initialY = y;
        int finalX = x + (cellsToCalculate.numElements - 1);
        int finalY = y - (cellsToCalculate.numElements - 1);
        // num elements to process (always minElements || minElements + 1)
        nElements = cellsToCalculate.numElements;
        minElements = cellsToCalculate.minElementsPerProcess;
        elementsInDiagonal = cellsToCalculate.elementsInDiagonal;

        // guaranteed to have the values we need at this point
        for (int j = 0; j < nElements; j++) {
            if (X[x-1] == Y[y - 1]) {
                dp[x][y] = dp[x - 1][y - 1] + 1;
            } else {
                dp[x][y] = max(dp[x - 1][y], dp[x][y - 1]);
            }
            x += 1;
            y -= 1;
        }

        // SYNCHRONIZE
        // Send first cell computed to previous process
        // Send last cell computed to next process

        // where receiving/sending cell data is stored (position and value)
        // index 0: row
        // index 1: col
        // index 2: DP value
        int *lastCellReceive  = new int[3];
        int *firstCellReceive = new int[3];
        int *lastCellSend     = new int[3];
        int *firstCellSend    = new int[3];

        if (world_rank % 2 == 0) {                                          // even rank
            if (world_rank < world_size - 1) {                              // not last process
                // MPI_Send last cell to world_rank + 1
                lastCellSend[0] = finalX;
                lastCellSend[1] = finalY;
                lastCellSend[2] = dp[finalX][finalY];
                MPI_Send(
                    lastCellSend,    //data
                    3,               //count
                    MPI_INT,         //datatype,
                    world_rank+1,    //destination,
                    0,               //tag
                    MPI_COMM_WORLD   //comm
                );
                // MPI_Recv first cell from world_rank + 1
                MPI_Recv(
                    firstCellReceive, //data,
                    3,                       //count,
                    MPI_INT,                 //datatype,
                    world_rank+1,            //source,
                    0,            //tag,
                    MPI_COMM_WORLD,   //communicator,
                    MPI_STATUS_IGNORE //status
                );
            }
            if (world_rank > 0) {                                           // not first process
                // MPI_Recv last cell from world_rank - 1
                MPI_Recv(
                    lastCellReceive,
                    3,
                    MPI_INT,
                    world_rank-1,
                    0,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                )
                // MPI_Send first cell to world_rank - 1
                firstCellSend[0] = initialX;
                firstCellSend[1] = initialY;
                firstCellSend[2] = dp[initialX][initialY];
                MPI_Send(
                    firstCellSend,   //data
                    3,               //count
                    MPI_INT,         //datatype,
                    world_rank-1,    //destination,
                    0,               //tag
                    MPI_COMM_WORLD   //comm
                );
            }
        } else {                                                            // odd rank
            if (world_rank > 0) {                                           // not first process
                // MPI_Recv last cell from world_rank - 1
                MPI_Recv(
                    lastCellReceive,
                    3,
                    MPI_INT,
                    world_rank-1,
                    0,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );
                // MPI_Send first cell to world_rank - 1
                firstCellSend[0] = initialX;
                firstCellSend[1] = initialY;
                firstCellSend[2] = dp[initialX][initialY];
                MPI_Send(
                    firstCellSend,   //data
                    3,               //count
                    MPI_INT,         //datatype,
                    world_rank-1,    //destination,
                    0,               //tag
                    MPI_COMM_WORLD   //comm
                );
            }
            if (world_rank < world_size - 1) {                              // not last process
                // MPI_Send last cell to world_rank + 1
                lastCellSend[0] = finalX;
                lastCellSend[1] = finalY;
                lastCellSend[2] = dp[finalX][finalY];
                MPI_Send(
                    lastCellSend,    //data
                    3,               //count
                    MPI_INT,         //datatype,
                    world_rank+1,    //destination,
                    0,               //tag
                    MPI_COMM_WORLD   //comm
                );
                // MPI_Recv first cell from world_rank + 1
                MPI_Recv(
                    firstCellReceive,
                    3,
                    MPI_INT,
                    world_rank+1,
                    0,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );                
            }
        }

        // need to check that we arent filling row 0 or col 0, those values are always supposed to equal 0
        if (firstCellReceive[0] != 0 && firstCellReceive[1] != 0) {
            int firstCellX = firstCellReceive[0];
            int firstCellY = firstCellReceive[1];
            dp[firstCellX][firstCellY] = firstCellReceive[2];
        }        

        if (lastCellReceive[0] != 0 && lastCellReceive[1] != 0) {
            int lastCellX = lastCellReceive[0];
            int lastCellY = lastCellReceive[1];
            dp[lastCellX][lastCellY] = lastCellReceive[2];
        }
        
        delete[] lastCellReceive;
        delete[] firstCellReceive;
        delete[] lastCellSend;
        delete[] firstCellSend;

        //MPI_Barrier(MPI_COMM_WORLD); probably dont need barrier with synchronized sends/recv

    } // END OF DIAGONALS LOOP

    // Send data back to process 0


    unordered_set<string> lcsSet;

    // --------------------------------------------------------------
    processData.timeTaken = t1.stop();

    // free dp table
    for (int i = 0; i < (n + 1); i++) {
        delete[] dp[i];
    }
    delete[] dp;

    return;
}

int main(int argc, char** argv) {
    // === BEGIN cli parsing ===
    cxxopts::Options options("lcs", "Find longest common subsequence.");
    options.add_options(
        "",
        { 
            { "t", "Number of processes", cxxopts::value<uint>()->default_value(DEFAULT_PROCESSES) },
            { "x", "1st sequence", cxxopts::value<string>()->default_value(DEFAULT_STRING_X) },
            { "y", "2nd sequence", cxxopts::value<string>()->default_value(DEFAULT_STRING_Y) },
        }
    );
    auto cli = options.parse(argc, argv);

    // -t ignored in distributed
    uint nProcesess = cli["t"].as<uint>();
    if (LCS_VERSION == "serial" && nProcesess != 1) { // error when serial version > 1 process
        printf("Error : Number of processes for serial version must be 1.\n");
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
        printf("Number of processes : %u\n", world_size);
        printf("Sequence X : %s\n", X.c_str());
        printf("Sequence Y : %s\n", Y.c_str());
        printf("Finding longest common subsequence...\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // ------------------------------------------------------------
    timer mainTimer;
    mainTimer.start();

    ProcessData processData(world_rank, 0.0);

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
