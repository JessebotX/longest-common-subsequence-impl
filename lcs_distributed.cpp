#include "core/utils.h"
#include <mpi.h>
#include <stdio.h>

#define LCS_VERSION         "distributed"
#define DEFAULT_STRING_X  "GACAT"
#define DEFAULT_STRING_Y  "ACCGATCG"
#define PRINT_TASK          0

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
};

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

vector<string> findLCS(const string& X, const string& Y, ProcessData& processData, int world_rank, int world_size) {
    timer t1;
    t1.start();
    // --------------------------------------------------------------
    int n = X.length();
    int m = Y.length();
    uint num_diagonals = n + m - 1;

    vector<vector<int>> dp(n + 1, vector<int>(m + 1, -100));
    // keep row 0 and column 0 equal to 0, the rest -100 (using impossible value as empty marker)
    std::fill(dp[0].begin(), dp[0].end(), 0);
    for(int i = 0; i < n + 1; i++) {
        dp[i][0] = 0;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    CellsToCalculate cellsToCalculate;
    int x, y, nElements;

    for (int i = 0; i < num_diagonals; i++) {   

        getCellsToCalculate(cellsToCalculate, n, m, world_rank, world_size, i);
        x = cellsToCalculate.startingRow;
        y = cellsToCalculate.startingColumn;
        int initialX = x;
        int initialY = y;
        int finalX = x + (cellsToCalculate.numElements - 1);
        int finalY = y - (cellsToCalculate.numElements - 1);
        // num elements to process (always minElements || minElements + 1)
        int nElements = cellsToCalculate.numElements;
        int minElements = cellsToCalculate.minPerProcess;
        int elementsInDiagonal = cellsToCalculate.elementsInDiagonal;

        // if (world_rank == PRINT_TASK) {
        //     printf("\n");
        //     printf("Process: %d\n", world_rank);
        //     printf("(diag, x, y): (%d, %d, %d)\n", i, x, y);
        //     printf("(finalX, finalY): (%d, %d)\n", finalX, finalY);
        //     printf("nElements: %d\n", nElements);
        //     printf("elementsInDiagonal: %d\n", elementsInDiagonal);
            
        // }

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

        // if (world_rank == PRINT_TASK) printf("SYNCHRONIZING\n");
        // SYNCHRONIZE
        // Send first cell computed to previous process
        // Send last cell computed to next process

        // where receiving/sending cell data is stored (position and value)
        // index 0: row
        // index 1: col
        // index 2: DP value
        // if (world_rank == PRINT_TASK) printf("CREATE INT[3]\n");
        int *lastCellReceive  = new int[3];
        int *firstCellReceive = new int[3];
        int *lastCellSend     = new int[3];
        int *firstCellSend    = new int[3];
        std::fill(lastCellReceive, lastCellReceive + 3, 0);
        std::fill(firstCellReceive, firstCellReceive + 3, 0);
        std::fill(lastCellSend, lastCellSend + 3, 0);
        std::fill(firstCellSend, firstCellSend + 3, 0);
        
        // This section of code is used to determine whether or not all the processes
        // have actual work to do. For instance, in diag 2 there may not be enough elements
        // for 8 processes, etc.
        //
        // In this scenario we skip trying to make these inactive processes send/recv
        int world_size_2;
        if (elementsInDiagonal < world_size) {
            world_size_2 = elementsInDiagonal;
        } else {
            world_size_2 = world_size;
        }

        if (world_rank < world_size_2) {
            if (world_rank % 2 == 0) {                                          // even rank
                if (world_rank < world_size_2 - 1) {                              // not last process
                    // MPI_Send last cell to world_rank + 1
                    lastCellSend[0] = finalX;
                    lastCellSend[1] = finalY;
                    lastCellSend[2] = dp[finalX][finalY];
                    // if(world_rank == PRINT_TASK) printf("SENDING lastCellSend TO WORLD_RANK + 1: (%d, %d, %d)\n", lastCellSend[0], lastCellSend[1], lastCellSend[2]);
                    MPI_Send(
                        lastCellSend,    //data
                        3,               //count
                        MPI_INT,         //datatype,
                        world_rank+1,    //destination,
                        0,               //tag
                        MPI_COMM_WORLD   //comm
                    );
                    // if (world_rank == PRINT_TASK) printf("DONE SENDING lastCellSend TO WORLD_RANK + 1\n");
                    // MPI_Recv first cell from world_rank + 1
                    // if(world_rank == PRINT_TASK) printf("RECEIVING firstCellReceive FROM WORLD_RANK + 1\n");
                    MPI_Recv(
                        firstCellReceive, //data,
                        3,                       //count,
                        MPI_INT,                 //datatype,
                        world_rank+1,            //source,
                        0,            //tag,
                        MPI_COMM_WORLD,   //communicator,
                        MPI_STATUS_IGNORE //status
                    );
                    // if(world_rank == PRINT_TASK) printf("DONE RECEIVING firstCellReceive FROM WORLD_RANK + 1: (%d, %d, %d)\n", firstCellReceive[0], firstCellReceive[1], firstCellReceive[2]);
                }
                if (world_rank > 0) {                                           // not first process
                    // MPI_Recv last cell from world_rank - 1
                    // if(world_rank == PRINT_TASK) printf("RECEIVING lastCellReceive FROM WORLD_RANK - 1\n");
                    MPI_Recv(
                        lastCellReceive,
                        3,
                        MPI_INT,
                        world_rank-1,
                        0,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                    );
                    // if(world_rank == PRINT_TASK) printf("DONE RECEIVING lastCellReceive FROM WORLD_RANK - 1: (%d, %d, %d)\n", lastCellReceive[0], lastCellReceive[1], lastCellReceive[2]);
                    
                    // MPI_Send first cell to world_rank - 1
                    firstCellSend[0] = initialX;
                    firstCellSend[1] = initialY;
                    firstCellSend[2] = dp[initialX][initialY];
                    // if(world_rank == PRINT_TASK) printf("SENDING firstCellSend TO WORLD_RANK - 1: (%d, %d, %d)\n", firstCellSend[0], firstCellSend[1], firstCellSend[2]);
                    MPI_Send(
                        firstCellSend,   //data
                        3,               //count
                        MPI_INT,         //datatype,
                        world_rank-1,    //destination,
                        0,               //tag
                        MPI_COMM_WORLD   //comm
                    );
                    // if(world_rank == PRINT_TASK) printf("DONE SENDING firstCellSend TO WORLD_RANK - 1\n");
                }
            } else {                                                            // odd rank
                if (world_rank > 0) {                                           // not first process
                    // MPI_Recv last cell from world_rank - 1
                    // if(world_rank == PRINT_TASK) printf("RECEIVING lastCellReceive FROM WORLD_RANK - 1\n");
                    MPI_Recv(
                        lastCellReceive,
                        3,
                        MPI_INT,
                        world_rank-1,
                        0,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                    );
                    // if(world_rank == PRINT_TASK) printf("DONE RECEIVING lastCellReceive FROM WORLD_RANK - 1: (%d, %d, %d)\n", lastCellReceive[0], lastCellReceive[1], lastCellReceive[2]);
                    // MPI_Send first cell to world_rank - 1
                    firstCellSend[0] = initialX;
                    firstCellSend[1] = initialY;
                    firstCellSend[2] = dp[initialX][initialY];
                    // if(world_rank == PRINT_TASK) printf("SENDING firstCellSend TO WORLD_RANK - 1: (%d, %d, %d)\n", firstCellSend[0], firstCellSend[1], firstCellSend[2]);
                    MPI_Send(
                        firstCellSend,   //data
                        3,               //count
                        MPI_INT,         //datatype,
                        world_rank-1,    //destination,
                        0,               //tag
                        MPI_COMM_WORLD   //comm
                    );
                    // if(world_rank == PRINT_TASK) printf("DONE SENDING firstCellSend TO WORLD_RANK - 1\n");
                    
                }
                if (world_rank < world_size_2 - 1) {                              // not last process
                    // MPI_Send last cell to world_rank + 1
                    lastCellSend[0] = finalX;
                    lastCellSend[1] = finalY;
                    lastCellSend[2] = dp[finalX][finalY];
                    // if(world_rank == PRINT_TASK) printf("SENDING lastCellSend TO WORLD_RANK + 1: (%d, %d, %d)\n", lastCellSend[0], lastCellSend[1], lastCellSend[2]);
                    MPI_Send(
                        lastCellSend,    //data
                        3,               //count
                        MPI_INT,         //datatype,
                        world_rank+1,    //destination,
                        0,               //tag
                        MPI_COMM_WORLD   //comm
                    );
                    // if(world_rank == PRINT_TASK) printf("DONE SENDING lastCellSend TO WORLD_RANK + 1\n");
                    // MPI_Recv first cell from world_rank + 1
                    // if(world_rank == PRINT_TASK) printf("RECEIVING firstCellReceive FROM WORLD_RANK + 1\n");
                    MPI_Recv(
                        firstCellReceive,
                        3,
                        MPI_INT,
                        world_rank+1,
                        0,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                    );
                    // if(world_rank == PRINT_TASK) printf("DONE RECEIVING firstCellReceive FROM WORLD_RANK + 1: (%d, %d, %d)\n", firstCellReceive[0], firstCellReceive[1], firstCellReceive[2]);             
                }
            }

            // if (world_rank == PRINT_TASK) {
            //     printf("UPDATE DP\n");
            //     printf("changing dp[%d][%d] to %d\n", firstCellReceive[0], firstCellReceive[1], firstCellReceive[2]);
            //     printf("changing dp[%d][%d] to %d\n", lastCellReceive[0], lastCellReceive[1], lastCellReceive[2]);
            // }

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
            // if (world_rank == PRINT_TASK) printf("DELETE INT[3]\n");
            delete[] lastCellReceive;
            delete[] firstCellReceive;
            delete[] lastCellSend;
            delete[] firstCellSend;
            // if (world_rank == PRINT_TASK) printf("DONE DELETE INT[3]\n");
        }
        MPI_Barrier(MPI_COMM_WORLD); //probably dont need barrier with synchronized sends/recv
        // if (world_rank == PRINT_TASK) printf("END DIAGONAL %d\n", i);

    } // END OF DIAGONALS LOOP

    // Send data back to process 0

    int* rowReceive = new int[m];
    int* rowSend = new int[m];

    // Reconstruct dp at process 0

    // for each row i
    for (int i = 2; i < n+1; i++) {
        
        // if we are the root process, receive all rows
        if(world_rank == 0) {
            for (int j = 1; j < world_size; j++) {
                // receive row i from that process j
                // printf("Receiving\n");
                MPI_Recv(
                    rowReceive,
                    m,
                    MPI_INT,
                    j,
                    0,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );
                // swap out blank cells in our DP table for the once we just received
                for (int k = 0; k < m; k++) {
                    //printf("looking at dp[%d][%d]\n", i, k+1);
                    if (dp[i][k+1] < rowReceive[k]) {
                        dp[i][k+1] = rowReceive[k];
                    }
                }
            }
        // if we arent the root process, send row i
        } else {
            
            std::copy(dp[i].begin() + 1, dp[i].end(), rowSend);

            MPI_Send( 
                rowSend,//dp is a vector, not an int arr anymore
                m,
                MPI_INT,
                0, // this is destination right?? it should be
                0,
                MPI_COMM_WORLD
            );
            
        }
    }
    delete[] rowSend;
    delete[] rowReceive;

    
    // PRINT DP TABLE
    if(world_rank == PRINT_TASK) {
        printf("        ");
        for(int i = 0; i < m; i++) {
            printf("%3c ", Y[i]);
        }
        printf("\n");
        for (int i = 0; i < dp.size(); i ++ ) {
            if(i > 0) {
                printf("%3c ", X[i-1]);
            } else {
                printf("    ");
            }
            for (int j = 0; j < dp[0].size(); j++ ) {
                printf("%3d ", dp[i][j]);
            }
            printf("\n");
            
        }
    }
    

    unordered_set<string> lcsSet;

    
    function<void(int, int, string)> backtrack = [&](int i, int j, string currentLCS) {
        printf("entering backtrack at %d %d DP: %d\n", i, j, dp[i][j]);
        if (i == 0 || j == 0) {
            // If we've found a valid LCS that is non-empty, insert it
            if (!currentLCS.empty()) {
                lcsSet.insert(currentLCS);
                printf("inserted %s\n", currentLCS.c_str());
            }
            return;
        }
        
        if (X[i - 1] == Y[j - 1]) {
            // If characters match, add the character to the LCS and move diagonally
            backtrack(i - 1, j - 1, X[i - 1] + currentLCS);
        } else {
            // If characters do not match, move in both directions (up and left)
            if (dp[i - 1][j] == dp[i][j]) {
                printf("backtrack %d %d %s\n", i-1, j, currentLCS.c_str());
                backtrack(i - 1, j, currentLCS);
            }
            if (dp[i][j - 1] == dp[i][j]) {
                printf("backtrack %d %d %s\n", i, j-1, currentLCS.c_str());
                backtrack(i, j - 1, currentLCS);
            }
        }
    };

    processData.timeTaken = t1.stop();

    if (world_rank == 0) {
        printf("backtrack %d %d %s\n", n, m, "");
        backtrack(n, m, "");
        if (lcsSet.empty()) {
            return {};  // Return an empty vector if no LCS exists
        }
        vector<string> result(lcsSet.begin(), lcsSet.end());
        return result;
    } else {
        return {};
    }

    

    
}

int main(int argc, char** argv) {
    // === BEGIN cli parsing ===
    cxxopts::Options options("lcs", "Find longest common subsequence.");
    options.add_options(
        "",
        {
            { "x", "1st sequence", cxxopts::value<string>()->default_value(DEFAULT_STRING_X) },
            { "y", "2nd sequence", cxxopts::value<string>()->default_value(DEFAULT_STRING_Y) },
        }
    );
    auto cli = options.parse(argc, argv);

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
        printf("Number of processes : %d\n", world_size);
        printf("Sequence X : %s\n", X.c_str());
        printf("Sequence Y : %s\n", Y.c_str());
        printf("Finding longest common subsequence...\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // ------------------------------------------------------------
    timer mainTimer;
    mainTimer.start();

    ProcessData processData(world_rank, 0.0);

    vector<string> lcsResults = findLCS(X, Y, processData, world_rank, world_size);

    // ------------------------------------------------------------
    double timeTaken = mainTimer.stop();

    MPI_Barrier(MPI_COMM_WORLD);
    // ------ BEGIN print ------

    if (world_rank == 0) {
        size_t nSubsequences = lcsResults.size();
        if (nSubsequences > 0) {
            size_t lcsLength = lcsResults[0].size();
            printf("LCS length : %zd\n", lcsLength);

            for (const string &lcs : lcsResults) {
                printf("Length %zd subsequence : %s\n", lcsLength, lcs.c_str());
            }
        } else {
            printf("No subsequence found\n");
        }
    }

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
