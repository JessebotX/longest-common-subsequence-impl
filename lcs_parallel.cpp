#include "core/cxxopts.h"
#include "core/get_time.h"
#include "core/utils.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <functional>
#include <thread>

#define LCS_VERSION       "parallel"
#define DEFAULT_THREADS   3
#define DEFAULT_STRING_X std::string(10000, 'A')
#define DEFAULT_STRING_Y std::string(10000, 'A')

using namespace std;

struct ThreadData {
    int id;
    double timeTaken;

    ThreadData(int _id, double _timeTaken)
        : id(_id), timeTaken(_timeTaken) {}
};

void getCellsToCalculate(vector<pair<int,int>>& vec, vector<vector<int>>& dp, int thread_id, uint num_threads, int diag) {
    // make vec empty
    vec.clear();

    // add indicies of cells to vec that thread_id needs to calculate
    int rows = dp.size() - 1;
    int cols = dp[0].size() - 1;  // Get the number of columns from first row
    
    // Find start position of the diagonal
    int start_row = std::max(0, diag - (cols - 1));
    int start_col = std::min(diag, cols - 1);
    
    // Calculate total number of elements in this diagonal
    int elements_in_diagonal = 0;
    int curr_row = start_row;
    int curr_col = start_col;
    while (curr_row < rows && curr_col >= 0) {
        elements_in_diagonal++;
        curr_row++;
        curr_col--;
    }
    if (elements_in_diagonal == 0) {
        return;
    }
    
    // Calculate the range for this thread
    int base_count = elements_in_diagonal / num_threads;  // Minimum each thread gets
    int extra = elements_in_diagonal % num_threads;       // Remaining elements
    
    // Calculate start and end indices for this thread
    int start_idx = thread_id * base_count + std::min(thread_id, extra);
    int end_idx = start_idx + base_count + (thread_id < extra ? 1 : 0);

    // Skip to this thread's starting position
    curr_row = start_row;
    curr_col = start_col;
    for (int i = 0; i < start_idx; i++) {
        curr_row++;
        curr_col--;
    }
    
    // Process elements for this thread
    for (int i = start_idx; i < end_idx && curr_row < rows && curr_col >= 0; i++) {
        vec.push_back({curr_row + 1, curr_col + 1});
        curr_row++;
        curr_col--;
    }
}

mutex print_mutex;

vector<string> findLCS(const string &X, const string &Y, vector<ThreadData> &threadData, int num_threads) {
    timer t1;
    t1.start();
    int n = X.length();
    int m = Y.length();
    
    CustomBarrier barrier(num_threads);

    uint num_diagonals = n + m - 1;

    // Create a DP table to store lengths of longest common subsequence
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));

    auto threadWorker = [&](int thread_id) {
        timer t;
        t.start();
        vector<pair<int, int>> cellsToCalculate;
        // fill dp table
        for (uint i = 0; i < num_diagonals; i++) {
            // fill cellsToCalculate with the cells this thread needs to calculate
            getCellsToCalculate(cellsToCalculate, dp, thread_id, num_threads, i);
            // calculate every cell in cellsToCalculate, updating dp
            for (const auto& pair : cellsToCalculate) {
                uint x = pair.first;
                uint y = pair.second;
                if (X[x-1] == Y[y-1]) {
                    dp[x][y] = dp[x - 1][y - 1] + 1;
                } else {
                    dp[x][y] = max(dp[x - 1][y], dp[x][y - 1]);
                }
            }

            barrier.wait();
        }
        threadData[thread_id].timeTaken = t.stop();

    };
    
    vector<thread> threads;

    // start and stop threads
    for (uint i = 0; i < num_threads; i++ ) {
        threads.emplace_back(threadWorker, i);
    }

    for (auto &t : threads) {
        t.join();
    }

    /*
    // Fill the DP table
    // Want to find dp[n][m]
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            if (X[i - 1] == Y[j - 1]) { 
                // If current characters equal, set DP table value to the previous diagonal DP entry + 1 (extending subsequence)
                dp[i][j] = dp[i - 1][j - 1] + 1;
            } 
            else { 
                // If current characters not equal, set DP table value to the maximum of:
                // The value from the cell directly above, ignoring the current character of X
                // The value from the cell directly to the left, ignoring the current character of Y
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
            }
        }
    }
    */

    // Set to store unique LCS
    unordered_set<string> lcsSet;
    
    // Helper function to backtrack and find all LCS
    function<void(int, int, string)> backtrack = [&](int i, int j, string currentLCS) {
        if (i == 0 || j == 0) {
            // If we reach the top or left boundary of DP table, we have found an LCS
            lcsSet.insert(currentLCS);
            return;
        }
        
        if (X[i - 1] == Y[j - 1]) {
            // If characters match, add the character to the LCS and move diagonally
            backtrack(i - 1, j - 1, X[i - 1] + currentLCS);
        } else {
            // If characters do not match, move in both directions (up and left)
            if (dp[i - 1][j] == dp[i][j]) {
                backtrack(i - 1, j, currentLCS);
            }
            if (dp[i][j - 1] == dp[i][j]) {
                backtrack(i, j - 1, currentLCS);
            }
        }
    };

    // Start at endpoint and work backwards
    backtrack(n, m, "");
    //threadData[thread_id].timeTaken = t1.stop();
    
    // Return a vector that contains all LCS
    vector<string> result(lcsSet.begin(), lcsSet.end());

    return result;
}

int main(int argc, char **argv) {
    // Begin CLI parsing
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

    uint nThreads = cli["t"].as<uint>();
    if (LCS_VERSION == "serial" && nThreads != 1) { // error when serial version > 1 thread
        printf("Error : Number of threads for serial version must be 1.\n");
        return 1;
    }

    string X = cli["x"].as<string>(); // sequence 1
    string Y = cli["y"].as<string>(); // sequence 2
    
    // End of CLI parsing

    printf("LCS Version : %s\n", LCS_VERSION);
    printf("Number of threads : %u\n", nThreads);
    printf("Sequence X : %s\n", X.c_str());
    printf("Sequence Y : %s\n", Y.c_str());
    printf("Finding longest common subsequence...\n");

    // --------------------------------------------------------------------------------
    timer mainTimer;
    mainTimer.start();

    // Initialize thread data storage
    vector<ThreadData> threadDataList;
    for (int i = 0; i < nThreads; i++) {
        threadDataList.push_back(ThreadData(i, 0.0));
    }

    
    // Find LCS
    vector<string> lcsResults = findLCS(X, Y, threadDataList, nThreads);

    // --------------------------------------------------------------------------------
    double timeTaken = mainTimer.stop();

    // Print all found LCS sequences
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

    // Print thread data
    printf("id, time_taken\n");
    for (int i = 0; i < nThreads; i++) {
        printf("%d, %.4f\n", threadDataList[i].id, threadDataList[i].timeTaken);
    }

    printf("Time taken (in seconds) : %.4f\n", timeTaken);

    return 0;
}
