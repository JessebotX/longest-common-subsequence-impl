#include "core/cxxopts.h"
#include "core/get_time.h"
#include "core/utils.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <functional>

#define LCS_VERSION       "serial"
#define DEFAULT_THREADS   1
#define DEFAULT_STRING_X  "abcd"
#define DEFAULT_STRING_Y  "acbadabc"

using namespace std;

struct ThreadData {
    int id;
    double timeTaken;

    ThreadData(int _id, double _timeTaken)
        : id(_id), timeTaken(_timeTaken) {}
};

vector<string> findLCS(const string &X, const string &Y, ThreadData &threadData) {
    timer t1;
    t1.start();
    int n = X.length();
    int m = Y.length();
    
    // Create a DP table to store lengths of longest common subsequence
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
    
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
    threadData.timeTaken = t1.stop();
    
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
    printf("Number of threads : %zd\n", nThreads);
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
    vector<string> lcsResults = findLCS(X, Y, threadDataList[0]);

    // --------------------------------------------------------------------------------
    double timeTaken = mainTimer.stop();

    // Print all found LCS sequences
    size_t nSubsequences = lcsResults.size();
    if (nSubsequences > 0) {
        printf("LCS length : %zd\n", nSubsequences);

        for (const string &lcs : lcsResults) {
            printf("Length %zd subsequence : %s\n", nSubsequences, lcs.c_str());
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
