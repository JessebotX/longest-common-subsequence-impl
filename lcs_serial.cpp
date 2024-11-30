#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_set>

using namespace std;

vector<string> findLCS(const string &X, const string &Y) {
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
    
    // Return a vector that contains all LCS
    vector<string> result(lcsSet.begin(), lcsSet.end());
    return result;
}

int main() {
    string X = "abcd";
    string Y = "acbad";
    
    // Find all LCS
    vector<string> lcsResults = findLCS(X, Y);
    
    // Print all found LCS sequences
    cout << "LCS sequences:" << endl;
    for (const string &lcs : lcsResults) {
        cout << lcs << endl;
    }

    return 0;
}