#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

pair<int, string> LCS(const string &X, const string &Y) {
    int n = X.length();
    int m = Y.length();
    
    // Create a DP table to store lengths of longest common subsequence.
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
    
    // Fill the DP table
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            if (X[i - 1] == Y[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1] + 1;
            } else {
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
            }
        }
    }
    
    // The length of the LCS is stored at dp[n][m]
    int lcsLength = dp[n][m];
    
    // Backtrack to find the actual LCS string
    string lcs = "";
    int i = n, j = m;
    while (i > 0 && j > 0) {
        if (X[i - 1] == Y[j - 1]) {
            lcs = X[i - 1] + lcs;  // If characters match, add to the result
            i--;
            j--;
        } else if (dp[i - 1][j] > dp[i][j - 1]) {
            i--;  // Move up if the value above is greater
        } else {
            j--;  // Move left if the value to the left is greater
        }
    }
    
    // Return both the length of the LCS and the LCS string
    return {lcsLength, lcs};
}

int main() {
    string X = "cadbrz";
    string Y = "asbz";
    
    auto result = LCS(X, Y);
    
    cout << "Length of LCS is " << result.first << endl;
    cout << "LCS is " << result.second << endl;
    
    return 0;
}