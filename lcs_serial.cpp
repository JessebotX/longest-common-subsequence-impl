#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

string findLongestCommonSubstring(const string& s1, const string& s2) {
    int m = s1.length();
    int n = s2.length();

    // Create a 2D Dynamic Programming table
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    int maxLength = 0;    // Length of longest common substring found
    int endIdx = 0;       // End index of the longest common substring in s1

    // Fill the DP table iteratively
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (s1[i - 1] == s2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1] + 1;
                if (dp[i][j] > maxLength) {
                    maxLength = dp[i][j];
                    endIdx = i - 1;  // Track the end index of the common substring in s1
                }
            } else {
                dp[i][j] = 0;
            }
        }
    }

    // If no common substring was found, return an empty string
    if (maxLength == 0) return "";

    // Extract the longest common substring
    return s1.substr(endIdx - maxLength + 1, maxLength);
}

int main() {
    string s1, s2;
    
    // Input strings
    cout << "Enter the first string: ";
    cin >> s1;
    cout << "Enter the second string: ";
    cin >> s2;
    
    string result = findLongestCommonSubstring(s1, s2);
    
    // Output the result
    if (result.empty()) {
        cout << "No common substring found.\n";
    } else {
        cout << "Longest common substring: " << result << endl;
    }
    
    return 0;
}