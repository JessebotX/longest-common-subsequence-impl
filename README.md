# CMPT 431 Project for Longest Common Subsequence (LCS) Problem

Building
--------

Requirements:
- g++
- mpic++
- make

To build executables run:

    make

Which will create the following:

- lcs_serial: serial version of the LCS problem
- lcs_parallel: thread-based implementation of LCS
- lcs_distributed: MPI implementation of LCS


Running
-------

All executables support the same flags.

    -t <number of threads>
    -x <1st sequence string>
    -y <2nd sequence string>

- When no flags are provided, default values will be used.
- lcs_serial -t only accepts 1
- lcs_distributed ignores -t
- x and y sequences accept symbols and any-case letters
  - lowercase and uppercase are treated as separate characters 

Examples:

    ./lcs_serial
    ./lcs_serial      -x abcd -y acbad
    ./lcs_serial      -t 1    -x hello
    ./lcs_parallel    -t 4    -x abcd  -y abcde
    ./lcs_distributed -x abde -y abcde

Testing
-------

Test scripts are provided for convenience.
