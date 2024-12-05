# CMPT 431 Project for Longest Common Subsequence (LCS) Problem

Testing is done on the SFU cluster/slurm workload manager.

Building
--------

Requirements:
- g++
- mpic++
- make
- python3

To build executables run:

    make

Which will create the following:

- lcs_serial: serial version of the LCS problem
- lcs_parallel: thread-based implementation of LCS
- lcs_distributed: MPI implementation of LCS

Running
-------

    -t <number of threads>
    -x <1st sequence string>
    -y <2nd sequence string>

- lcs_distributed does not implement -t
- lcs_serial -t only accepts 1
- x and y sequences accept symbols and any-case letters
  - lowercase and uppercase are treated as separate characters 
- When no flags are provided, default values will be used.

Examples:

    ./lcs_serial
    ./lcs_serial      -x abcd -y acbad
    ./lcs_serial      -t 1    -x hello
    ./lcs_parallel    -t 4    -x abcd  -y abcde

Testing
-------

Test scripts are provided for convenience.

Run on the cluster with sbatch <batch file name>

i.e. sbatch test-slurm-mpi-2.sh
