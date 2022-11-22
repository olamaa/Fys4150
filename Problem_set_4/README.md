# Fys4150

## Project 4
We study the Ising model using the Markov Chains Monte Carlo method with and without parallellization.

The main program must be compiled and run to produce the necessary txt files used for plotting in the python program. The important parameters for saving time are "num_t" on line 236 (standard = 100) and "operations(L(i),T(j),1e6,"True","False","True");" on line 258. If you do not want to run the code for 1.5 hours, change "1e6" to "1e5", and "num_t" to 50. This should shorten the time used considerably.

Step 1) Compile with $ g++ -o working_para.exe working_para.cpp -std=c++11 -O2 -fopenmp -larmadillo
Step 2) Run with $ ./working_para.exe

These two commands should produce 5 different .txt files, and results from both the numerical 2x2 lattice, and the time spent on the 20x20 lattice with and without parallelization.

Step 3) Run the python program with $ python3 plots.py



