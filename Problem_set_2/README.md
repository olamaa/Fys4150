# Fys4150

We solve the differential equation for the buckeling beam using the Jacobi rotation method.

In the cpp file named problem 2, the bools (line 12-16) make sure we create all the necessary txt files for the entire project. When you run the program, you will get an output of corresponding eigenvalues and eigenvectors both for the numerical and analytical in the terminal. Compile and run the cpp file named problem2 to obtain the .txt files and output (this works for mac):

Compile with: g++ problem2.cpp -o problem2 -std=c++11 -O2 -larmadillo

Run with: ./problem2

The Python file must be run to produce the plots. Use the following command to compile and run:

python3 plots.py

The Jacobi's method is partially borrowed from Mortens lecture notes
