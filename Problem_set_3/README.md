# Fys4150

## Project 3

We study the Penning trap with the Euler and Runge-Kutta methods.

The main program must be compiled and run to produce the necessary txt files used for plotting in the python program. The main program can take multiple hours to run, but the majority of files will be made in a matter of 20 seconds. When the program tells you to abort (when the file 0.100000_test.txt is made), you should cancel the program. From there we rather use data produced by the Hyades compute cluster. The txt files we add in the GitHub repo is in order to save time.

Step 1) Compile with $ g++ -o main.exe penningtrap.cpp particle.cpp main.cpp -std=c++11 -O2 -larmadillo -I ./
Step 2) Run with $ ./main.exe

You should now have 32 files, where 16 are on the form "Forward_Euler" and "RK4" respectivly.

Step 3) Run  the Python program with $ python3 problem_set_3.py

From there on you will get a list of outputs you might want to see. just type
 $ y
 if it is so. These are all the results discussed in the report.
