import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interpolate
from tqdm import trange
plt.rcParams["text.usetex"] =True

def readfile(filename):
    x = []
    #v = []
    u = []
    with open(filename) as infile:
        for line in infile:
            columns = line.split()
            #v.append(float(columns[2]))
            x.append(float(columns[0]))
            u.append(float(columns[1]))
    infile.close()
    return x, u#, v
x, u = readfile('exercise_2.txt')
x2, v2 = readfile('exercise_7_10.txt')
x3, v3 = readfile('exercise_7_100.txt')
x4, v4 = readfile('exercise_7_1000.txt')
x5, v5 = readfile('exercise_7_10000.txt')
x6, v6 = readfile('exercise_7_100000.txt')
x7, v7 = readfile('exercise_7_1000000.txt')
x8, v8 = readfile('exercise_7_10000000.txt')


def u_(x_def):
    x_def = np.array(x_def)
    return 1-(1-np.exp(-10))*x_def - np.exp(-10*x_def)
#print(x2)
u2 = u_(x2)
u3 = u_(x3)
u4 = u_(x4)
u5 = u_(x5)
u6 = u_(x6)
u7 = u_(x7)
u8 = u_(x8)
"""plt.figure()
plt.plot(x, u, 'b-', label = 'Exact solution')
plt.xlabel('x')
plt.ylabel('u')
plt.show()"""
plt.figure()
plt.plot(x, u, 'b-', label = 'Exact solution')
plt.plot(x2, v2, 'r--', label = 'n = 1e1')
plt.plot(x3, v3, 'y--', label = 'n = 1e2')
plt.plot(x4, v4, 'g--', label = 'n = 1e3')
plt.plot(x5, v5, 'c:', label = 'n = 1e4')
plt.plot(x6, v6, 'k:', label = 'n = 1e5')
plt.plot(x7, v7, ':', label = 'n = 1e6')
plt.plot(x8, v8, ':', label = 'n = 1e7')
plt.xlabel('x')
plt.ylabel('u')
plt.legend()
#plt.show()
n_all = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
v_all = [v2, v3, v4, v5, v6, v7, v8]
u_all = [u2, u3, u4, u5, u6, u7, u8]
x_all = [x2, x3, x4, x5, x6, x7, x8]
max_err = np.zeros(len(u_all))
colors = ['b', 'c', 'r', 'g', 'm', 'k', 'y']
fig, ax = plt.subplots(2, 1, sharex=True)
ax[1].set_ylabel(r'$\epsilon$')
ax[0].set_yscale('log')
#ax[0].set_ylim(-10, 1)
ax[0].set_ylabel(r'$\Delta$')
#ax[1].set_ylim(-10, 1)
ax[1].set_yscale('log')
ax[1].set_xlabel('X')
ax[1].set_xlim(0, 1.5)
for i in range(len(u_all)):
    abs_error, rel_error, max_err = np.zeros(len(u_all[i])), np.zeros(len(u_all[i])), np.zeros(len(u_all))
    for j in range(1, len(u_all[i])-1):
        abs_error[j] = np.abs(u_all[i][j] - v_all[i][j])
        rel_error[j] = abs_error[j]/u_all[i][j]
    max_err[i] = np.amax(rel_error[1:-1])
    print(f'At n = {10**(i+1):.1e} the maximum relative error is {max_err[i]:.5e}')
    ax[0].plot(x_all[i][1:-1], abs_error[1:-1], label = f'n = {10**(i+1):.1e}', color=colors[i])
    ax[1].plot(x_all[i][1:-1], rel_error[1:-1], label = f'n = {10**(i+1):.1e}', color=colors[i])
plt.legend()
plt.show()
