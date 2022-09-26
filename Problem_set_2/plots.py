import matplotlib.pyplot as plt
import numpy as np


def readfile(filename):
    x = []
    v_1 = []
    v_2 = []
    v_3 = []
    with open(filename) as infile:
        for line in infile:
            columns = line.split()
            x.append(float(columns[0]))
            v_1.append(float(columns[1]))
            v_2.append(float(columns[2]))
            v_3.append(float(columns[3]))

    infile.close()
    return x, v_1, v_2, v_3

x,v_1,v_2,v_3 = readfile('exercise_6.txt')
plt.plot(x,v_1/np.amax(np.abs(v_1)),label=r'$\lambda_1$', linewidth=2, color='blue')
plt.plot(x,v_2/np.amax(np.abs(v_2)),label=r'$\lambda_2$', linewidth=2, color='green')
plt.plot(x,-1*np.array(v_3/np.amax(np.abs(v_3))),label=r'$\lambda_1$', linewidth=2, color='black')
plt.legend()
plt.savefig("eigenstates.pdf")


def test(n):
    h = 1./n
    a = -1./h**2
    d = 2./h**2
    lambda_ = np.zeros(n-1)
    eigenvecs = np.zeros((n-1,n-1))
    for i in range(n):
        lambda_[i-1] = d + 2.*a*np.cos(i*np.pi/n)
        for element in range(n):
            eigenvecs[i-1,element-1] = np.sin(element*i*np.pi/n)
    #print(lambda_)

    saved_indicies = np.zeros(3) #the three lowest eigenvalues
    lambda_val = 0
    for i in range(len(saved_indicies)):
        lambda_val = np.argmin(lambda_)
        saved_indicies[i] = lambda_val
        lambda_[lambda_val] = np.inf
    plt.ylabel('Normalized displacement')
    plt.xlabel('x/L')
    plt.plot(x[1:-1],np.transpose(eigenvecs)[int(saved_indicies[0])], '--', color='yellow', linewidth=1)
    plt.plot(x[1:-1],np.transpose(eigenvecs)[int(saved_indicies[1])], '--', color='yellow', linewidth=1)
    plt.plot(x[1:-1],np.transpose(eigenvecs)[int(saved_indicies[2])], '--', color='yellow', linewidth=1)
    plt.savefig("eigenstates_anal.pdf")
    plt.show()
test(100)
