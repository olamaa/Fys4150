import numpy as np
import matplotlib.pyplot as plt


def readfile(filename):
    x = []
    u = []
    with open(filename) as infile:
        for line in infile:
            columns = line.split()
            x.append(float(columns[0]))
            u.append(float(columns[1]))
    infile.close()
    return x,u

x,u = readfile('0.100000_with_hyades_full_resolved.txt')
x1,u1 = readfile('0.400000_with_hyades_full_resolved.txt')
x2,u2 = readfile('0.700000_with_hyades_full_resolved.txt')
plt.plot(x,u)
plt.plot(x1,u1)
plt.plot(x2,u2)
#plt.savefig("resolved.pdf")

plt.show()

#x,u = readfile('x_u.txt')
#x_v,v = readfile('x_v.txt')
#plt.plot(x,u,label='u')
#plt.plot(x_v,v,'--')
#plt.xlabel('x')
#plt.ylabel('u(x)') 
#plt.savefig("u(x).pdf")
##plt.plot(x,v,'--')
#plt.legend()
#plt.show()
