import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np


def readfile(filename):
    x = []
    y = []
    z = []
    t = []
    with open(filename) as infile:
        for line in infile:
            columns = line.split()
            x.append(float(columns[0]))
            y.append(float(columns[1]))
            z.append(float(columns[2]))
            t.append(float(columns[3]))

    infile.close()
    return np.array(x), np.array(y), np.array(z), np.array(t)



def trajectory_plotter(filename, ax):
    x, y, z, t = readfile(filename)
    filename = filename.replace("_", " ").replace(".txt", " ").replace("p", "P")
    ax.plot3D(x, y, z, label=f"{filename}")
    size=30
    beginning = "green"
    end = "red"
    ax.scatter(x[0], y[0], z[0], s=size, color=beginning)
    ax.scatter(x[-1], y[-1], z[-1], s=size, color=end)

def pair_plotter(filename1, filename2):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_xlabel(r"x [$\mu$m]")
    ax.set_ylabel(r"y [$\mu$m]")
    ax.set_zlabel(r"z [$\mu$m]")
    trajectory_plotter(filename1, ax)
    trajectory_plotter(filename2, ax)
    plt.legend(loc="upper left")



filename1 = "particle_0_without_interactions.txt"
filename2 = "particle_0_with_interactions.txt"
filename3 = "particle_1_without_interactions.txt"
filename4 = "particle_1_with_interactions.txt"


filename5 = "particle_0_without_interactions_RK4.txt"
filename6 = "particle_0_with_interactions_RK4.txt"
filename7 = "particle_1_without_interactions_RK4.txt"
filename8 = "particle_1_with_interactions_RK4.txt"


if str(input("Do you want to see how particle 0 without interactions evolves in the z-direction over time? (y/n): ")) == "y":
    x, y, z, t = readfile(filename1)
    plt.plot(t, z)
    plt.xlabel(r"t [$\mu$s]")
    plt.ylabel(r"z [$\mu$m]")
    plt.show()


pair_plotter(filename2, filename6)
pair_plotter(filename4, filename8)
plt.show()




