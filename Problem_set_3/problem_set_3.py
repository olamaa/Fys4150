import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np


def readfile(filename):
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []
    t = []
    with open(filename) as infile:
        for line in infile:
            columns = line.split()
            x.append(float(columns[0]))
            y.append(float(columns[1]))
            z.append(float(columns[2]))
            vx.append(float(columns[3]))
            vy.append(float(columns[4]))
            vz.append(float(columns[5]))
            t.append(float(columns[6]))
    infile.close()
    return np.array(x), np.array(y), np.array(z), np.array(vx), np.array(vy), np.array(vz), np.array(t)


"""
def trajectory_plotter_xy(filename, ax):
    x, y, z, t = readfile(filename)
    filename = filename.replace("_", " ").replace(".txt", " ").replace("p", "P")
    plt.title(f"Trajectorty in the x-y plane for 50$\mu$s")
    plt.plot(x, y, label=f"{filename}")
    beginning = "green"
    end = "red"
    plt.scatter(x[0], y[0], s=size, color=beginning)
    plt.scatter(x[-1], y[-1], s=size, color=end)"""


def trajectory_plotter_3d(filename, ax):
    x, y, z, vx, vy, vz, t = readfile(filename)
    filename = filename.replace('1_','').replace('0_','').replace('particle_', '').replace('n_steps_400', '').replace("_", " ").replace(".txt", " ")
    ax.plot3D(x, y, z, label=f"{filename}")
    size=30
    beginning = "green"
    end = "red"
    ax.scatter(x[0], y[0], z[0], s=size, color=beginning)
    ax.scatter(x[-1], y[-1], z[-1], s=size, color=end)




def pair_plotter_3d(filename1, filename2, ax):
    ax.set_xlabel(r"x [$\mu$m]")
    ax.set_ylabel(r"y [$\mu$m]")
    ax.set_zlabel(r"z [$\mu$m]")
    trajectory_plotter_3d(filename1, ax)
    trajectory_plotter_3d(filename2, ax)
    plt.legend(loc="upper right")


def phase_space(r, v, ax, argument, title=False, points=True):
    if title != False:
        ax.set_title(title)
    ax.set_xlabel(fr"{argument} [$\mu$m]")
    ax.set_ylabel(fr"v$_{argument}$ [$\mu$m/$\mu$s]")
    ax.plot(r, v)
    if points==True:
        ax.scatter(r[0], v[0], color='green')
        ax.scatter(r[-1], v[-1], color='red')



p0wo = []
p0w = []
p1wo = []
p1w = []

p0wo_RK4 = []
p0w_RK4 = []
p1wo_RK4 = []
p1w_RK4 = []

for i in range(2, 6):
    p0wo.append("Euler_Cromer_particle_0_n_steps_" + str(2**i * 1000) + "_without_interactions.txt")
    p0w.append("Euler_Cromer_particle_0_n_steps_" + str(2**i * 1000) + "_with_interactions.txt")
    p1wo.append("Euler_Cromer_particle_1_n_steps_" + str(2**i * 1000) + "_without_interactions.txt")
    p1w.append("Euler_Cromer_particle_1_n_steps_" + str(2**i * 1000) + "_with_interactions.txt")

    p0wo_RK4.append("RK4_particle_0_n_steps_" + str(2**i * 1000) + "_without_interactions.txt")
    p0w_RK4.append("RK4_particle_0_n_steps_" + str(2**i * 1000) + "_with_interactions.txt")
    p1wo_RK4.append("RK4_particle_1_n_steps_" + str(2**i * 1000) + "_without_interactions.txt")
    p1w_RK4.append("RK4_particle_1_n_steps_" + str(2**i * 1000) + "_with_interactions.txt")


x0_wo_4000, y0_wo_4000, z0_wo_4000, vx0_wo_4000, vy0_wo_4000, vz0_wo_4000, t0_wo_4000 = readfile(p0wo_RK4[0])
x0_w_4000, y0_w_4000, z0_w_4000, vx0_w_4000, vy0_w_4000, vz0_w_4000, t0_w_4000 = readfile(p0w_RK4[0])
x1_wo_4000, y1_wo_4000, z1_wo_4000, vx1_wo_4000, vy1_wo_4000, vz1_wo_4000, t1_wo_4000 = readfile(p1wo_RK4[0])
x1_w_4000, y1_w_4000, z1_w_4000, vx1_w_4000, vy1_w_4000, vz1_w_4000, t1_w_4000 = readfile(p1w_RK4[0])


## Using RungeKutta, bulletpoint 1:
if str(input("Do you want to see how particle 0 without interactions evolves in the z-direction over time? (y/n): ")) == "y":
    plt.plot(t0_wo_4000, z0_wo_4000)
    plt.xlabel(r"t [$\mu$s]")
    plt.ylabel(r"z [$\mu$m]")
    #plt.show()

# bulletpoint 2:
input = str(input("Next, do you want to see with and without particle interactions? (y/n): "))
if input == "y":
    fig, ax = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    particle_x = [[x0_wo_4000, x0_w_4000], [x1_wo_4000, x1_w_4000]]
    particle_y = [[y0_wo_4000, y0_w_4000], [y1_wo_4000, y1_w_4000]]
    for i in range(len(particle_x[0])):
        ax[i].set_aspect('equal', 'box')
        ax[i].plot(particle_x[i][0], particle_y[i][0], label="Without interaction")
        ax[i].plot(particle_x[i][1], particle_y[i][1], label="With interaction")
        #ax[i].set_xlim(-80, 80)
        #ax[i].set_ylim(-80, 80)
        ax[i].set_xlabel(r"x [$\mu$m]")
        ax[i].set_ylabel(r"y [$\mu$m]")
        ax[i].set_title('Particle ' + str(i))
        ax[i].legend(loc="lower left")
    plt.axis('equal')
    #plt.show()


# bulletpoint 3:
input2 = 'y' #str(input("Next, do you want to the phase space plots? (y/n): "))
if input2 == 'y':
    fig, ax = plt.subplots(2, 2, figsize=(10, 8))
    beginning='green'
    end='red'
    phase_space(x0_wo_4000, vx0_wo_4000, ax[0,0], argument='x', title='Phase space without interaction')
    phase_space(x1_wo_4000, vx1_wo_4000, ax[0,0], argument='x')
    phase_space(x0_w_4000, vx0_w_4000, ax[0,1], argument='x', title='Phase space with interaction')
    phase_space(x1_w_4000, vx1_w_4000, ax[0,1], argument='x')
    phase_space(z0_wo_4000, vz0_wo_4000, ax[1,0], argument='z')
    phase_space(z1_wo_4000, vz1_wo_4000, ax[1,0], argument='z')
    phase_space(z0_w_4000, vz0_w_4000, ax[1,1], argument='z')
    phase_space(z1_w_4000, vz1_w_4000, ax[1,1], argument='z')
    #plt.show()

# bulletpoint 4:

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(2, 2, 1, projection='3d')
ax.set_title('Particle 0')
pair_plotter_3d(p0wo[0], p0wo_RK4[0], ax=ax) # particle 0 without interaction euler+rk4
ax = fig.add_subplot(2, 2, 3, projection='3d')
pair_plotter_3d(p0w[0], p0w_RK4[0], ax=ax) # particle 0 with interaction euler+rk4
ax = fig.add_subplot(2, 2, 2, projection='3d')
ax.set_title('Particle 1')
pair_plotter_3d(p1wo[0], p1wo_RK4[0], ax=ax) # particle 1 without interaction euler+rk4
ax = fig.add_subplot(2, 2, 4, projection='3d')
pair_plotter_3d(p1w[0], p1w_RK4[0], ax=ax) # particle 1 with interaction euler+rk4
#plt.show()



x0_wo_8000, y0_wo_8000, z0_wo_8000, vx0_wo_8000, vy0_wo_8000, vz0_wo_8000, t0_wo_8000 = readfile(p0wo_RK4[1])
x0_w_8000, y0_w_8000, z0_w_8000, vx0_w_8000, vy0_w_8000, vz0_w_8000, t0_w_8000 = readfile(p0w_RK4[1])
x1_wo_8000, y1_wo_8000, z1_wo_8000, vx1_wo_8000, vy1_wo_8000, vz1_wo_8000, t1_wo_8000 = readfile(p1wo_RK4[1])
x1_w_8000, y1_w_8000, z1_w_8000, vx1_w_8000, vy1_w_8000, vz1_w_8000, t1_w_8000 = readfile(p1w_RK4[1])


x0_wo_16000, y0_wo_16000, z0_wo_16000, vx0_wo_16000, vy0_wo_16000, vz0_wo_16000, t0_wo_16000 = readfile(p0wo_RK4[2])
x0_w_16000, y0_w_16000, z0_w_16000, vx0_w_16000, vy0_w_16000, vz0_w_16000, t0_w_16000 = readfile(p0w_RK4[2])
x1_wo_16000, y1_wo_16000, z1_wo_16000, vx1_wo_16000, vy1_wo_16000, vz1_wo_16000, t1_wo_16000 = readfile(p1wo_RK4[2])
x1_w_16000, y1_w_16000, z1_w_16000, vx1_w_16000, vy1_w_16000, vz1_w_16000, t1_w_16000 = readfile(p1w_RK4[2])


x0_wo_32000, y0_wo_32000, z0_wo_32000, vx0_wo_32000, vy0_wo_32000, vz0_wo_32000, t0_wo_32000 = readfile(p0wo_RK4[3])
x0_w_32000, y0_w_32000, z0_w_32000, vx0_w_32000, vy0_w_32000, vz0_w_32000, t0_w_32000 = readfile(p0w_RK4[3])
x1_wo_32000, y1_wo_32000, z1_wo_32000, vx1_wo_32000, vy1_wo_32000, vz1_wo_32000, t1_wo_32000 = readfile(p1wo_RK4[3])
x1_w_32000, y1_w_32000, z1_w_32000, vx1_w_32000, vy1_w_32000, vz1_w_32000, t1_w_32000 = readfile(p1w_RK4[3])



### With interaction
Ex0_wo_4000, Ey0_wo_4000, Ez0_wo_4000, Evx0_wo_4000, Evy0_wo_4000, Evz0_wo_4000, Et0_wo_4000 = readfile(p0wo[0])
Ex0_wo_8000, Ey0_wo_8000, Ez0_wo_8000, Evx0_wo_8000, Evy0_wo_8000, Evz0_wo_8000, Et0_wo_8000 = readfile(p0wo[1])
Ex0_wo_16000, Ey0_wo_16000, Ez0_wo_16000, Evx0_wo_16000, Evy0_wo_16000, Evz0_wo_16000, Et0_wo_16000 = readfile(p0wo[2])
Ex0_wo_32000, Ey0_wo_32000, Ez0_wo_32000, Evx0_wo_32000, Evy0_wo_32000, Evz0_wo_32000, Et0_wo_32000 = readfile(p0wo[3])

q = 1
m = 40.0775
v0_d_2 = 9.65
B0 = 9.65*10
omega_0 = q*B0/m
omega_z_square = 2*q*v0_d_2/m

def z_func(t):
    return 20*np.cos(np.sqrt(omega_z_square)*t)
def relative_error(x_num, y_num, z_num, x, y, z, t):
    r_numerical = np.sqrt(x_num**2 + y_num**2 + z_num**2)
    r_analytical = np.sqrt(x**2 + y**2 + z**2)
    return np.abs(r_analytical - r_numerical)/r_analytical

def f(initial_position, initial_velocity, t):
    omega_plus = (omega_0 + np.sqrt(omega_0**2 - 2*omega_z_square))/2
    omega_minus = (omega_0 - np.sqrt(omega_0**2 - 2*omega_z_square))/2
    A_plus = (initial_velocity + omega_minus*initial_position)/(omega_minus - omega_plus)
    A_minus = -(initial_velocity + omega_plus*initial_position)/(omega_minus - omega_plus)
    return A_plus*np.exp(-1.j*(omega_plus*t)) + A_minus*np.exp(-1.j*(omega_minus*t))

"""plt.figure()
time = np.linspace(0, 50, 4001)
x = np.real(f(20, 25, time))
y = np.imag(f(20, 25, time))
z = z_func(time)
plt.semilogy(t0_wo_4000, relative_error(x0_wo_4000, y0_wo_4000, z0_wo_4000, x, y, z, t0_wo_4000))
plt.show()"""

#str(input('Particle 1 wo interactions?'))
# Relative error for RungeKutta
fig, ax = plt.subplots(2, 2, figsize=(10,7), sharex=True)

ax[0][0].set_title('RungeKutta4')
ax[0][0].set_ylim(-0.1, 1.3)
ax[0][0].set_ylabel(r"Relative error [$r_{err}$]")
ax[0][0].plot(t0_wo_4000, relative_error(x0_wo_4000, y0_wo_4000, z0_wo_4000, np.real(f(20, 25, t0_wo_4000)), np.imag(f(20, 25, t0_wo_4000)), z_func(t0_wo_4000),  t0_wo_4000), label=r'$n=4\cdot 10^3$')
ax[0][0].plot(t0_wo_8000, relative_error(x0_wo_8000, y0_wo_8000, z0_wo_8000, np.real(f(20, 25, t0_wo_8000)), np.imag(f(20, 25, t0_wo_8000)), z_func(t0_wo_8000),  t0_wo_8000), label=r'$n=8\cdot 10^3$')
ax[0][0].plot(t0_wo_16000, relative_error(x0_wo_16000, y0_wo_16000, z0_wo_16000, np.real(f(20, 25, t0_wo_16000)), np.imag(f(20, 25, t0_wo_16000)), z_func(t0_wo_16000),  t0_wo_16000), label=r'$n=1.6\cdot 10^4$')
ax[0][0].plot(t0_wo_32000, relative_error(x0_wo_32000, y0_wo_32000, z0_wo_32000, np.real(f(20, 25, t0_wo_32000)), np.imag(f(20, 25, t0_wo_32000)), z_func(t0_wo_32000),  t0_wo_32000), label=r'$n=3.2\cdot 10^4$', linestyle='dashed')

ax[0][1].set_title('Forward Euler')
ax[0][1].set_ylim(-0.1, 1.3)
ax[0][1].plot(Et0_wo_4000, relative_error(Ex0_wo_4000, Ey0_wo_4000, Ez0_wo_4000, np.real(f(20, 25, Et0_wo_4000)), np.imag(f(20, 25, Et0_wo_4000)), z_func(Et0_wo_4000),  Et0_wo_4000), label=r'$n=4\cdot 10^3$')
ax[0][1].plot(Et0_wo_8000, relative_error(Ex0_wo_8000, Ey0_wo_8000, Ez0_wo_8000, np.real(f(20, 25, Et0_wo_8000)), np.imag(f(20, 25, Et0_wo_8000)), z_func(Et0_wo_8000),  Et0_wo_8000), label=r'$n=8\cdot 10^3$')
ax[0][1].plot(Et0_wo_16000, relative_error(Ex0_wo_16000, Ey0_wo_16000, Ez0_wo_16000, np.real(f(20, 25, Et0_wo_16000)), np.imag(f(20, 25, Et0_wo_16000)), z_func(Et0_wo_16000),  Et0_wo_16000), label=r'$n=1.6\cdot 10^4$')
ax[0][1].plot(Et0_wo_32000, relative_error(Ex0_wo_32000, Ey0_wo_32000, Ez0_wo_32000, np.real(f(20, 25, Et0_wo_32000)), np.imag(f(20, 25, Et0_wo_32000)), z_func(Et0_wo_32000),  Et0_wo_32000), label=r'$n=3.2\cdot 10^4$')

ax[1][0].set_ylim(10**(-15), 10**(1))
ax[1][0].set_ylabel(r"Relative error [$\log_{10}(r_{err})$]")
ax[1][0].set_xlabel(r"Time [$\mu$s]")
ax[1][0].semilogy(t0_wo_4000, relative_error(x0_wo_4000, y0_wo_4000, z0_wo_4000, np.real(f(20, 25, t0_wo_4000)), np.imag(f(20, 25, t0_wo_4000)), z_func(t0_wo_4000),  t0_wo_4000), label='n=4000')
ax[1][0].semilogy(t0_wo_8000, relative_error(x0_wo_8000, y0_wo_8000, z0_wo_8000, np.real(f(20, 25, t0_wo_8000)), np.imag(f(20, 25, t0_wo_8000)), z_func(t0_wo_8000),  t0_wo_8000), label='n=8000')
ax[1][0].semilogy(t0_wo_16000, relative_error(x0_wo_16000, y0_wo_16000, z0_wo_16000, np.real(f(20, 25, t0_wo_16000)), np.imag(f(20, 25, t0_wo_16000)), z_func(t0_wo_16000),  t0_wo_16000), label='n=16000')
ax[1][0].semilogy(t0_wo_32000, relative_error(x0_wo_32000, y0_wo_32000, z0_wo_32000, np.real(f(20, 25, t0_wo_32000)), np.imag(f(20, 25, t0_wo_32000)), z_func(t0_wo_32000),  t0_wo_32000), label='n=32000', linestyle='dashed')

ax[1][1].set_ylim(10**(-15), 10**(1))
ax[1][1].set_xlabel(r"Time [$\mu$s]")
ax[1][1].semilogy(Et0_wo_4000, relative_error(Ex0_wo_4000, Ey0_wo_4000, Ez0_wo_4000, np.real(f(20, 25, Et0_wo_4000)), np.imag(f(20, 25, Et0_wo_4000)), z_func(Et0_wo_4000),  Et0_wo_4000), label='n=4000')
ax[1][1].semilogy(Et0_wo_8000, relative_error(Ex0_wo_8000, Ey0_wo_8000, Ez0_wo_8000, np.real(f(20, 25, Et0_wo_8000)), np.imag(f(20, 25, Et0_wo_8000)), z_func(Et0_wo_8000),  Et0_wo_8000), label='n=8000')
ax[1][1].semilogy(Et0_wo_16000, relative_error(Ex0_wo_16000, Ey0_wo_16000, Ez0_wo_16000, np.real(f(20, 25, Et0_wo_16000)), np.imag(f(20, 25, Et0_wo_16000)), z_func(Et0_wo_16000),  Et0_wo_16000), label='n=16000')
ax[1][1].semilogy(Et0_wo_32000, relative_error(Ex0_wo_32000, Ey0_wo_32000, Ez0_wo_32000, np.real(f(20, 25, Et0_wo_32000)), np.imag(f(20, 25, Et0_wo_32000)), z_func(Et0_wo_32000),  Et0_wo_32000), label='n=32000')

ax[0][0].legend(loc='upper left')
ax[0][1].legend(loc='upper left')
#ax[0][1].legend(loc='lower right')
#ax[1][1].legend(loc='lower right')
plt.show()



"""
def r_error(z, n):
    sum = np.zeros(len(z))
    h = np.zeros(len(z))

    delta_max = np.zeros(len(z))
    for k in range(1, len(z)):
        h[k-1] = 50/n[k-1]
        h[k] = 50/n[k]
        t0 = np.linspace(0, 50, len(z[k-1]))
        t = np.linspace(0, 50, len(z[k]))
        delta_max[k-1] = np.max(z_func(t0, 20) - np.array(z[:][k-1]))
        delta_max[k] = np.max(z_func(t, 20) - z[:][k])
        sum[k] += np.log10(np.abs(delta_max[k]/delta_max[k-1]))/(np.log10(np.abs(h[k]/h[k-1])))
    return 1/3 * sum
#print(z0_wo_4000)
plt.figure()
plt.plot([4000, 8000, 16000, 32000], r_error(z = [z0_wo_4000, z0_wo_8000, z0_wo_16000, z0_wo_32000], n=[4000, 8000, 16000, 32000]))
plt.show()"""
