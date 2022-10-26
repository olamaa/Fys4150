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
    p0wo.append("Forward_Euler_particle_0_n_steps_" + str(2**i * 1000) + "_without_interactions.txt")
    p0w.append("Forward_Euler_particle_0_n_steps_" + str(2**i * 1000) + "_with_interactions.txt")
    p1wo.append("Forward_Euler_particle_1_n_steps_" + str(2**i * 1000) + "_without_interactions.txt")
    p1w.append("Forward_Euler_particle_1_n_steps_" + str(2**i * 1000) + "_with_interactions.txt")

    p0wo_RK4.append("RK4_particle_0_n_steps_" + str(2**i * 1000) + "_without_interactions.txt")
    p0w_RK4.append("RK4_particle_0_n_steps_" + str(2**i * 1000) + "_with_interactions.txt")
    p1wo_RK4.append("RK4_particle_1_n_steps_" + str(2**i * 1000) + "_without_interactions.txt")
    p1w_RK4.append("RK4_particle_1_n_steps_" + str(2**i * 1000) + "_with_interactions.txt")


x0_wo_4000, y0_wo_4000, z0_wo_4000, vx0_wo_4000, vy0_wo_4000, vz0_wo_4000, t0_wo_4000 = readfile(p0wo_RK4[0])
x0_w_4000, y0_w_4000, z0_w_4000, vx0_w_4000, vy0_w_4000, vz0_w_4000, t0_w_4000 = readfile(p0w_RK4[0])
x1_wo_4000, y1_wo_4000, z1_wo_4000, vx1_wo_4000, vy1_wo_4000, vz1_wo_4000, t1_wo_4000 = readfile(p1wo_RK4[0])
x1_w_4000, y1_w_4000, z1_w_4000, vx1_w_4000, vy1_w_4000, vz1_w_4000, t1_w_4000 = readfile(p1w_RK4[0])


## Using RungeKutta, bulletpoint 1:
if str(input("movement of particle 0 without interactions in the zt-plane (y/n): "))== "y":
    plt.figure(figsize=(7,5))
    plt.plot(t0_wo_4000, z0_wo_4000)
    plt.xlabel(r"t [$\mu$s]")
    plt.ylabel(r"z [$\mu$m]")
    plt.show()

# bulletpoint 2:
if str(input("Particles in the xy-plane with and without interaction (y/n): "))== "y":
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
    plt.show()


# bulletpoint 3:
if str(input("Phase spaces with and without interaction (y/n): "))== "y":
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
    plt.show()

# bulletpoint 4:
if str(input("3D paths of both particles with and without interaction? (y/n): "))== "y":
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
    plt.show()



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
mass_Ca_I = 40.078
c = 299792458
atomic_mass = 1.66053906660*10**(-27)
elementary_charge = 1.60217663*1e-19
m = (11.87172*elementary_charge/(c**2*atomic_mass) + mass_Ca_I)-5.4858*10**(-4)
V = 9.64852558*10**7
V0 = 25*10**(-3)*V
T = 9.64852558*10
d = 500
v0_d_2 = V0/d**2
B0 = 1*T
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


#str(input('Particle 1 wo interactions?'))
# Relative error for RungeKutta
if str(input("Relative errors for the RK4 and forward Euler (y/n): "))== "y":
    fig, ax = plt.subplots(2, 2, figsize=(10,7), sharex=True)

    ax[0][0].set_title('RungeKutta4')
    ax[0][0].set_ylim(-0.1, 1.3)
    ax[0][0].set_ylabel(r"Relative error [$r_{err}$]")
    ax[0][0].plot(t0_wo_4000, relative_error(x0_wo_4000, y0_wo_4000, z0_wo_4000, np.real(f(20, 25, t0_wo_4000)), np.imag(f(20, 25, t0_wo_4000)), z_func(t0_wo_4000),  t0_wo_4000), label=r'$n_1=4\cdot 10^3$')
    ax[0][0].plot(t0_wo_8000, relative_error(x0_wo_8000, y0_wo_8000, z0_wo_8000, np.real(f(20, 25, t0_wo_8000)), np.imag(f(20, 25, t0_wo_8000)), z_func(t0_wo_8000),  t0_wo_8000), label=r'$n_2=8\cdot 10^3$')
    ax[0][0].plot(t0_wo_16000, relative_error(x0_wo_16000, y0_wo_16000, z0_wo_16000, np.real(f(20, 25, t0_wo_16000)), np.imag(f(20, 25, t0_wo_16000)), z_func(t0_wo_16000),  t0_wo_16000), label=r'$n_3=1.6\cdot 10^4$')
    ax[0][0].plot(t0_wo_32000, relative_error(x0_wo_32000, y0_wo_32000, z0_wo_32000, np.real(f(20, 25, t0_wo_32000)), np.imag(f(20, 25, t0_wo_32000)), z_func(t0_wo_32000),  t0_wo_32000), label=r'$n_4=3.2\cdot 10^4$', linestyle='dashed')

    ax[0][1].set_title('Forward Euler')
    ax[0][1].set_ylim(-0.1, 1.3)
    ax[0][1].plot(Et0_wo_4000, relative_error(Ex0_wo_4000, Ey0_wo_4000, Ez0_wo_4000, np.real(f(20, 25, Et0_wo_4000)), np.imag(f(20, 25, Et0_wo_4000)), z_func(Et0_wo_4000),  Et0_wo_4000), label=r'$n_1=4\cdot 10^3$')
    ax[0][1].plot(Et0_wo_8000, relative_error(Ex0_wo_8000, Ey0_wo_8000, Ez0_wo_8000, np.real(f(20, 25, Et0_wo_8000)), np.imag(f(20, 25, Et0_wo_8000)), z_func(Et0_wo_8000),  Et0_wo_8000), label=r'$n_2=8\cdot 10^3$')
    ax[0][1].plot(Et0_wo_16000, relative_error(Ex0_wo_16000, Ey0_wo_16000, Ez0_wo_16000, np.real(f(20, 25, Et0_wo_16000)), np.imag(f(20, 25, Et0_wo_16000)), z_func(Et0_wo_16000),  Et0_wo_16000), label=r'$n_3=1.6\cdot 10^4$')
    ax[0][1].plot(Et0_wo_32000, relative_error(Ex0_wo_32000, Ey0_wo_32000, Ez0_wo_32000, np.real(f(20, 25, Et0_wo_32000)), np.imag(f(20, 25, Et0_wo_32000)), z_func(Et0_wo_32000),  Et0_wo_32000), label=r'$n_4=3.2\cdot 10^4$')

    ax[1][0].set_ylim(10**(-15), 10**(1))
    ax[1][0].set_ylabel(r"Relative error [$\log_{10}(r_{err})$]")
    ax[1][0].set_xlabel(r"Time [$\mu$s]")
    ax[1][0].semilogy(t0_wo_4000, relative_error(x0_wo_4000, y0_wo_4000, z0_wo_4000, np.real(f(20, 25, t0_wo_4000)), np.imag(f(20, 25, t0_wo_4000)), z_func(t0_wo_4000),  t0_wo_4000), label='n=4000')
    ax[1][0].semilogy(t0_wo_8000, relative_error(x0_wo_8000, y0_wo_8000, z0_wo_8000, np.real(f(20, 25, t0_wo_8000)), np.imag(f(20, 25, t0_wo_8000)), z_func(t0_wo_8000),  t0_wo_8000), label='n=8000')
    ax[1][0].semilogy(t0_wo_16000, relative_error(x0_wo_16000, y0_wo_16000, z0_wo_16000, np.real(f(20, 25, t0_wo_16000)), np.imag(f(20, 25, t0_wo_16000)), z_func(t0_wo_16000),  t0_wo_16000), label='n=16000')
    ax[1][0].semilogy(t0_wo_32000, relative_error(x0_wo_32000, y0_wo_32000, z0_wo_32000, np.real(f(20, 25, t0_wo_32000)), np.imag(f(20, 25, t0_wo_32000)), z_func(t0_wo_32000),  t0_wo_32000), label='n=32000')

    ax[1][1].set_ylim(10**(-6), 10**(1))
    ax[1][1].set_xlabel(r"Time [$\mu$s]")
    ax[1][1].semilogy(Et0_wo_4000, relative_error(Ex0_wo_4000, Ey0_wo_4000, Ez0_wo_4000, np.real(f(20, 25, Et0_wo_4000)), np.imag(f(20, 25, Et0_wo_4000)), z_func(Et0_wo_4000),  Et0_wo_4000), label='n=4000')
    ax[1][1].semilogy(Et0_wo_8000, relative_error(Ex0_wo_8000, Ey0_wo_8000, Ez0_wo_8000, np.real(f(20, 25, Et0_wo_8000)), np.imag(f(20, 25, Et0_wo_8000)), z_func(Et0_wo_8000),  Et0_wo_8000), label='n=8000')
    ax[1][1].semilogy(Et0_wo_16000, relative_error(Ex0_wo_16000, Ey0_wo_16000, Ez0_wo_16000, np.real(f(20, 25, Et0_wo_16000)), np.imag(f(20, 25, Et0_wo_16000)), z_func(Et0_wo_16000),  Et0_wo_16000), label='n=16000')
    ax[1][1].semilogy(Et0_wo_32000, relative_error(Ex0_wo_32000, Ey0_wo_32000, Ez0_wo_32000, np.real(f(20, 25, Et0_wo_32000)), np.imag(f(20, 25, Et0_wo_32000)), z_func(Et0_wo_32000),  Et0_wo_32000), label='n=32000')

    ax[0][0].legend(loc='upper left')
    ax[0][1].legend(loc='upper left')
    plt.show()


n = [4000, 8000, 16000, 32000]
delta_max = np.zeros(len(n))
delta_max_E = np.zeros(len(n))
h = np.zeros(len(n))
t = [t0_wo_4000, t0_wo_8000, t0_wo_16000, t0_wo_32000]
x_RK4 = [x0_wo_4000, x0_wo_8000, x0_wo_16000, x0_wo_32000]
y_RK4 = [y0_wo_4000, y0_wo_8000, y0_wo_16000, y0_wo_32000]
z_RK4 = [z0_wo_4000, z0_wo_8000, z0_wo_16000, z0_wo_32000]

x_E = [Ex0_wo_4000, Ex0_wo_8000, Ex0_wo_16000, Ex0_wo_32000]
y_E = [Ey0_wo_4000, Ey0_wo_8000, Ey0_wo_16000, Ey0_wo_32000]
z_E = [Ez0_wo_4000, Ez0_wo_8000, Ez0_wo_16000, Ez0_wo_32000]

for k in range(len(n)):
    h[k] = 50/n[k]
    norms = np.zeros(len(x_RK4[k]))
    norms_E = np.zeros(len(x_E[k]))
    r_numerical = np.array([x_RK4[k], y_RK4[k], z_RK4[k]])
    r_numerical_E = np.array([x_E[k], y_E[k], z_E[k]])
    r_analytical = [np.real(f(20, 25, t[k])), np.imag(f(20, 25, t[k])),  z_func(t[k])]
    for i in range(len(r_numerical[:])):
        norms[:] = np.linalg.norm(r_analytical[:][i] - r_numerical[:][i])
        norms_E[:] = np.linalg.norm(r_analytical[:][i] - r_numerical_E[:][i])
    delta_max[k] = np.amax(norms)
    delta_max_E[k] = np.amax(norms_E)

if str(input("Maximum errors and rates of convergence? (y/n): ")) == "y":
    sum = np.zeros(len(n))
    sum_E = np.zeros(len(n))
    print(' ')
    for k in range(1, len(n)):
        sum[k] += (1/3)*np.log(delta_max[k]/delta_max[k-1])/np.log(h[k]/h[k-1])
        sum_E[k] += (1/3)*np.log(delta_max_E[k]/delta_max_E[k-1])/np.log(h[k]/h[k-1])
        print(f'Maximum error for RK4 n = {n[k]}: {(delta_max[k]):.4e}, and for forward Euler: {(delta_max_E[k]):.4e}')
    print(f'The rates of convergence for RK4 is {np.sum(sum):.2f}, and forward Euler: {np.sum(sum_E):.2f}')
    print(' ')


def readfile_2(filename):
    omega = []
    trapped_particles = []
    with open(filename) as infile:
        for line in infile:
            columns = line.split()
            omega.append(float(columns[0]))
            trapped_particles.append(float(columns[1]))
    infile.close()
    return np.array(omega), np.array(trapped_particles)



f_values = [0.100000, 0.400000, 0.700000]
colors = ['blue', 'red', 'green']
if str(input("Particles in trap for differnt amplitued? (y/n): "))== "y":
    fig, ax = plt.subplots(3, 1, figsize=(12, 7))
    ax[0].set_title(fr'Particles in the trap for different amplitudes ($f$)')
    ax[0].set_ylabel('Fractional percentage')
    ax[1].set_ylabel('Fractional percentage')
    ax[2].set_ylabel('Fractional percentage')
    plt.xlabel(r'Frequency $\omega$ [MHz]')
    ax[0].set_ylim(-5, 110)
    ax[1].set_ylim(-5, 50)
    ax[2].set_ylim(-5, 50)
    for i in range(len(f_values)):
        omega, trapped_particles = readfile_2('/textfiles/' + str(f_values[i]) + '00000.txt')
        omega_full, trapped_particles_full = readfile_2(str(f_values[i]) + '00000_with_hyades_full.txt')
        omega_resolved, trapped_particles_resolved = readfile_2(str(f_values[i]) + '00000_with_hyades_full_resolved.txt')
        ax[0].plot(omega, trapped_particles, label=f'f = {str((f_values[i]))}', color=colors[i], alpha=0.4)
        ax[1].plot(omega_full, trapped_particles_full, label=f'f = {str((f_values[i]))}', color=colors[i], alpha=0.4)
        xbot = np.linspace(omega_resolved[0], omega_resolved[-1], len(omega))
        ybot = np.linspace(-1, -1, len(omega))
        xleft = np.linspace(omega_resolved[0], omega_resolved[0], len(omega))
        yleft = np.linspace(-1, 45, len(omega))
        ytop = np.linspace(45, 45, len(omega))
        xright = np.linspace(omega_resolved[-1], omega_resolved[-1], len(omega))


        ax[1].plot(xbot, ybot, color='red')
        ax[1].plot(xleft, yleft, color='red')
        ax[1].plot(xbot, ytop, color='red')
        ax[1].plot(xright, yleft, color='red')
        ax[2].plot(omega_resolved, trapped_particles_resolved, label=f'f = {str((f_values[i]))}', color=colors[i], alpha=0.4)
    ax[0].legend(loc='lower left')
    ax[1].legend(loc='lower left')
    ax[2].legend(loc='lower left')
    plt.show()
