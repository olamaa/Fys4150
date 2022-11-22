import numpy as np
import matplotlib.pyplot as plt

def readfile(filename):
    monte_carlo_cycles = []
    avg_eps = []
    avg_abs_m = []
    hist_eps = []
    with open(filename) as infile:
        infile.readline()
        for line in infile:
            columns = line.split()
            monte_carlo_cycles.append(float(columns[0]))
            avg_eps.append(float(columns[1]))
            avg_abs_m.append(float(columns[2]))
            hist_eps.append(float(columns[3]))
    infile.close()
    return np.array(monte_carlo_cycles), np.array(avg_eps), np.array(avg_abs_m), np.array(hist_eps)

"""

EXERCISE 5:

"""

L_20_T_10_rand = 'L_20_T_1.000000_random.txt'
L_20_T_24_rand = 'L_20_T_2.400000_random.txt'
L_20_T_10_nrand = 'L_20_T_1.000000_not_random.txt'
L_20_T_24_nrand = 'L_20_T_2.400000_not_random.txt'

cycles_20_10_rand, avg_eps_20_10_rand, avg_abs_m_20_10_rand, eps_hist_20_10_rand = readfile(L_20_T_10_rand)
cycles_20_24_rand, avg_eps_20_24_rand, avg_abs_m_20_24_rand, eps_hist_20_24_rand = readfile(L_20_T_24_rand)
cycles_20_10_nrand, avg_eps_20_10_nrand, avg_abs_m_20_10_nrand, eps_hist_20_10_nrand = readfile(L_20_T_10_nrand)
cycles_20_24_nrand, avg_eps_20_24_nrand, avg_abs_m_20_24_nrand, eps_hist_20_24_nrand = readfile(L_20_T_24_nrand)

fig, ax = plt.subplots(2, 2, figsize=(10, 7))
ax[0,0].set_title('Random initial matrix')
ax[0,0].plot(cycles_20_10_rand/1000, avg_eps_20_10_rand, label=r'$<\epsilon>$, T=1')
ax[0,0].plot(cycles_20_24_rand/1000, avg_eps_20_24_rand, label=r'$<\epsilon>$, T=2.4')
ax[0,0].set_ylim(-2.1, 0.1)
ax[0,0].legend()
ax[1,0].set_title('Not random initial matrix')
ax[1,0].plot(cycles_20_10_nrand/1000, avg_eps_20_10_nrand, label=r'$<\epsilon>$, T=1')
ax[1,0].plot(cycles_20_24_nrand/1000, avg_eps_20_24_nrand, label=r'$<\epsilon>$, T=2.4')
ax[1,0].set_ylim(-2.1, 0.1)
ax[1,0].legend()
ax[0,1].set_title('Random initial matrix')
ax[0,1].plot(cycles_20_10_rand/1000, avg_abs_m_20_10_rand, label=r'$<|m|>$, T=1')
ax[0,1].plot(cycles_20_24_rand/1000, avg_abs_m_20_24_rand, label=r'$<|m|>$, T=2.4')
ax[0,1].set_ylim(-0.1, 1.1)
ax[0,1].legend()
ax[1,1].set_title('Not random initial matrix')
ax[1,1].plot(cycles_20_10_nrand/1000, avg_abs_m_20_10_nrand, label=r'$<|m|>$, T=1')
ax[1,1].plot(cycles_20_24_nrand/1000, avg_abs_m_20_24_nrand, label=r'$<|m|>$, T=2.4')
ax[1,1].set_ylim(-0.1, 1.1)
ax[1,1].legend()
ax[1,0].set_xlabel(r'Monte Carlo cycles $\times 10^3$')
ax[1,1].set_xlabel(r'Monte Carlo cycles $\times 10^3$')
ax[1,0].set_ylabel(r'$<\epsilon>$, <|m|>')
ax[0,0].set_ylabel(r'$<\epsilon>$, <|m|>')
plt.show()

"""

EXERCISE 6:

"""

after_burn_in = np.where(cycles_20_10_rand>15000)[0] #burning time ended after approx 15000 monte carlo cycles
fig, ax = plt.subplots(2, 2, figsize=(12, 7))
ax[0,0].hist(eps_hist_20_10_rand[after_burn_in],  bins=6,  label=r'$\epsilon_{T=1}$, random',  density=True, log=True)
ax[0,0].set_xlim(-2.01,-1.88)
ax[1,0].hist(eps_hist_20_10_nrand[after_burn_in], bins=6,  label=r'$\epsilon_{T=1}$, not random', density=True, log=True)
ax[1,0].set_xlim(-2.01,-1.88)
#ax[0,1].hist(eps_hist_20_24_rand[after_burn_in],  bins=105, label=r'$\epsilon_{T=2.4}$, random',  density=True, log=False)
#ax[1,1].hist(eps_hist_20_24_nrand[after_burn_in], bins=110, label=r'$\epsilon_{T=2.4}$, not random', density=True, log=False)

from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy.stats import norm

n, bins, patches = ax[0,1].hist(eps_hist_20_24_rand[after_burn_in],  bins=105, label=r'$\epsilon_{T=2.4}$, random',  density=True, log=False)
n2, bins2, patches2 = ax[1,1].hist(eps_hist_20_24_nrand[after_burn_in], bins=110, label=r'$\epsilon_{T=2.4}$, not random', density=True, log=False)
(mu, sigma) = norm.fit(eps_hist_20_24_rand[after_burn_in])
(mu2, sigma2) = norm.fit(eps_hist_20_24_nrand[after_burn_in])
y = norm.pdf(bins, mu, sigma)
y2 = norm.pdf(bins2, mu2, sigma2)

ax[0,1].plot(bins, y, 'r--', linewidth=2, label=rf'$\sigma^2=${(sigma**2):.2f}')
ax[1,1].plot(bins2, y2, 'r--', linewidth=2, label=rf'$\sigma^2=${(sigma2**2):.2f}')


ax[0,0].set_ylabel('PDF')
ax[1,0].set_ylabel('PDF')
ax[0,1].set_ylabel('PDF')
ax[1,1].set_ylabel('PDF')
ax[0,0].set_xlabel('Energy [J]')
ax[1,0].set_xlabel('Energy [J]')
ax[0,1].set_xlabel('Energy [J]')
ax[1,1].set_xlabel('Energy [J]')
ax[0,0].legend()
ax[0,1].legend()
ax[1,0].legend()
ax[1,1].legend()
plt.show()

"""

EXERCISE 8:

"""

def readfile_2(filename):
    T   = []
    L   = []
    eps = []
    m   = []
    C_v = []
    chi = []
    with open(filename) as infile:
        infile.readline()
        for line in infile:
            columns = line.split()
            T.append(float(columns[0]))
            L.append(float(columns[1]))
            eps.append(float(columns[2]))
            m.append(float(columns[3]))
            C_v.append(float(columns[4]))
            chi.append(float(columns[5]))
    infile.close()
    return np.array(T), np.array(L), np.array(eps), np.array(m), np.array(C_v), np.array(chi)

fig, ax = plt.subplots(2, 2, sharex=True, figsize=(14, 7))
size = 15
def plotter(filename):
    T, L, eps, m, C_v, chi = readfile_2(filename)
    names = ['T', 'eps', 'm', 'C_v', 'chi']
    arrays = [T, eps, m, C_v, chi]
    T_new = np.zeros((len(arrays), int(len(T)/4 -3)))
    eps_new = np.zeros((len(arrays), int(len(T)/4 -3)))
    m_new = np.zeros((len(arrays), int(len(T)/4 -3)))
    C_v_new = np.zeros((len(arrays), int(len(T)/4 -3)))
    chi_new = np.zeros((len(arrays), int(len(T)/4 -3)))
    new_array = np.zeros((len(arrays), int(len(T)/4 -3)))

    for array in range(len(arrays)):
        for i in range(4): #number of lattices we save arrays for (40, 60, 80, 100, 200)
            new_array[i] = arrays[array][np.where(L==int(40 + 20*i))]
            if names[array] == 'T':
                T_new[i,:] =  new_array[i,:]
            elif names[array] == 'eps':
                eps_new[i,:] =  new_array[i,:]
            elif names[array] == 'm':
                m_new[i,:] =  new_array[i,:]
            elif names[array] == 'C_v':
                C_v_new[i,:] =  new_array[i,:]
            elif names[array] == 'chi':
                chi_new[i,:] =  new_array[i,:]

    arrays = [eps_new, m_new, C_v_new, chi_new]
    colors = ['blue', 'green', 'red', 'orange']
    for array in range(len(arrays)+1):
        for i in range(len(arrays)): #number of lattices we save arrays for (40, 60, 80, 100, 1000)
            if names[array] == 'eps':
                ax[0,0].set_title(r'$<\epsilon>$')
                ax[0,0].scatter(T_new[i], eps_new[i], label=rf'L={int(40 + 20*i)}', color=colors[i], s=size)
            elif names[array] == 'm':
                ax[1,0].set_title(r'$<|m|>$')
                ax[1,0].scatter(T_new[i], m_new[i], label=rf'L={int(40 + 20*i)}', color=colors[i], s=size)
            elif names[array] == 'C_v':
                ax[0,1].set_title(r'C$_v$')
                ax[0,1].scatter(T_new[i], C_v_new[i], label=rf'L={int(40 + 20*i)}', color=colors[i], s=size)
            elif names[array] == 'chi':
                ax[1,1].set_title(r'$\chi$')
                ax[1,1].scatter(T_new[i], chi_new[i], label=rf'L={int(40 + 20*i)}', color=colors[i], s=size)

    ax[0,0].set_ylabel('Energy')
    ax[1,0].set_ylabel('Net magnetization')
    ax[0,1].set_ylabel(r'Specific heat capacity')
    ax[1,1].set_ylabel(r'Susceptibility')

    ax[1,1].set_xlabel(r'T [J/k$_B$]')
    ax[1,0].set_xlabel(r'T [J/k$_B$]')
    ax[0,0].legend()
    ax[0,1].legend()
    ax[1,0].legend()
    ax[1,1].legend()
    return arrays, T_new, C_v_new

filename = '8a_beehive.txt'
arrays, T_new, C_v_new = plotter(filename)
plt.show()

"""

EXERCISE 9:

"""

plt.figure()
plt.scatter(T_new[3][:], C_v_new[3][:], alpha=0.1, color='orange')
plt.scatter(T_new[2][:], C_v_new[2][:], alpha=0.1, color='red')
plt.scatter(T_new[1][:], C_v_new[1][:], alpha=0.1, color='green')
plt.scatter(T_new[0][:], C_v_new[0][:], alpha=0.1, color='blue')

exs0 = np.where(C_v_new[0][:] == np.amax(C_v_new[0][:]))[0][0]
exs1 = np.where(C_v_new[1][:] == np.max(C_v_new[1][:]))[0][0]
exs2 = np.where(C_v_new[2][:] == np.max(C_v_new[2][:]))[0][0]
exs3 = np.where(C_v_new[3][:] == np.max(C_v_new[3][:]))[0][0]

T_c = [T_new[3][exs3], T_new[2][exs2], T_new[1][exs1], T_new[0][exs0]]
L_new = [40, 60, 80, 100]

plt.scatter(T_new[3][exs3], C_v_new[3][exs3], color='orange', label='L=100')
plt.scatter(T_new[2][exs2], C_v_new[2][exs2], color='red', label='L=80')
plt.scatter(T_new[1][exs1], C_v_new[1][exs1], color='green', label='L=60')
plt.scatter(T_new[0][exs0], C_v_new[0][exs0], color='blue', label='L=40')
plt.ylabel(r'Specific heat capacity')
plt.xlabel(r'T [J/k$_B$]')
plt.legend()
plt.show()

a = []
c = []
for i in range(len(L_new)-1):
    a.append((T_c[i] - T_c[i+1]) / (L_new[i]**(-1) - L_new[i+1]**(-1)))
a = np.average(a)
for i in range(len(arrays)):
    c.append(T_c[i] - a*L_new[i]**(-1))
c = np.average(c)
print(f'{c} [J/k_B], {np.abs(100*(c-2.269)/2.269):.2f} % relative error from Onsagers 2.269')
