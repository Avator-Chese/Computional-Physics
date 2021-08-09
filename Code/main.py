import numpy as np
import matplotlib.pyplot as plt
from random import randint, random
#==============================================================

N = 20
steps = 10 * N
N_T = steps
J = 1
T = np.linspace(0.01, 10, N_T)


def initial_state(start='Low'):
    if start == 'High':
        state = np.zeros(N)
        for i in range(N):
            state[i] = 2 * randint(0, 1) - 1
    if start == 'Low':
        state = np.full(N, 1)
    return state


def cal_tol_energy(some_con, J):
    energy = 0
    for i in range(len(some_con)):
        spin = some_con[i]
        neighbor_spin = some_con[(i + 1) % N] + some_con[(i - 1) % N]
        energy += - J * spin * neighbor_spin
    return energy / 2


def cal_tol_magnetization(some_con):
    return np.sum(some_con)


def con_update(some_con, beta, J):
    for i in range(N):
        a = randint(0, N - 1)
        spin = some_con[a]

        neighbor_spin = some_con[(a + 1) % N] + some_con[(a - 1) % N]
        dE = 2 * J * (spin * neighbor_spin)

        if dE <= 0:
            some_con[a] *= - 1
        else:
            if np.exp(-dE * beta) >= random():
                some_con[a] *= -1
    return some_con


def collect_spin_data(some_con, beta, J):
    spin_record = []
    for i in range(steps):
        spin_record.append(np.array(some_con))
        some_con = con_update(some_con, beta, J)
    spin_record = np.transpose(spin_record)
    return spin_record


def get_lattice_plot(spin_data):
    x = np.arange(steps)
    y = np.arange(N)
    X, Y = np.meshgrid(x, y)
    colors = np.zeros(N * steps)
    for i in range(N):
        for j in range(steps):
            if spin_data[i][j] == 1:
                colors[i * steps + j] = 0
            elif spin_data[i][j] == -1:
                colors[i * steps + j] = 1
    plt.figure()
    plt.scatter(X, Y, c=colors, s=13, marker='.')
    plt.xlabel("Steps", fontsize=15)
    plt.ylabel("Position", fontsize=15)
    plt.title('Configuration vs Steps')
    # plt.show()


con = initial_state('Low')
data1 = collect_spin_data(con, 1 / 1, J)
get_lattice_plot(data1)

# T_ch = 1
# M_single_spin = []
# con = initial_state('High')
# for i in range(steps):
#     M_single_spin.append(cal_tol_magnetization(con) / N)
#     con_update(con, 1 / T_ch, 1)

# plt.figure()
# plt.plot(np.arange(steps), M_single_spin)
# plt.xlabel('Steps',fontsize=15)
# plt.ylabel(" Average Magnetization per Spin")
# plt.grid(linestyle='-.')
# plt.title(' Average Magnetization every Time Step')


E, M, C, X = np.zeros(N_T), np.zeros(N_T), np.zeros(N_T), np.zeros(N_T)
n_1, n_2 = (1.0) / (steps * N), (1.0) / (steps * steps * N)

for i in range(N_T):
    E_tol = M_tol = E_tol_2 = M_tol_2 = 0
    con = initial_state('High')
    T_i = 1.0 / T[i]
    T_i2 = T_i * T_i

    for j in range(steps):
        con_update(con, T_i, J)
        en = cal_tol_energy(con, J)
        ma = cal_tol_magnetization(con)

        M_tol = M_tol + ma
        E_tol = E_tol + en
        M_tol_2 = M_tol_2 + ma**2
        E_tol_2 = E_tol_2 + en**2

    E[i] = n_1 * E_tol  # Energy
    M[i] = n_1 * M_tol  # Magnetization
    C[i] = (n_1 * E_tol_2 - n_2 * E_tol * E_tol) * T_i2  # Specific Heat
    X[i] = (n_1 * M_tol_2 - n_2 * M_tol * M_tol) * T_i  # Susceptibility


plt.figure()
plt.plot(T, E, 'r')
plt.scatter(T, E, s=50, marker='*', color='g')
plt.xlabel("Temperature ($J/k_B$)", fontsize=15)
plt.ylabel(" Average Energy per Spin($J$)", fontsize=15)
plt.title('Average Energy vs Temperature')
plt.grid(linestyle='-.')

plt.figure()
plt.plot(T, M, 'r')
plt.scatter(T, M, s=50, marker='*', color='g')
plt.xlabel("Temperature ($J/k_B$)", fontsize=15)
plt.ylabel(" Average Magnetization per Spin", fontsize=15)
plt.title('Average Magnetization vs Temperature')
plt.grid(linestyle='-.')


plt.figure()
plt.scatter(T, C, s=50, marker='*', color='g')
plt.xlabel("Temperature ($J/k_B$)", fontsize=15)
plt.ylabel(" Specific Heat ", fontsize=15)
plt.title('Specific Heat vs Temperature')
plt.grid(linestyle='-.')

plt.figure()
plt.plot(T, X, 'r')
plt.scatter(T, X, s=50, marker='*', color='g')
plt.xlabel("Temperature ($J/k_B$)", fontsize=15)
plt.ylabel(" Magnetic Susceptibility", fontsize=15)
plt.title('Magnetic Susceptibility vs Temperature')
plt.grid(linestyle='-.')


plt.show()
