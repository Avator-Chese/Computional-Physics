import numpy as np
import matplotlib.pyplot as plt
from random import randint, random
#==============================================================

N=20
steps = 100
N_T=100

def initial_state(start='Low'):
    if start=='High':
        state=np.zeros(N)
        for i in range(N):
            state[i]=2*randint(0,1)-1
    if start=='Low':
        state=np.full(N,1)
    return state

def cal_tol_energy(some_con,J):
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
        a=randint(0,N-1)
        spin=some_con[a]

        neighbor_spin=some_con[(a+1)%N]+some_con[(a-1)%N]
        dE=2*J*(spin*neighbor_spin)

        if dE<=0:
            some_con[a] *= - 1
        else:
            if np.exp(-dE*beta)>random():
                some_con[a] *= -1
    return some_con


def collect_spin_data(some_con,beta,J):
    spin_record=[]
    for i in range(steps):
        spin_record.append(np.array(some_con))
        some_con=con_update(some_con,beta,J)
    spin_record=np.transpose(spin_record)
    return spin_record

def get_lattice_plot(spin_data):
    x=np.arange(steps)
    y=np.arange(N)
    X,Y=np.meshgrid(x,y)
    colors=np.zeros(N*steps)
    for i in range(N):
        for j in range(steps):
            if spin_data[i][j]==1:
                colors[i*steps+j]=0
            elif spin_data[i][j]==-1:
                colors[i*steps+j]=1
    plt.scatter(X,Y,c=colors,s=13,marker='.')
    plt.show()


T=np.linspace(0.01,5,N_T)
E,M,C,X = np.zeros(N_T), np.zeros(N_T), np.zeros(N_T), np.zeros(N_T) 
n_1, n_2 = (1.0)/(steps*N), (1.0)/(steps*steps*N)

for i in range(N_T):
    E_1 = M_1 = E_2 = M_2 = 0
    con = initial_state('Low')
    i_T = 1.0 / T[i]
    i_T_2 = i_T * i_T

    for j in range(steps): 
        con_update(con, i_T, 1)
        Energy = cal_tol_energy(con, 1) 
        Magnet= cal_tol_magnetization(con) 

        M_1 = M_1 + Magnet
        E_1 = E_1 + Energy
        M_2 = M_2 + Magnet*Magnet
        E_2 = E_2 + Energy*Energy

    E[i] = n_1*E_1 #Energy
    M[i] = n_1*M_1 #Magnetization
    C[i] = (n_1*E_2 - n_2*E_1*E_1)*i_T_2 #Specific Heat
    X[i] = (n_1*M_2 - n_2*M_1*M_1)*i_T #Susceptibility


    
plt.plot(T, X)
plt.show()

