import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import scipy.stats as scst


class markov_chain_FORT:
    """
    this class defines a markov chain used to describe the dynamics of the number of atoms in a dipolar trap.
    """

    def __init__(self, R, gamma, beta, truncation_size, method="model_1") -> None:
        self.R = R
        self.gamma = gamma
        self.beta = beta
        self.truncation_size = truncation_size
        self.method = method
        self.transition_matrix = self.get_transition_matrix(method)
        self.evos = []
        pass

    def get_transition_matrix(self, method):
        if method == "model_1":
            return self.get_transition_matrix_model_1()
        if method == "model_2":
            return self.get_transition_matrix_model_2()
        print("method must be model_1 or model_2")
        return
    
    def get_transition_matrix_model_1(self):
        P = np.zeros([self.truncation_size, self.truncation_size])
        P[0, 0] = - self.R
        P[0, 1] = self.R
        for i in range(1, self.truncation_size-1):
            p_i = i*(self.gamma+self.beta*(i-1))
            P[i, i+1] = self.R
            P[i, i] = - self.R - p_i
            P[i, i-1] = p_i
        P[self.truncation_size-1, self.truncation_size-1] = - (self.truncation_size-1)*(self.gamma+self.beta*((self.truncation_size-1)-1)) - self.R
        P[self.truncation_size-1, self.truncation_size-2] = (self.truncation_size-1)*(self.gamma+self.beta*((self.truncation_size-1)-1))
        return P
    
    def get_transition_matrix_model_2(self):
        P = np.zeros([self.truncation_size, self.truncation_size])
        P[0, 0] = - self.R
        P[0, 1] = self.R
        P[1, 0] = self.gamma
        P[1, 1] = - self.R - self.gamma
        P[1, 2] = self.R
        for i in range(2, self.truncation_size-1):
            p_i_1 = i*self.gamma
            p_i_2 = self.beta/2*(i-1)*i
            P[i, i+1] = self.R
            P[i, i] = - self.R - p_i_1 - p_i_2
            P[i, i-1] = p_i_1
            P[i, i-2] = p_i_2
        P[self.truncation_size-1, self.truncation_size-1] = - (self.truncation_size-1)*(self.gamma+self.beta/2*((self.truncation_size-1)-1)) - self.R
        P[self.truncation_size-1, self.truncation_size-2] = (self.truncation_size-1)*self.gamma
        P[self.truncation_size-1, self.truncation_size-3] = (self.truncation_size-1)*self.beta/2*((self.truncation_size-1)-1)
        return P

    def function_(self, t, v):
        return np.dot(v, self.transition_matrix)

    def function_BURST(self, t, v):
        if self.FIRST == True and t > self.burst_time:
            self.FIRST = False
            self.R = 0
            #self.beta = 0
            self.transition_matrix = self.get_transition_matrix(self.method)
        return np.dot(v, self.transition_matrix)

    def evolve(self, v0, time, steps, burst_time = None, scale="linear"):
        if scale == "log":
            ts = 10**np.linspace(time[0], time[1], steps)
        else:
            ts = np.linspace(0, time, steps)

        if burst_time is None:
            self.evo = solve_ivp(self.function_, t_span=(ts[0],ts[-1]), y0=v0, t_eval=ts)
            return
        else:
            aus_R = self.R
            aus_beta = self.beta
            self.FIRST = True
            self.burst_time = burst_time
            self.evo = solve_ivp(self.function_BURST, t_span=(ts[0],ts[-1]), y0=v0, t_eval=ts)
            self.R = aus_R
            self.beta = aus_beta
            return
        """
    def evolve(self, v0, time, steps, pars=None, scale="linear"):
         
        evolve the distribution, from a initial probability vector "v0",
        for a time "time", and compuding "steps" point 
        
        if scale == "log":
            ts = 10**np.linspace(time[0], time[1], steps)
        else:
            ts = np.linspace(0, time, steps)

        if pars is None:
            self.evos.append(solve_ivp(self.function_, t_span=(ts[0],ts[-1]), y0=v0, t_eval=ts))
            return self.evos[-1]
        else:
            aus_R = self.R
            aus_beta = self.beta
            aus_gamma = self.gamma
            aus_transition_matrix = self.transition_matrix 
            self.R = pars[0]
            self.beta = pars[1]
            self.gamma = pars[2]
            self.transition_matrix = self.get_transition_matrix(self.method)
            self.evos.append(solve_ivp(self.function_, t_span=(ts[0],ts[-1]), y0=v0, t_eval=ts))
            self.R = aus_R
            self.beta = aus_beta
            self.gamma = aus_gamma
            self.transition_matrix = aus_transition_matrix
            return self.evos[-1]
        """
    def show_matrix(self):
        from sympy import Matrix
        return Matrix(self.transition_matrix)
    
    def plot_mandel_Q(self, scale="lin"):
        means_ = []
        vars_ = []
        Qs_ = []
        N_nm1s_ = []
        N_2s_ = []
        for i in range(len(self.evo.y[0])):
            mean_, var_, mandel_Q_, N_nm1_, N_2_ = mean_var_mandel_Q(self.evo.y[:,i])
            means_.append(mean_)
            vars_.append(var_)
            Qs_.append(mandel_Q_)
            N_nm1s_.append(N_nm1_)
            N_2s_.append(N_2_)

        def fun_paper(t, N):
            return self.R - self.gamma*N - self.beta*N*(N-1)

        sol_paper = solve_ivp(fun=fun_paper, t_span=(0, self.evo.t[-1]), y0=[0], t_eval=np.linspace(0, self.evo.t[-1], len(self.evo.t))) 

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        #print(sol_paper.y[-1], sol_paper.t[-1])
        #print(means_[-1], self.evo.t[-1])

        ax1.plot(self.evo.t, Qs_, linestyle="--", label="mandel Q", color="blue")
        ax1.set_ylabel("mandel Q", fontsize=15)
        
        ax1.set_xlabel(r"time [s]", fontsize=15)
        #ax2.plot([0, self.evo.t[-1]],[1,1],c="gray",linestyle="--",alpha=0.7, label=r"$\langle N \rangle=1$")
        ax2.plot(self.evo.t, means_, color="red",label=r"$\langle N \rangle$")
        #ax2.plot(self.evo.t, N_nm1s_, color="black",label=r"$\langle N \rangle \langle (N-1) \rangle$")
        #ax2.plot(self.evo.t, N_2s_, color="orange",label=r"$\langle N(N-1) \rangle$")
        #ax2.plot(self.evo.t, sol_paper.y[0], label="paper")
        #vars_ = np.array(vars_)#*self.beta
        #ax2.plot(self.evo.t, vars_, color="blue", label="var")
        #ax2.plot(self.evo.t, (sol_paper.y[0]-means_)/self.beta, label="diff")

        ax2.set_ylabel(r"$\langle N \rangle$", fontsize=15)
        if scale == "log":
            ax1.set_xscale("log")
            ax2.set_xscale("log")

        fig.legend(loc=(0.15,0.7), fontsize=15)

        plt.show()
        return
"""
    def plot_mandel_Q(self, evo_y, t_range, scale="lin"):
        means_ = []
        vars_ = []
        Qs_ = []
        N_nm1s_ = []
        N_2s_ = []
        for i in range(len(evo_y[0])):
            mean_, var_, mandel_Q_, N_nm1_, N_2_ = mean_var_mandel_Q(evo_y[:,i])
            means_.append(mean_)
            vars_.append(var_)
            Qs_.append(mandel_Q_)
            N_nm1s_.append(N_nm1_)
            N_2s_.append(N_2_)

        def fun_paper(t, N):
            return self.R - self.gamma*N - self.beta*N*(N-1)

        sol_paper = solve_ivp(fun=fun_paper, t_span=(t_range[0], t_range[1]), y0=[0], t_eval=np.linspace(t_range[0], t_range[1], len(evo_y[0]))) 

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        #print(sol_paper.y[-1], sol_paper.t[-1])
        #print(means_[-1], self.evo.t[-1])
        ts = np.linspace(t_range[0], t_range[1], len(evo_y[0]))
        ax1.plot(ts, Qs_, linestyle="--", label="mandel Q", color="blue")
        ax1.set_ylabel("mandel Q", fontsize=15)
        
        ax1.set_xlabel(r"time [s]", fontsize=15)
        #ax2.plot([0, self.evo.t[-1]],[1,1],c="gray",linestyle="--",alpha=0.7, label=r"$\langle N \rangle=1$")
        ax2.plot(ts, means_, color="red",label=r"$\langle N \rangle$")
        #ax2.plot(self.evo.t, N_nm1s_, color="black",label=r"$\langle N \rangle \langle (N-1) \rangle$")
        #ax2.plot(self.evo.t, N_2s_, color="orange",label=r"$\langle N(N-1) \rangle$")
        #ax2.plot(self.evo.t, sol_paper.y[0], label="paper")
        #vars_ = np.array(vars_)#*self.beta
        #ax2.plot(self.evo.t, vars_, color="blue", label="var")
        #ax2.plot(self.evo.t, (sol_paper.y[0]-means_)/self.beta, label="diff")

        ax2.set_ylabel(r"$\langle N \rangle$", fontsize=15)
        if scale == "log":
            ax1.set_xscale("log")
            ax2.set_xscale("log")

        fig.legend(loc=(0.15,0.7), fontsize=15)

        plt.show()
        return
"""
def mean_var_mandel_Q(probability_vec):
    mean = 0
    men_2 = 0
    for i in range(len(probability_vec)):
        mean += i*probability_vec[i]
        men_2 += i**2*probability_vec[i]
    variance = men_2 - mean**2
    Q = variance/mean - 1
    N_nm1 = mean**2 - mean
    N_nm12 = men_2 - mean
    return mean, variance, Q, N_nm1, N_nm12