import qutip as qtp
import numpy as np
from arc import Rubidium87

class single_atom:
    def __init__(self, power_blue, waste_blue, power_IR, waste_IR, detuning_blue, power_depumping=None, waste_depumping=None, detuning_IR=0) -> None:
        """ 
        Here a single rubidium87 atom is defined.
        the considered levels are 
            - the ground state 5s(1/2)
            - the intermediate state 6p(3/2)
            - the rydberg state 70s(1/2)
        
        The hamiltonian is defined as follows:
                             0        Omega_b/2        0             0 
            H =         Omega_b/2    - \Delta      Omega_IR/2        0
                             0        Omega_IR/2    -delta       Omega_dep/2
                             0           0         Omega_dep/2     -delta 
        where the columns (rows) are respectivelly the following dressed states:
                    - |g, N_IR+1, N_b+1, N_d>   with energy   E_0 = 0.
                    - |i, N_IR+1, N_b, N_d>     with energy   -Delta_blue.
                    - |r, N_IR, N_b, N_d>       with energy   -delta_IR.
                    - |i, N_IR, N_b, N_d+1>     with energy   -delta_IR.
        """
        self.two_photon_rabi_freq = get_two_photon_rabi_freq(power_blue, waste_blue, power_IR, waste_IR, detuning_blue)
        self.lower_transition_rabi_freq = Rubidium87().getRabiFrequency(n1=5,  l1=0, j1=0.5, mj1=0.5,  n2=6,  l2=1, j2=3/2, q=-1, laserPower=power_blue, laserWaist=waste_blue, s=0.5)/(2*np.pi)
        self.higher_transition_rabi_freq = Rubidium87().getRabiFrequency(n1=70,  l1=0, j1=0.5, mj1=0.5,  n2=6,  l2=1, j2=3/2, q=-1, laserPower=power_IR, laserWaist=waste_IR, s=0.5)/(2*np.pi)
        if power_depumping is not None and waste_depumping is not None:
            self.depumping_rabi_freq = Rubidium87().getRabiFrequency(n1=70,  l1=0, j1=0.5, mj1=0.5,  n2=6,  l2=1, j2=3/2, q=-1, laserPower=power_depumping, laserWaist=waste_depumping, s=0.5)/(2*np.pi)
        else:
            self.depumping_rabi_freq = 0
        self.Gamma_6p = 1/Rubidium87().getStateLifetime(n=6, l=1, j=1.5, temperature=150*10**-6, includeLevelsUpTo=0, s=0.5)*10**6 # in MHz
        self.Gamma_70s = 1/Rubidium87().getStateLifetime(n=70, l=0, j=0.5, temperature=150*10**-6, includeLevelsUpTo=0, s=0.5)*10**6 # in MHz
        self.H = np.array([
            [0, self.lower_transition_rabi_freq/2, 0, 0],
            [self.lower_transition_rabi_freq/2, -detuning_blue, self.higher_transition_rabi_freq/2, 0],
            [0, self.higher_transition_rabi_freq/2, -detuning_IR, self.depumping_rabi_freq/2],
            [0, 0, self.depumping_rabi_freq/2, -detuning_IR]
        ])
        pass

    def evolve(self, time, steps, initial_DM, include_spontaneous_emission=True, include_defasing=False):
        times = np.linspace(0, time, steps)
        # jump operators
        def_ry = 0
        def_i = 0

        S_minus_i = qtp.Qobj(self.Gamma_6p**(1/2)*(np.array([[0,1,0,1],[0,0,0,0],[0,0,0,0],[0,0,0,0]])))
        S_minus_ry = qtp.Qobj(self.Gamma_70s**(1/2)*np.array([[0,0,1,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]))
        J_def_i = qtp.Qobj(def_i**(1/2)*np.array([[0,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,1]]))
        J_def_ry = qtp.Qobj(def_ry**(1/2)*np.array([[0,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,0]]))

        Jump_ops = [S_minus_i, S_minus_ry, J_def_ry, J_def_i]

        result = qtp.mesolve(qtp.Qobj(H), qtp.Qobj(psi0), times, c_ops = Jump_ops)

class two_atoms:
    def __init__(self) -> None:
        pass

def get_two_photon_rabi_freq(power_blue, waste_blue, power_IR, waste_IR, detuning_blue):

    atom=Rubidium87()
    # transition from 5s to 6p
    O5s_6p = atom.getRabiFrequency(n1=5,  l1=0, j1=0.5, mj1=0.5,  n2=6,  l2=1, j2=3/2, q=-1, laserPower=power_blue, laserWaist=waste_blue, s=0.5)/(2*np.pi)
    # transition from 6p to 70s
    O6p_70s = atom.getRabiFrequency(n1=70,  l1=0, j1=0.5, mj1=0.5,  n2=6,  l2=1, j2=3/2, q=-1, laserPower=power_IR, laserWaist=waste_IR, s=0.5)/(2*np.pi)

    return (O5s_6p*O6p_70s)/(detuning_blue)
