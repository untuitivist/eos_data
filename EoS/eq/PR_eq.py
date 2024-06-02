import numpy as np
import pandas as pd
import scipy as sc
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

class PR:
    def __init__(self, y_list, V_c_list, omega_list, k_PR_matrix, T_r_list, P_r_list, T_c_list, P_c_list, T_list, P_list, debug, is_column_r, is_index_r) -> None:
        self.R = 8.314472 * (10 ** 6)
        self.T_r_list = T_r_list
        self.P_r_list = P_r_list
        self.T_c_list = np.array(T_c_list)
        self.P_c_list = np.array(P_c_list)
        self.V_c_list = V_c_list
        self.T_list = T_list
        self.P_list = P_list
        self.omega_list = omega_list
        self.y_list = y_list
        self.k_matrix = k_PR_matrix
        self.debug = debug
        self.is_column_r = is_column_r
        self.is_index_r = is_index_r

    def a_list(self):
        P_c_list = np.array(self.P_c_list)
        T_c_list = np.array(self.T_c_list)
        a_list = 0.45724 * (self.R ** 2) * ((T_c_list ** 2) / P_c_list)
        if self.debug:
            print('a_list:', len(a_list))
            print(a_list)
        return a_list

    def a_mix(self) -> float:
        a_list = self.a_list()
        k_matrix = np.array(self.k_matrix)
        y_list = np.matrix(self.y_list)
        y_matrix = np.array(np.dot(y_list.T, y_list))
        a_list = np.matrix(a_list)
        a_matrix = np.sqrt(np.array(np.dot(a_list.T, a_list)))
        a_mix = np.sum(np.sum(-1 * y_matrix * a_matrix * (k_matrix - 1)))
        if self.debug:
            print('a_mix:', a_mix)
        return float(a_mix)

    def b_list(self):
        P_c_list = np.array(self.P_c_list)
        T_c_list = np.array(self.T_c_list)
        b_list = 0.0778 * self.R * (T_c_list / P_c_list)
        if self.debug:
            print('b_list:', len(b_list))
            print(b_list)
        return b_list

    def b_mix(self) -> float:
        b_list = self.b_list()
        b_mix = float(np.dot(np.matrix(self.y_list), np.matrix(b_list).T))
        if self.debug:
            print('b_mix:', b_mix)
        return b_mix

    def alpha_T_r(self, omega, T_r):
        alpha_T_r = (1 + ((0.37464 + 1.54226 * omega - 0.2699 * (omega ** 2)) * (1-np.sqrt(T_r)))) ** 2
        return alpha_T_r
    
    def Theta_list(self):
        a_mix = self.a_mix()
        omega = sum(np.array(self.y_list)*np.array(self.omega_list))
        T_r_list = self.T_r_list
        Thate_list = []
        for i in range(len(T_r_list)):
            alpha_T_r = self.alpha_T_r(omega, T_r_list[i])
            # print(a_mix, omega, T_r, alpha_T_r, a_mix * alpha_T_r)
            Thate_list.append(a_mix * alpha_T_r)
        if self.debug:
            print('Thate_list:', len(Thate_list))
            print(Thate_list)
        return Thate_list
    
    def varepsilon(self):
        b_mix = self.b_mix()
        varepsilon = -1 * (b_mix ** 2)
        if self.debug:
            print('varepsilon:', varepsilon)
        return varepsilon
    def delta(self):
        b_mix = self.b_mix()
        delta = 2 * b_mix
        if self.debug:
            print('delta:', delta)
        return delta
    def eta(self):
        b_mix = self.b_mix()
        eta = b_mix
        if self.debug:
            print('eta:', eta)
        return eta

    def P_equation(self, V, P, T, b, Theta, varepsilon, delta, eta) -> float:
        return float(P - ((self.R * T) / (V - b) - (Theta * (V - eta)) / ((V - b) * (V ** 2 + delta * V + varepsilon))))

    def V_DF(self) -> pd.DataFrame:
        b_mix = self.b_mix()
        Theta_list = self.Theta_list()
        V_matrix = []
        varepsilon = self.varepsilon()
        delta = self.delta()
        eta = self.eta()
        num = 1
        for T in self.T_list:
            V_t_list = []
            for P in self.P_list:
                # 设置合适的初始猜测值，例如初始猜测值为 V_c_list 中的某个值
                initial_guess = self.V_c_list[0]
                # 调用 scipy.optimize.fsolve 函数求解方程组
                V = sc.optimize.fsolve(self.P_equation, initial_guess, maxfev= 400, args=(P, T, b_mix, Theta_list[list(self.T_list).index(T)], varepsilon, delta, eta))[0]
                V_t_list.append(V)
                num+=1
            V_matrix.append(V_t_list)
        V_DataFrame = pd.DataFrame(np.array(V_matrix), index=self.T_r_list if self.is_index_r else self.T_list, columns=self.P_r_list if self.is_column_r else self.P_list)
        return V_DataFrame

    def Z(self, V, T, b_mix, Theta, varepsilon, delta, eta) -> float:
        return float(V / (V - b_mix) - \
                ((Theta / (self.R * T)) * V * (V - eta)) / ((V - b_mix) * (V ** 2 + delta * V + varepsilon)))

    def Z_DF(self) -> pd.DataFrame:
        Theta_list = self.Theta_list()
        Z_matrix = []
        b_mix = self.b_mix()
        V_DataFrame = np.matrix(self.V_DF())
        varepsilon = self.varepsilon()
        delta = self.delta()
        eta = self.eta()
        for t in range(self.T_list.__len__()):
            Z_t_list = []
            T = self.T_list[t]
            Theta = Theta_list[t]
            for p in range(self.P_list.__len__()):
                V = V_DataFrame[t, p]
                Z_t_list.append(self.Z(V, T, b_mix, Theta, varepsilon, delta, eta))
            Z_matrix.append(Z_t_list)
        return pd.DataFrame(Z_matrix, index=self.T_r_list if self.is_index_r else self.T_list, columns=self.P_r_list if self.is_column_r else self.P_list)
