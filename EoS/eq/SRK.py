import numpy as np
import pandas as pd
import scipy as sc
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

class SRK:
    def __init__(self, y_dict, omega_list, k_matrix, T_r_list, P_r_list, T_c_list, P_c_list, T_list, P_list, debug, is_column_r, is_index_r) -> None:
        self.R = 8.314472 * (10 ** 6)
        self.omega_list = omega_list
        self.k_matrix = k_matrix
        self.T_r_list = T_r_list
        self.P_r_list = P_r_list
        self.T_c_list = np.array(T_c_list)
        self.P_c_list = np.array(P_c_list)
        self.T_list = T_list
        self.P_list = P_list
        self.gas_list = list(y_dict.keys())
        self.y_list = list(y_dict.values())
        self.debug = debug
        self.is_column_r = is_column_r
        self.is_index_r = is_index_r

    def alpha_T_r_matrix(self) -> np.ndarray:
        T_r = list(self.T_r_list)
        omega = list(self.omega_list)
        T_r_matrix = np.array([T_r] * len(omega))
        omega_matrix = np.array([omega] * len(T_r)).T
        T_r = T_r_matrix
        omega = omega_matrix
        if T_r.shape != omega.shape:
            raise ValueError('T_r_matrix和omega_matrix的维度不一致')
        if self.debug:
            print('T_r_matrix:', T_r.shape)
            print(T_r)
            print('omega_matrix:', omega.shape)
            print(omega)
        alpha_T_r_matrix = (1 + (1 - (T_r ** 0.5)) * (0.48508 + 1.55171 * omega - 0.15613 * (omega ** 2))) ** 2
        if 'H_2' in self.gas_list:
            i = self.gas_list.index('H_2')
            alpha_T_r_matrix[i] = 1.096 * np.exp(-0.15114 * self.T_r_list)
        if self.debug:
            print('alpha_T_r_matrix:', alpha_T_r_matrix.shape)
            print(alpha_T_r_matrix)
        return alpha_T_r_matrix

    def a_i_alpha_T_r_matrix(self) -> np.matrix:
        alpha_T_r_matrix = self.alpha_T_r_matrix()
        R = self.R
        T_c = self.T_c_list
        P_c = self.P_c_list
        a_c = np.array([0.42748 * ((R ** 2) * (T_c ** 2)) / P_c] * len(alpha_T_r_matrix[0])).T
        a_i_alpha_T_r_matrix = np.matrix(alpha_T_r_matrix * a_c)
        if self.debug:
            print('a_i_alpha_T_r_matrix:', a_i_alpha_T_r_matrix.shape)
            print(a_i_alpha_T_r_matrix)
        return a_i_alpha_T_r_matrix

    def a_mix_list(self) -> np.ndarray:
        y_list = np.matrix(self.y_list)
        y_ij = np.array(np.dot(y_list.T, y_list))
        a_i_alpha_T_r_matrix = self.a_i_alpha_T_r_matrix()
        a_ij = np.array([np.dot(i.T, i) for i in a_i_alpha_T_r_matrix.T])
        a_mix_matrix = y_ij * (a_ij ** 0.5)
        a_mix_list = np.array([np.sum(sum(i)) for i in a_mix_matrix])
        if self.debug:
            print('a_mix_list:', len(a_mix_list))
            print(a_mix_list)
        return a_mix_list
    
    def A_matrix(self) -> np.matrix:
        R = self.R
        T = self.T_list
        a = self.a_mix_list()
        a_RT = np.matrix(a / ((R ** 2) * (T ** 2))).T
        P = np.matrix(self.P_list)
        A_matrix = np.dot(a_RT, P)
        if self.debug:
            print('A_matrix:', A_matrix.shape)
            print(A_matrix)
        return A_matrix
    
    def b_i_T_r_matrix(self) -> np.matrix:
        R = self.R
        T_c = self.T_c_list
        P_c = self.P_c_list
        b_i_list = (0.08664 * R * T_c) / P_c
        b_i_T_r_matrix = np.matrix([b_i_list for i in range(self.T_r_list.__len__())]).T
        if self.debug:
            print('b_i_T_r_matrix:', b_i_T_r_matrix.shape)
            print(b_i_T_r_matrix)
        return b_i_T_r_matrix

    def b_mix_list(self) -> np.ndarray:
        y_list = np.matrix(self.y_list)
        b_i_T_r_matrix = self.b_i_T_r_matrix()
        b_mix_list = np.dot(y_list, b_i_T_r_matrix)
        if self.debug:
            print('b_mix_list:', len(b_mix_list))
            print(b_mix_list)
        return b_mix_list
    
    def B_matrix(self) -> np.matrix:
        R = self.R
        T = self.T_list
        b = self.b_mix_list()
        a_RT = np.matrix(b / (R * T)).T
        P = np.matrix(self.P_list)
        B_matrix = np.dot(a_RT, P)
        if self.debug:
            print('B_matrix:', B_matrix.shape)
            print(B_matrix)
        return B_matrix
    
    def Z_eq(self, Z, A, B) -> float:
        return (Z ** 3) - (Z ** 2) + ((A - B - (B ** 2)) * Z) - (A * B)
    
    def Z_matrix(self) -> np.matrix:
        A_matrix = self.A_matrix()
        B_matrix = self.B_matrix()
        if A_matrix.size == B_matrix.size:
            size = A_matrix.size
        else:
            ValueError('A_matrix和B_matrix的维度不一致')
        A_list = np.array(A_matrix).reshape(size)
        B_list = np.array(B_matrix).reshape(size)
        Z_matrix = []
        for i in range(size):
            A = A_list[i]
            B = B_list[i]
            Z = sc.optimize.fsolve(self.Z_eq, 1, args=(A, B))[0]
            Z_matrix.append(Z)
        Z_matrix = np.array(Z_matrix).reshape(A_matrix.shape)
        if self.debug:
            print('Z_matrix:', Z_matrix.shape)
            print(Z_matrix)
            print('size:', size)
        return np.matrix(Z_matrix)

    def Z_DF(self) -> pd.DataFrame:
        Z_DF = pd.DataFrame(self.Z_matrix(), index=self.T_r_list if self.is_index_r else self.T_list, columns=self.P_r_list if self.is_column_r else self.P_list)
        if self.debug:
            print('Z_DF:', Z_DF.shape)
            print(Z_DF)
        return Z_DF