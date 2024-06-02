import numpy as np
import pandas as pd
import scipy as sc
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

class RK:
    def __init__(self, y_list, Omega_a_list, Omega_b_list, T_r_list, P_r_list, T_c_list, P_c_list, T_list, P_list, debug, is_column_r, is_index_r) -> None:
        self.R = 8.314472 * (10 ** 6)
        self.T_r_list = T_r_list
        self.P_r_list = P_r_list
        self.T_c_list = T_c_list
        self.P_c_list = P_c_list
        self.T_list = T_list
        self.P_list = P_list
        self.Omega_a_list = Omega_a_list
        self.Omega_b_list = Omega_b_list
        self.y_list = y_list
        self.debug = debug
        self.is_column_r = is_column_r
        self.is_index_r = is_index_r


    def alpha_T_r_list(self) -> np.ndarray:
        T_r = np.array(self.T_r_list)
        alpha_T_r_list = T_r ** 0
        if self.debug:
            print('alpha_T_r_list:', len(alpha_T_r_list))
            print(alpha_T_r_list)
        return alpha_T_r_list

    def a_i_alpha_T_r_matrix(self) -> np.matrix:
        alpha_T_r_list = self.alpha_T_r_list()
        R = self.R
        T_c = self.T_c_list
        P_c = self.P_c_list
        Omega_a = self.Omega_a_list
        a_i_list = (Omega_a *( R ** 2) * (T_c ** 2.5)) / P_c
        a_i_alpha_T_r_matrix = np.dot(np.matrix([i for i in a_i_list]).T, [alpha_T_r_list])
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
        a_RT = np.matrix(a / ((R ** 2) * (T ** 2.5))).T
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
        Omega_b = self.Omega_b_list
        b_i_list = (Omega_b * R * T_c) / P_c
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
        return np.matrix(Z_matrix)

    def Z_DF(self) -> pd.DataFrame:
        Z_DF = pd.DataFrame(self.Z_matrix(), index=self.T_r_list if self.is_index_r else self.T_list, columns=self.P_r_list if self.is_column_r else self.P_list)
        if self.debug:
            print('Z_DF:', Z_DF.shape)
            print(Z_DF)
        return Z_DF
