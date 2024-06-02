import numpy as np
import pandas as pd
from plotly import graph_objects as go
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
import os
print(os.getcwd())

from .eq.PR_eq import PR
from .eq.RK_eq import RK
from .eq.SRK_eq import SRK

Gas_data =  pd.read_json(r'./EoS/Gas_data.json')
Gas_data.index = pd.Index(list(Gas_data['data']))
del Gas_data['data']
Gas_data = Gas_data.T
del Gas_data['k_PR_eq']
Gas_data['P_c'] *=  (10 ** 5)
Gas_data = Gas_data.T.apply(lambda row: row.fillna(row['common']), axis=1)

def k_matrix(k_list: np.ndarray) -> np.ndarray:
    k_matrix = np.zeros(shape=(len(k_list), len(k_list)))
    k_matrix[0, :] = k_list
    k_matrix[:, 0] = k_list
    return k_matrix

class EoS_data:
    def __init__(self, y_dict: dict, T_r_list: list[float]|float|None= None, P_r_list: list[float]|float|None= None, T_list: list[float]|float|None= None, P_list: list[float]|float|None= None, decimals: int= 4, debug: bool= False) -> None:
        '''
        y_ dict: 气体组成
        T_r_list: 对比态温度
        P_r_list: 对比态压力
        T_list: 绝对温度(K)
        P_list: 绝对压力(Pa)
        decimals: 保留小数位数
        debug: 是否调试
        '''
        if y_dict.keys() - set(Gas_data.columns) != set():
            raise KeyError('k_dict中包含未定义的气体', y_dict.keys() - set(Gas_data.columns))
        if round(sum(y_dict.values()), 6) != 1:
            raise ValueError('k_dict中k的值之和不为1', y_dict.values(), sum(y_dict.values()))
        self.y_dict = y_dict
        self.y_list = np.array(list(y_dict.values()), ndmin= 1)
        self.gas_list = np.array(list(y_dict.keys()), ndmin= 1)
        self.gas_data = Gas_data[self.gas_list]

        self.omega_list = np.array(np.around(self.gas_data.T['omega'], decimals), ndmin= 1)
        self.V_c_list = np.array(np.around(self.gas_data.T['V_c'], decimals), ndmin= 1)
        self.P_c_list = np.array(np.around(self.gas_data.T['P_c'], decimals), ndmin= 1)
        self.T_c_list = np.array(np.around(self.gas_data.T['T_c'], decimals), ndmin= 1)
        self.k_PR_list = np.array(np.around(self.gas_data.T['k_PR'], decimals), ndmin= 1)
        self.k_RK_list = np.array(np.around(self.gas_data.T['k_RK'], decimals), ndmin= 1)
        self.k_SRK_list = np.array(np.around(self.gas_data.T['k_SRK'], decimals), ndmin= 1)
        self.Omega_a_list = np.array(np.around(self.gas_data.T['Omega_a'], decimals), ndmin= 1)
        self.Omega_b_list = np.array(np.around(self.gas_data.T['Omega_b'], decimals), ndmin= 1)
        
        self.R = 8.314472 * (10 ** 6)

        self.T_sum = sum(self.y_list * self.T_c_list)
        self.P_sum = sum(self.y_list * self.P_c_list)


        self.T_r_list, self.T_list = self._validate_temperature(T_r_list, T_list, decimals)
        self.P_r_list, self.P_list = self._validate_pressure(P_r_list, P_list, decimals)

        self.k_PR_matrix = k_matrix(self.k_PR_list)
        self.k_RK_matrix = k_matrix(self.k_RK_list)
        self.k_SRK_matrix = k_matrix(self.k_SRK_list)

        self.debug = debug
        if self.debug:
            self._debug_info()

    def _validate_inputs(self, y_dict, T_r_list, P_r_list, T_list, P_list):
        if not isinstance(y_dict, dict):
            raise TypeError('y_dict应为字典类型')
        if not all(isinstance(k, (int, float)) for k in y_dict.values()):
            raise ValueError('y_dict的值应为数字类型')
        if y_dict.keys() - set(Gas_data.columns):
            raise KeyError('y_dict中包含未定义的气体', y_dict.keys() - set(Gas_data.columns))
        if round(sum(y_dict.values()), 6) != 1:
            raise ValueError('y_dict中组成比例之和不为1', y_dict.values(), sum(y_dict.values()))
        if T_r_list is not None and not isinstance(T_r_list, (list, float)):
            raise TypeError('T_r_list应为列表或浮点数类型')
        if P_r_list is not None and not isinstance(P_r_list, (list, float)):
            raise TypeError('P_r_list应为列表或浮点数类型')
        if T_list is not None and not isinstance(T_list, (list, float)):
            raise TypeError('T_list应为列表或浮点数类型')
        if P_list is not None and not isinstance(P_list, (list, float)):
            raise TypeError('P_list应为列表或浮点数类型')

    def _validate_temperature(self, T_r_list, T_list, decimals):
        if (T_list is None) and (T_r_list is not None):
            T_r_list = np.array(np.around(T_r_list, decimals), ndmin=1)
            T_list = T_r_list * self.T_sum
        elif (T_r_list is None) and (T_list is not None):
            T_list = np.array(np.around(T_list, decimals), ndmin=1)
            T_r_list = T_list / self.T_sum
        elif (T_list is None) and (T_r_list is None):
            raise ValueError('T_r_list和T_list不能同时为空')
        elif (T_list is not None) and (T_r_list is not None):
            T_r_list = np.array(np.around(T_r_list, decimals), ndmin=1)
            T_list = np.array(np.around(T_list, decimals), ndmin=1)
            if T_r_list.shape != T_list.shape:
                raise ValueError('T_r_list和T_list的形状不一致')
            if not np.allclose(T_list, T_r_list * self.T_sum):
                raise ValueError('T_r_list和T_list的值不匹配或其他原因')
        else:
            raise KeyError('T_r_list和T_list出错')
        return T_r_list, T_list

    def _validate_pressure(self, P_r_list, P_list, decimals):
        if (P_list is None) and (P_r_list is not None):
            P_r_list = np.array(np.around(P_r_list, decimals), ndmin=1)
            P_list = P_r_list * self.P_sum
        elif (P_r_list is None) and (P_list is not None):
            P_list = np.array(np.around(P_list, decimals), ndmin=1)
            P_r_list = P_list / self.P_sum
        elif (P_list is None) and (P_r_list is None):
            raise ValueError('P_r_list和P_list不能同时为空')
        elif (P_list is not None) and (P_r_list is not None):
            P_r_list = np.array(np.around(P_r_list, decimals), ndmin=1)
            P_list = np.array(np.around(P_list, decimals), ndmin=1)
            if P_r_list.shape != P_list.shape:
                raise ValueError('P_r_list和P_list的形状不一致')
            if not np.allclose(P_list, P_r_list * self.P_sum):
                raise ValueError('P_r_list和P_list的值不匹配或其他原因')
        else:
            raise KeyError('P_r_list和P_list出错')
        return P_r_list, P_list

    def _debug_info(self):
        print('y_dict:', len(self.y_dict))
        print(self.y_dict)
        print('gas_data:', len(self.gas_data))
        print(self.gas_data)
        print('T_r_list:', len(self.T_r_list))
        print(self.T_r_list)
        print('T_list:', len(self.T_list))
        print(self.T_list)
        print('P_r_list:', len(self.P_r_list))
        print(self.P_r_list)
        print('P_list:', len(self.P_list))
        print(self.P_list)

    def Z_DF(self, Eq: str= 'PR', is_P_r: bool = True, is_T_r: bool = True) -> dict:
        '''
        "DF": 返回Eq求解的Z_DF, 其中Z_DF的列为 P or P_r, 行为 T or T_r, 值为Z
        Eq: 求解的方程, 默认为'PR', 可选'RK', 'SRK'
        is_P_r: 是否以P_r为列, 默认为True
        is_T_r: 是否以T_r为行, 默认为True
        '''
        if Eq == 'PR':
            eq_class = PR(self.y_list, self.V_c_list,
                          self.omega_list, self.k_PR_matrix,
                          self.T_r_list, self.P_r_list, 
                          self.T_c_list, self.P_c_list, 
                          self.T_list, self.P_list, 
                          self.debug, is_P_r, is_T_r)
        elif Eq == 'RK':
            eq_class = RK(self.y_list, 
                          self.Omega_a_list, self.Omega_b_list, 
                          self.T_r_list, self.P_r_list, 
                          self.T_c_list, self.P_c_list, 
                          self.T_list, self.P_list, 
                          self.debug, is_P_r, is_T_r)
        elif Eq == 'SRK':
            eq_class = SRK(self.y_dict, 
                          self.omega_list, self.k_SRK_matrix,
                          self.T_r_list, self.P_r_list, 
                          self.T_c_list, self.P_c_list, 
                          self.T_list, self.P_list, 
                          self.debug, is_P_r, is_T_r)
        else:
            raise ValueError('计算公式错误, 可选"PR", "RK"或"SRK", 实为'+Eq)
        zDF = eq_class.Z_DF()
        return {'DF': zDF, 'Eq': Eq, 'is_P_r': is_P_r, 'is_T_r': is_T_r}
    
    def Z_figure(self, Eq: str, is_P_r: bool = True, is_T_r: bool = True, linemode: str= 'lines') -> go.Figure:
        '''返回Eq求解的Z_figure, 其中Z_figure的横坐标为 P or P_r, 纵坐标为Z, 线段为为 T or T_r'''
        zDF = self.Z_DF(Eq, is_P_r, is_T_r)['DF'].T
        linename = lambda name: 'T'+('_r' if is_T_r else '')+'='+str(round(name, 3))
        Z_Scatter_list = []
        for i in zDF.columns:
            Z_Scatter_list.append(go.Scatter(x=zDF[i].index, y=zDF[i].values, mode= linemode, name= linename(i)))
        Z_fig = go.Figure(Z_Scatter_list)
        Z_fig.update_layout(title=Eq, xaxis=dict(title='P_r' if is_P_r else 'P'), yaxis=dict(title='Z'))
        return Z_fig





