import numpy as np

from EoS.EoS import EoS_data

T_r_list = list(np.arange(1, 2.0 + 0.1, 0.1))
P_r_list = list(np.arange(0, 6.5 + 0.01, 0.1))
EoS_data_test = EoS_data(y_dict={'CO_2': 1}, P_r_list=P_r_list, T_r_list=T_r_list)
print('验证纯CO_2在P_r: 0.0-6.5, T_r: 1.0-2.0时的PR方程计算结果:')
zDF = EoS_data_test.Z_DF('PR')['DF']
zfig = EoS_data_test.Z_figure('PR', is_P_r= False, is_T_r= False)
print('在P_r=6.5时, P应为47931000.0: ', zfig.data[0].x[-1] == 47931000.0)
print('在T_r=2.0时, T应为608.24: ', eval(zfig.data[-1].name[2:]) == 608.24)
print('在P_r=1.5, T_r=1.2时, 计算得到的Z应为0.9057152117451431: ', zDF[1.2][1.5] == 0.9057152117451431)
print('显示图像, PR方程存在特殊解为正常情况, 示例图片见readme')
zfig.show()