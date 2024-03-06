
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 13:41:48 2023

@author: ASUS
"""
#(C, H)策略
import pyforest
from pyforest import *
import numpy
import pandas as pd
import math

p = np.arange(0.01,1.00,0.01) 
D = np.arange(3000,3100,100) 

# 创建DataFrame的列名，第一列为'p', 后面的列为D的值
columns = ['p'] + [f'D{d}' for d in D]
df = pd.DataFrame(columns=columns)

for i in range(len(p)):
    row_data = [p[i]] + [0] * len(D)  # 初始化数据
    for j in range(len(D)):
        N = 3                         # 路段总车道数
        Nd = 1                        # CDL数量
        K = 6                         # 最大车队规模
        Hhh = Hch = 2.0
        Hhc = 1.5
        Hcsc = 0.6
        Hcc = 1                       # 车头时距
        vf = 33.3                     # 自由流速度
        Dd = p[i] * D[j]              # CAV的需求
        Dg = (1 - p[i]) * D[j]        # HDV的需求
        
        
        #锥形容量延迟函数计算速度
        def vel(Q, C, vf):
            t_0 = 1/33.3                                     # 自由流速度行驶
            a = 4                                            # BPR函数中的alpha
            b = (2*a-1)/(2*a-2)                              # BPR函数中的beta
            A = math.sqrt(a*a*(1 - Q / C)*(1 - Q / C)+b*b)
            t = t_0 * (2 + A - a*(1 - Q / C) - b)
            v = 1 / t
            return v                                         # m/s
        
        
        #通行能力函数
        def cap(m, p, K):               # i为车道类型， p为渗透率，k为最大排规模，vf为自由流速度
            if m == 1: 
                T = Hcc/K+(K-1)*Hcsc/K  # CDL
            else:  
                T = Hhh                 # HDL
            c = 1/T
            return c*3600
        
        
        pd=1
        cap_d = Nd*cap(1, pd, K)                 # CDL通行能力
        cap_g = (N-Nd)*cap(2, 0, K)              # HDL通行能力
        v_d = vel(Dd, Nd*cap(1, 1, K), vf)       # CDL速度
        v_g = vel(Dg, (N-Nd)*cap(2, 0, K), vf)   # HDL车道
        
        
        #油耗函数
        def fuel_emission(i, v):  # i为指标类, 0——fuel, 1——HC, 2——NOx, 3——CO, 4——CO2; v为交通流速度
            fe_mtx = [[1560, 35.4, -0.0388, 0.000776],
                      [10.8, -0.00711, 0.000376, 0.0000363],
                      [2.0, -0.0449, -0.000336, 0.0000349],
                      [80.8, 1.16, 0.00503, 0.0000535],
                      [4780, 111, -1.24, 0.0237]]
            fe_mtx = np.array(fe_mtx)
            a = fe_mtx[i, 0]
            b = fe_mtx[i, 1]
            c = fe_mtx[i, 2]
            d = fe_mtx[i, 3]
            v1 = v *3.6                                 # km/h
            res = a / v1 + b + c * v1 + d * (v1 ** 2)   # 单位为g
            return res

        v = (v_d * Dd + v_g * Dg) / (Dd + Dg)
        fuel_ave = fuel_emission(0, v)
        
        row_data[j + 1] = fuel_ave  # j+1 是因为第一列是 'p'
        df.loc[i] = row_data        # 将数据添加到DataFrame

df.to_excel('lane3_f_CH.xlsx', index=False)