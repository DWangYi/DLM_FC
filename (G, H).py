# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 15:53:19 2023

@author: ASUS
"""
#(G, H)策略
import pyforest
from pyforest import *
import numpy
import pandas as pd
import decimal
import math

p = np.arange(0.01,1.00,0.01) 
D = np.arange(3000,7100,50) 

# 创建DataFrame的列名，第一列为'p', 后面的列为D的值
columns = ['p'] + [f'D{d}' for d in D]
df = pd.DataFrame(columns=columns)

for i in range(len(p)):
    row_data = [p[i]] + [0] * len(D)  # 初始化数据
    for j in range(len(D)):
        N = 2                         # 路段总车道数
        Ng = 1                        # GL数量
        K = 6                         # 最大车队规模
        Hhh = Hch = 2.0
        Hhc = 1.5
        Hcsc = 0.6
        Hcc = 1                       # 车头时距
        vf = 33.3                     # 自由流速度
        
        Dd1 = p[i] * D[j]             # CAV需求
        Dg1 = (1 - p[i]) * D[j]       # HDV需求
        
        
        #锥形容量延迟函数计算速度
        def vel(Q, C, vf):
            t_0 = 1/33.3                                     # 自由流速度行驶
            a = 4                                            # BPR函数中的alpha
            b = (2*a-1)/(2*a-2)                              # BPR函数中的beta
            A = math.sqrt(a*a*(1 - Q / C)*(1 - Q / C)+b*b)
            t = t_0 * (2 + A - a*(1 - Q / C) - b)
            v = 1 / t
            return v
        
        
        #通行能力函数
        def cap(m, p, K):  # i为车道类型， p为渗透率，k为最大排规模，vf为自由流速度
            if m == 1:                      # CDL
                T = Hcc/K+(K-1)*Hcsc/K
            else:                           # GL
                A = 1-p**K
                B = 1-p
                C = 1-pow(p,K-1)
                D = p*(1-p)        
                T = B**2*Hhh + D*Hch + (p**2*C)*Hcsc/A + pow(p,K+1)*B*Hcc/A + D*Hhc
                if p == 1:
                    T = Hcc/K+(K-1)*Hcsc/K
            c = 1/T
            return c*3600
        
        vg1 = vel(Dd1, Ng*cap(1, 1, K), vf)                  # GL虚拟速度
        vd1 = vel(Dg1, (N-Ng)*cap(2, 0, K), vf)              # HDL虚拟速度
        
    
    #(G, H)的流量分布
        if vg1 < vd1:                                        # 当通用车道上的虚拟速度低于专用车道，此时HDV不会选择通用车道
            Dd = p[i] * D[j]                                 # GL流量
            Dg = (1 - p[i]) * D[j]
            pd = 0
            pg = 1        
            cap_g = Ng*cap(1, 1, K)
        else:
            pd = 0
            pg = 1
            rate_high = 1 - p[i]
            rate_low = 0
            rate = 0
            e = 0.01
            diff = vel((p[i] + rate)*D[j], Ng*cap(1, pg, K),vf) - vel((1 - p[i] - rate)*D[j],(N-Ng)*cap(2, 0, K),vf)
            while np.all((abs(diff)) > e):
                rate = (rate_high - rate_low) * 0.5 + rate_low
                pg1 = p[i] / (p[i] + rate)
                pd1 = 0
                diff = vel((p[i] + rate)*D[j],Ng*cap(2, pg1, K),vf) - vel((1 - p[i] - rate)*D[j],(N-Ng)*cap(2, 0, K),vf)
                if np.all(diff >= 0):
                    rate_low = rate
                else:
                    rate_high = rate
            beta = rate
            Dd = (p[i] + beta) * D[j]  #g的流量
            Dg = (1 - p[i] - beta) * D[j]
            pg = p[i] / (p[i] + beta)
            pd = 0
            cap_g = Ng*cap(2, pg, K)   
        cap_d = (N-Ng)*cap(2, pd, K)
            
        
        v_g = vel(Dd, cap_g, vf) #g的均衡速度
        v_d = vel(Dg, cap_d, vf)
    
    
        #油耗函数
        def fuel_emission(i, v):  # i为指标类, 0——fuel, 1——HC, 2——NOx, 3——CO, 4——CO2; v为交通流速度;
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
             v1 = v *3.6                               # 换km/h
             res = a / v1 + b + c * v1 + d * (v1 ** 2) # 单位为g
             return res
         
        
        v = (v_d * Dg + v_g * Dd) / (Dd + Dg)
        fuel_ave = fuel_emission(0, v)
        
        row_data[j + 1] = fuel_ave                     # j+1 是因为第一列是 'p'
        df.loc[i] = row_data                           # 将数据添加到DataFrame
     
df.to_excel('lane3_f_GH.xlsx', index=False)