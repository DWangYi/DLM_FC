from pyforest import *
# ================================================================
## 变量说明：D——路段交通需求；N——路段总车道数；Nd——CAV专用道数量；p——路段CAV渗透率；K——车队规模; Dd——专用车道上的流量;\
## Dg——通用车道上的流量; pd——专用车道上的渗透率; pg——通用车道上的渗透率; vd——专用车道上的速度; vg——通用车道上的速度;\
## vd1——专用车道上的虚拟速度; vg1——通用车道上的虚拟速度

# ================================================================
## 单条车道通行能力计算函数，输入渗透率和最大车队规模，输出对应的通行能力
def cap(i, p, K):  # i为车道类型， p为渗透率，k为最大排规模，l为车身长度，S0为最小停车间距，vf为自由流速度，T为车辆安全车头时距
    if i == 1: #专用车道
        T = Hcc/K+(K-1)*Hcsc/K
    else: #混合车道
        A = 1-pow(p,K)
        B = 1-p
        C = 1-pow(p,K-1)
        D = p*(1-p)        
        T = B**2*Hhh + D*Hch + (p**2*C)*Hcsc/A + pow(p,K+1)*B*Hcc/A + D*Hhc
    E = (l+S0)/vf
    c = 1/(E+T)
    return c*3600

# ================================================================
## 基于BPR的速度计算函数，输入交通需求和通行能力，输出对应的车道速度
def vel(Q, C, vf):
    t_0 = 1/vf  # 自由流速度行驶
    a = 1.0122  # BPR函数中的alpha
    b = 4.1856  # BPR函数中的beta
    t = t_0 * (1 + a * ((Q / C) ** b))
    v = 1 / t
    return v

# ================================================================
## 油耗和排放计算函数，输入速度和指标编码，输出对应的油耗或者污染物排放
def fuel_emission(i, v):  # i为指标类, 0——fuel, 1——HC, 2——NOx, 3——CO, 4——CO2; v为交通流速度;
    fe_mtx = [[1560, 35.4, -0.388, 0.00776],
              [10.8, -0.00711, 0.000376, 0.0000363],
              [2.0, -0.0449, -0.000336, 0.0000349],
              [80.8, 1.16, 0.00503, 0.0000535],
              [4780, 111, -1.24, 0.0237]]
    fe_mtx = np.array(fe_mtx)
    a = fe_mtx[i, 0]
    b = fe_mtx[i, 1]
    c = fe_mtx[i, 2]
    d = fe_mtx[i, 3]
    res = a / v + b + c * v + d * (v ** 2)
    return res

# ================================================================
# 速度均衡函数，输出车辆转移比例
def get_trans_2():  # 针对策略2的二分法;
    pd = 1
    pg = 0
    rate_high = p
    rate_low = 0
    rate = 0
    diff = vel((p - rate) * D,Nd*cap(1, pd, K),vf) - vel((1 - p + rate) * D,(N-Nd)*cap(2, pg, K),vf)
    while abs(diff) > 0.1:
        rate = (rate_high - rate_low) * 0.5 + rate_low
        pd1 = 1
        pg1 = rate / (1 - p + rate)
        diff = vel((p - rate) * D,Nd*cap(1, pd1, K),vf) - vel((1 - p + rate) * D,(N-Nd)*cap(2, pg1, K),vf)
        if diff >= 0:
            rate_high = rate
        else:
            rate_low = rate
    return rate

def get_trans_3():  # 针对策略3的二分法;
    pd = 1
    pg = 0
    rate_high = 1 - p
    rate_low = 0
    rate = 0
    diff = vel((p + rate)*D,Nd*cap(1, pd, K),vf) - vel((1 - p - rate)*D,(N-Nd)*cap(2, pg, K),vf)
    while abs(diff) > 0.1:
        rate = (rate_high - rate_low) * 0.5 + rate_low
        pd1 = p / (p + rate)
        pg1 = 0
        diff = vel((p + rate)*D,Nd*cap(2, pd1, K),vf) - vel((1 - p - rate)*D,(N-Nd)*cap(2, pg1, K),vf)
        if diff >= 0:
            rate_low = rate
        else:
            rate_high = rate
    return rate

# ================================================================
## 流量分布计算函数，输入路段需求、路段车道设置方案、路段CAV渗透率，输出不同车道上的CAV渗透率
def flow_dis(i):
    if i == 1:  # 策略1
        Dd = p * D
        Dg = (1 - p) * D
        pd = 1
        pg = 0
    elif i == 2:  # 策略2
        if vd1 >= vg1:  # 当专用车道上的虚拟速度高于通用车道，此时CAV不会选择通用车道
            Dd = p * D
            Dg = (1 - p) * D
            pd = 1
            pg = 0
        else:
            alpha = get_trans_2()
            Dd = (p - alpha) * D
            Dg = (1 - p + alpha) * D
            pd = 1
            pg = alpha / (1 - p + alpha)
    elif i == 3:  # 策略3
        if vd1 < vg1:  # 当专用车道上的虚拟速度低于通用车道，此时HDV不会选择专用车道
            Dd = p * D
            Dg = (1 - p) * D
            pd = 1
            pg = 0
        else:
            beta = get_trans_3()
            Dd = (p + beta) * D
            Dg = (1 - p - beta) * D
            pd = p / (p + beta)
            pg = 0
    else:  # 基准策略，全部为混合车道
        Dd = Nd / N * D
        Dg = (1 - Nd) / N * D
        pd = p
        pg = p
    res = [Dd, Dg, pd, pg]
    return res

# ================================================================
##已知参数为路段交通需求、路段总车道数、CAVs专用道数量、路段CAV渗透率、最大车队规模
D = 5000 # 路段交通需求
N = 3  # 路段总车道数
Nd = 1  # CAV专用道数量
p = 0.8  # 路段CAV渗透率
K = 5  # 最大车队规模
Hhh = Hch = 1.52
Hhc = 1.35
Hcsc = 0.65
Hcc = 0.85 #车头时距
l = 5 #车身长度
S0 = 2 #最小停车间距
vf = 33.3 #自由流速度

Dd1 = p * D
Dg1 = (1 - p) * D
vd1 = vel(Dd1, Nd*cap(1, 1, K), vf) #虚拟速度
vg1 = vel(Dg1, (N-Nd)*cap(2, 0, K), vf)
dis = flow_dis(2) #策略i的流量分布
Dd = dis[0] #专用车道的需求
Dg = dis[1] #混合车道的需求
pd = dis[2] #专用车道的渗透率
pg = dis[3] #混合车道的渗透率
cap_d = Nd*cap(1, pd, K)   #专用车道通行能力
cap_g = (N-Nd)*cap(2, pg, K)   #通用车道通行能力
v_d = vel(Dd, cap_d, vf) #均衡速度
v_g = vel(Dg, cap_g, vf) #均衡速度

# fuel_d = fuel_emission(0, v_d)
# fuel_g = fuel_emission(0, v_g)
# fuel_ave = (fuel_d * Dd + fuel_g * Dg) / (Dd + Dg)

print('虚拟速度:'"{:.3f}".format(vd1),"{:.3f}".format(vg1))
print('均衡速度:'"{:.3f}".format(v_d),"{:.3f}".format(v_g))
