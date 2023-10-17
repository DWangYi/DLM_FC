import pyforest

# ================================================================
### 变量说明：D——路段交通需求；N——路段总车道数；Nd——CAV专用道数量；p——路段CAV渗透率；K——车队规模; Dd——专用车道上的流量;\
### Dg——通用车道上的流量; pd——专用车道上的渗透率; pg——通用车道上的渗透率; vd——专用车道上的速度; vg——通用车道上的速度;\
### vd1——专用车道上的虚拟速度; vg1——通用车道上的虚拟速度
D = 1000  # 路段交通需求
N = 3  # 路段总车道数
Nd = 1  # CAV专用道数量
p = 0.1  # 路段CAV渗透率
K = 5  # 最大车队规模
Dd1 = p * D
Dg1 = (1 - p) * D
vd1 = vel(Dd1, cap(1, K))
vg1 = vel(Dg1, cap(0, K))


# ================================================================
## 单条车道通行能力计算函数，输入渗透率和最大车队规模，输出对应的通行能力
def cap(p, k):  # p为渗透率，k为最大排规模

    return 1800

# ================================================================
## 基于BPR的速度计算函数，输入交通需求和通行能力，输出对应的车道速度
def vel(Q, C):
    t_0 = 1  # 自由流速度行驶
    a = 1  # BPR函数中的alpha
    b = 1  # BPR函数中的beta
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
## 速度均衡函数，输出车辆转移比例
def get_trans_2():  # 针对策略2的二分法;
    pd = 1
    pg = 0
    rate_high = p
    rate_low = 0
    rate = 0
    diff = (p - rate) * D / cap(pd, K) - (1 - p + rate) * D / cap(pg, K)
    while abs(diff) > 0.001:
        rate = (rate_high - rate_low) * 0.5 + rate_low
        pd1 = 1
        pg1 = rate / (1 - p + rate)
        diff = (p - rate) / cap(pd1, K) - (1 - p + rate) / cap(pg1, K)
        if diff >= 0:
            rate_low = rate
        else:
            rate_high = rate
    return rate

def get_trans_3():  # 针对策略3的二分法;
    pd = 1
    pg = 0
    rate_high = 1 - p
    rate_low = 0
    rate = 0
    diff = (p + rate) / cap(pd, K) - (1 - p - rate) / cap(pg, K)
    while abs(diff) > 0.001:
        rate = (rate_high - rate_low) * 0.5 + rate_low
        pd1 = p / (p + rate)
        pg1 = 0
        diff = (p + rate) / cap(pd1, K) - (1 - p - rate) / cap(pg1, K)
        if diff >= 0:
            rate_high = rate
        else:
            rate_low = rate
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
D = 1000  # 路段交通需求
N = 3  # 路段总车道数
Nd = 1  # CAV专用道数量
p = 0.1  # 路段CAV渗透率
K = 5  # 最大车队规模

dis = flow_dis(2)
Dd = dis[0]
Dg = dis[1]
pd = dis[2]
pg = dis[3]
cap_d = cap(pd, K)   #专用车道通行能力
cap_g = (N-Nd)*cap(pg, K)   #通用车道通行能力
v_d = vel(Dd, cap_d)
v_g = vel(Dg, cap_g)
fuel_d = fuel_emission(0, v_d)
fuel_g = fuel_emission(0, v_g)
fuel_road = (fuel_d*Dd + fuel_g*Dg)/D