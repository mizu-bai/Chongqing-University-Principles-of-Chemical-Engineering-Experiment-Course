## 调包
using CurveFit
using CSV
using Interpolations
using PyCall
using PyPlot
using DataFrames

## 调 python 包
@pyimport numpy as np
pygui(true)

## 原始数据
"""
原始数据输入
(进料组成, 塔顶产品组成, 塔釜组成) 的质量分数 (ω_f, ω_d, ω_w)
进料温度 t_f
回流流量 L
产品流量 D
"""

# (进料组成, 塔顶产品组成, 塔釜组成) 的质量分数
ω_f = 16.2419e-2
ω_d = 89.8466e-2
ω_w = 6.0228e-2

# 进料温度 / C
t_f = 33.9

# 回流流量 / (L * h ^ -1)
L = 3.6
# 产品流量 / (L * h ^ -1)
D = 1.5

# 质量分数 -> 摩尔分数 换算函数
f(ω_i) = (ω_i / 46.07) / (ω_i / 46.07 + (1 - ω_i) / 18.02)

# (进料组成, 塔顶产品组成, 塔釜组成) 摩尔分数
x_f = f(ω_f)
x_d = f(ω_d)
x_w = f(ω_w)

## 对角线
x = y = np.linspace(0, 1)

## 平衡线
# EtOH-H2O 相平衡数据
EtOH = CSV.read("EtOH-H2O.csv", DataFrame)
# EtOH-H2O 相平衡组成插值
interp_y = Interpolations.LinearInterpolation(EtOH.x, EtOH.y)
# 相平衡线
yy = interp_y(x)

## 精馏段操作线
R = L / D
f_r(x) = (R * x + x_d) / (R + 1)

## q 线
# 泡点温度插值
interp_t = Interpolations.LinearInterpolation(EtOH.x, EtOH.t)
# 泡点温度
t_s = interp_t(x_f)

# 乙醇摩尔热容
c_pm_EtOH(t) = (1.56 * 10e-2 * t + 2.012) * 46.07
# 水摩尔热容
c_pm_H2O(t) = (2.143 * 10e-4 * t + 4.198) * 18.02
# 溶液热容 / (kJ * kmol ^ -1 * C ^ -1)
t_m = (t_f + t_s) / 2
c_pm = x_f * c_pm_EtOH(t_m) + (1 - x_f) * c_pm_H2O(t_m)

# 乙醇摩尔气化潜热
r_m_EtOH(t) = 113 * (243 - t) ^ 0.4218 * 46.07
# 水的摩尔气化潜热
r_m_H2O(t) = 445.6 * (374 - t) ^ 0.3003 * 18.02
# 溶液气化潜热 / (kJ * kmol ^ -1)
r_m = x_f * r_m_EtOH(t_f) + (1 - x_f) * r_m_EtOH(t_f)

# q 值
q = 1 + c_pm * (t_s - t_f) / r_m

# q 线方程
f_q(x) = (q * x - x_f) / (q - 1)

## 交点
# 精馏段操作线和对角线交点
point_a = [x_d, x_d]
# 提馏段操作线和对角线的交点
point_b = [x_w, x_w]
# 精馏段操作线和 y 轴的交点
point_c = [0, f_r(0)]
# q 线和精馏段操作线的交点
point_d = [((R + 1) * x_f + (q - 1) * x_d) / (R + q), (R * x_f + q * x_d) / (R + q)]
# q 线和对角线的交点
point_f = [x_f, x_f]

## 提馏段操作线方程
p_d = CurveFit.linear_fit([point_b[1], point_d[1]], [point_b[2], point_d[2]])
f_d(x) = p_d[1] + p_d[2] * x

## 绘图
# 对角线
PyPlot.plot(x, y, "-")
# 相平衡线
PyPlot.plot(x, yy, "-")
# 精馏段操作线 a-c
PyPlot.plot([point_a[1], point_c[1]], [point_a[2], point_c[2]], "-")
# q 线 d-f
PyPlot.plot([point_d[1], point_f[1]], [point_d[2], point_f[2]], "-")
# 提馏段操作线 b-d
PyPlot.plot([point_b[1], point_d[1]], [point_b[2], point_d[2]], "-")
# x_w 垂直线
PyPlot.plot([x_w, x_w], [0, x_w], "b-")
# x_d 垂直线
PyPlot.plot([x_d, x_d], [0, x_d], "b-")
# x_f 垂直线
PyPlot.plot([x_f, x_f], [0, x_f], "b-")

# 绘图设置
PyPlot.xlabel("\$ x \$")
PyPlot.ylabel("\$ y \$")
PyPlot.xlim([0, 1])
PyPlot.ylim([0, 1])

## 塔板数计算
# 迭代点
point_1 = [x_d, x_d]
point_2 = zeros(2)
# 塔板数
num = 0
# 迭代记数
n = 1
# 迭代计算 & 绘图
println("开始迭代")
while true
    global point_1, point_2, num, n
    if point_1[1] < x_w
        break
    end
    println("迭代第 $(n) 次")
    n += 1
    for x_i in np.linspace(0, 1, 10000000)
        if abs(interp_y(x_i) - point_1[2]) < 0.00001 * point_1[2]
            point_temp = [x_i, point_1[2]]
            point_2 = zeros(2)
            if point_d[1] < x_i
                # 在精馏段操作线上
                point_2 = [x_i, f_r(x_i)]
            elseif point_b[1] < x_i < point_d[1]
                # 在提馏段操作线上
                point_2 = [x_i, f_d(x_i)]
            elseif x_i < point_b[1]
                # 在对角线上
                point_2 = [x_i, x_i]
            end
            num += 1
            PyPlot.plot([point_1[1], point_temp[1]], [point_1[2], point_temp[2]], "b-")
            PyPlot.plot([point_2[1], point_temp[1]], [point_2[2], point_temp[2]], "b-")
            # 迭代后新点
            point_1 = point_2
            point_2 = zeros(2)
            break
        end
    end
end

# 理论塔板数结果
println("理论塔板数 = $(num)")
