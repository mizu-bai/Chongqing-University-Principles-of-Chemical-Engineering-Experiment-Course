## 调包
using CSV
using CurveFit
using Interpolations
using DataFrames
using PyCall
using PyPlot

## 调 Python 的包
@pyimport numpy as np
pygui(true)

## 读原始数据
"""
原始数据格式
|  q_v  |  t_1  |  t_2  |  T_1  |  T_2  |

变量名  变量含义              变量单位
q_v    换热管内流体流量       m ^ 3 * L ^ -1
t_1    冷空气进换热器处的温度  C
t_2    冷空气出换热器处的温度  C
T_1    冷空气进口处蒸汽的温度  C
T_1    冷空气出口处蒸汽的温度  C
"""

source_file_name = "data_of_rough_pipe.csv"
source_data = CSV.read(source_file_name)

## 温度计算
# 换热器壳程蒸汽的温度
source_data.T = [
    np.mean([source_data.T_1[i], source_data.T_2[i]])
    for i = 1:size(source_data)[1]
]
# 进口处温差 K
Δt_1 = source_data.T - source_data.t_1
# 出口处温差 K
Δt_2 = source_data.T - source_data.t_2
# 冷空气与蒸汽的对数平均温差计算 K
Δt_m = @. (Δt_1 - Δt_2) / log(Δt_1 / Δt_2)
# 定性温度 C
t_m = @. (source_data.t_1 + source_data.t_2) / 2

## 管路参数计算
# 管内径 m
d = (21 - 2.5 * 2) * 1e-3
# 管程长 m
l = 1000e-3
# 流速 m / s
u = @. source_data.q_v / (π / 4 * d ^ 2 * 3600)
# 表面积 m * s ^ (-1)
A = π * d * l + 0.03987

## 气体物性计算
# 空气数据
air = CSV.read("air.csv")
# 密度插值
interp_ρ = Interpolations.LinearInterpolation(air.t, air.rho)
# 空气密度 kg * m ^ (-3)
ρ = interp_ρ(t_m)
# 热导率插值
interp_λ = Interpolations.LinearInterpolation(air.t, air.lambda)
# 空气热导率 W * m ^ (-1) * K ^ (-1)
λ = interp_λ(t_m)
# 比热容插值
interp_c_p = Interpolations.LinearInterpolation(air.t, air.c_p)
# 空气比热容 J * (kg * K) ^ (-1)
c_p = @. interp_c_p(t_m)
# 空气粘度插值
interp_μ = Interpolations.LinearInterpolation(air.t, air.mu)
# 空气粘度 Pa * s
μ = @. interp_μ(t_m)

## 传热参数埏
# 传热量
Q = @. c_p * ρ * source_data.q_v * (source_data.t_2 - source_data.t_1) / 3600
# 总传热系数
K = @. Q / (A * Δt_m)

## 准数计算
Nu = @. K * d / λ
Re = @. d * u * ρ / μ
Pr = @. c_p * μ / λ

## 线性拟合
x = @. log(Re)
y = @. log(Nu / Pr ^ 0.4)
p = linear_fit(x, y)
yy = @. p[1] + p[2] * x
println("y = $(p[1]) + $(p[2]) * x")

## 决定系数
y_avg = np.mean(y)
SSreg = sum(@. (yy - y_avg) ^ 2)
SStot = sum(@. (y - y_avg) ^ 2)
R2 = SSreg / SStot
println("R2 = $R2")

## 绘图
PyPlot.plot(x, y, "o", x, yy, "-")
PyPlot.xlabel("\$ \\ln{Re} \$")
PyPlot.ylabel("\$ \\ln \\frac{Nu}{Pr ^ {0.4}} \$")
PyPlot.title("\$ \\ln \\frac{Nu}{Pr ^ {0.4}} - \\ln{Re} \$ figure")

## 输出结果到文件
result = DataFrame(
    :q_v => source_data.q_v,
    :Q => Q,
    :K => K,
    :Nu => Nu,
    :Re => Re,
    :Pr => Pr,
)
target_file_name = "result_of_rough_pipe.csv"
CSV.write(target_file_name, result)
