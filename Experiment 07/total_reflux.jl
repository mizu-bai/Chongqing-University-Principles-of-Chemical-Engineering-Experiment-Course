## 调包
using CSV
using Interpolations
using PyCall
using PyPlot

## 调 python 包
@pyimport numpy as np
pygui(true)

## 原始数据
"""
原始数据输入
(塔顶产品组成, 塔釜组成) 的质量分数 (x_d, x_w)
"""

# (塔顶产品组成, 塔釜组成) 的质量分数
ω_d = 93.5107e-2
ω_w = 12.4097e-2

# 质量分数 -> 摩尔分数 换算函数
f(ω_i) = (ω_i / 46.07) / (ω_i / 46.07 + (1 - ω_i) / 18.02)

# (塔顶产品组成, 塔釜组成) 摩尔分数
x_d = f(ω_d)
x_w = f(ω_w)

## 对角线
x = y = np.linspace(0, 1, 1000)

## 平衡线
# EtOH-H2O 相平衡数据
EtOH = CSV.read("EtOH-H2O.csv")
# EtOH-H2O 相平衡组成插值
interp_y = Interpolations.LinearInterpolation(EtOH.x, EtOH.y)
## 相平衡线
yy = interp_y(x)

## 绘图
# 对角线
PyPlot.plot(x, y, "-")
# 相平衡线
PyPlot.plot(x, yy, "-")
# x_w 垂直线
PyPlot.plot([x_w, x_w], [0, x_w], "b-")
# x_d 垂直线
PyPlot.plot([x_d, x_d], [0, x_d], "b-")

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
            pointTemp = [x_i, point_1[2]]
            point_2 = [x_i, x_i]
            num += 1
            PyPlot.plot([point_1[1], pointTemp[1]], [point_1[2], pointTemp[2]], "b-")
            PyPlot.plot([point_2[1], pointTemp[1]], [point_2[2], pointTemp[2]], "b-")
            # 迭代后新点
            point_1 = point_2
            point_2 = zeros(2)
            println(point_1[1])
            break
        end
    end
end

# 理论塔板数结果
println("理论塔板数 = $(num)")
