# Exercises from the Ramon Moore book

f(x) = [x[1]^2 + x[2]^2 - 1, x[1] - x[2]^2]
a = [Interval(0.5, 0.8), Interval(0.6, 0.9)]
krawczyk2d(f, a, 64)
trueroots = [(sqrt(5) - 1)/2, sqrt((sqrt(5) - 1)/2)]


f(x) = [20 - 20x[1] - x[2], x[1] - x[2]/20 - 1e-9*exp(x[2]/0.052)]
a = [Interval(0.5, 1.2), Interval(0.6, 1.2)]
krawczyk2d(f, a, 64)
trueroots = [0.9464142468335176, 1.0717150633296477]


f(x) = [0.5*(-(17.76x[1] - 103.79x[1]^2 + 229.62x[1]^3 - 226.31x[1]^4 + 83.72x[1]^5) + x[2]), 0.2*(-x[1] - 1.5x[2] + 1.2)]

a = [Interval(0.01, 0.1), Interval(0.7, 0.9)]
a = [Interval(0.7, 0.9), Interval(0.1, 0.3)]
a = [Interval(0.11, 0.3), Interval(0.5, 0.69)]

a = [Interval(0.01, 1), Interval(0.01, 1)]

trueroots = [0.06263595920972119, 0.758242693860186]
trueroots = [0.8844295888702942, 0.21038027408647064]
trueroots = [0.2853687241300338, 0.6097541839133109]
