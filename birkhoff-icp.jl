include("birkhoffVN.jl")
using IntervalConstraintProgramming
f(x) = path_general(x, c, [1,2,3,1])
g(x) = f(x) - x

X = IntervalBox(Interval(-0.99, 0.99), Interval(0.0001, float(pi)/3))

orbit1231 = IntervalBox(Interval(-0.5001, -0.4998), Interval(float(pi)/6-0.004, float(pi)/6+0.007))

orbit1231exact = IntervalBox(Interval(-0.5), Interval(float(pi)/6))

S = @constraint abs(g(x)[1]) < 1e-3 && abs(g(x)[2]) < 1e-3
