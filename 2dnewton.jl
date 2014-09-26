using AutoDiff

center(x) = Interval(mid(x))
N(x) = x - inv(J(x))*f(x)



