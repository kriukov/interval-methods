using NewtonMethod
using KrawczykMethod

f(x) = 5040 - 3828x - 2356x^2 + 1009x^3 + 200x^4 - 62x^5 - 4x^6 + x^7
a = Interval(-20, 20)
@time newton(f, a, 128)
@time krawczyk(f, a, 128)

f(x) = x^2 - 4
@time newton(f, a, 128)
@time krawczyk(f, a, 128)

f(x) = exp(x^2) - 4
a = Interval(-4, 4)
@time newton(f, a, 128)
@time krawczyk(f, a, 128)

f(x) = log(x^2 - 1) + x + 1
a = Interval(-10, -1.1)
@time newton(f, a, 128)
@time krawczyk(f, a, 128)

a = Interval(1.1, 10)
@time newton(f, a, 128)
@time krawczyk(f, a, 128)
