using KrawczykMethod2D
#=
f_p(x) = [arcsin(x[1] - x[2]) - float(pi)/2, sqrt1(x[1] + x[2]) - 3]

f_d1(x) = arcsin_d(x[1] - x[2])
f_d2(x) = sqrt_d(x[1] + x[2]) - 3
f_d(x) = [f_d1(x), f_d2(x)]
=#
bigprec = 64
set_bigfloat_precision(bigprec)


f(x) = [arcsin_d(x[1] - x[2]) - float(pi)/2, sqrt_d(x[1] + x[2]) - 3]
a = [Interval(3, 8), Interval(3, 7)]

krawczyk2d_loose(f, a)

# True root: [5, 4]

# A simpler function but it is uncoupled, which may complicate things

#f0(x) = [sqrt_d(x[1]) - 1, sqrt_d(2x[2]) - 2]
# True root: [1, 2]


