using IntervalArithmetic
using AutoDiff

intdet(M) = M[1]*M[4] - M[2]*M[3]

bigprec=64
roots_array = MultiDimInterval[]

set_bigfloat_precision(bigprec)
tol = 1e-6

I = [Interval(1) Interval(0); Interval(0) Interval(1)]
mdelta = [Interval(1e-5) Interval(1e-5); Interval(1e-5) Interval(1e-5)]

# If the jacobian is non-invertible, the SingularException error is returned for Y. We need to choose a slightly different Y then.
function Y(f, x)
    if intdet(jacobian(f, mid(x))) == Interval(0)
	    #return make_intervals(inv(jacobian(f, mid(x) + 0.0001*norm(diam(x)))))
	    return make_intervals(inv(jacobian(f, mid(x)) + mdelta))
    else
	    return make_intervals(inv(jacobian(f, mid(x))))
    end
end

M(f, x) = I - Y(f, x)*jacobian(f, x)
K(f, x) = make_intervals(mid(x)) - Y(f, x)*make_intervals(f(mid(x))) + M(f, x)*(x - make_intervals(mid(x)))


bigprec = 64
set_bigfloat_precision(bigprec)


f(x) = [arcsin_d(x[1] - x[2]) - float(pi)/2, sqrt_d(x[1] + x[2]) - 3]
a = [Interval(3, 8), Interval(3, 7)]

a2 = bisect(a)[2]
Y(f, a2)
jacobian(f, mid(a2))
det(jacobian(f, mid(a2)))

# det() blows up with 0-interval matrices
B = [Interval(0, 0) Interval(5,6); Interval(7, 8) Interval(9,11)]
det(B)

# det() produces semi-inf interval if there is a 0
B = [Interval(0, 3) Interval(5,6); Interval(7, 8) Interval(9,11)]
det(B)

B11 = Interval(0, 0); B12 = Interval(5,6); B21 = Interval(7, 8); B22 = Interval(9,11)

using KrawczykMethod2D; f(x) = [asin(x[1] - x[2]) - float(pi)/2, sqrt(x[1] + x[2]) - 3]; a = [Interval(3, 8), Interval(3, 7)]; krawczyk2d_purity(f, a)


#krawczyk2d detects "unique zero" in a = [IntervalArithmetic.Interval(4.00146484375000000000,4.00170898437500000000),IntervalArithmetic.Interval(3.00170898437500000000,3.00195312500000000000)]. The interval is clean, but it is not a zero!
