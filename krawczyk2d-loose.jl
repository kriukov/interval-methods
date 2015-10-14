using IntervalArithmetic
using AutoDiff

function all_inside(x::MultiDimInterval, y::MultiDimInterval)
	k = 0
	for i = 1:length(x)
		if !inside(x[i], y[i])
			k += 1
		end
	end
	if k > 0
		return false
	else
		return true
	end
end

function bisect(xx::MultiDimInterval)
 	if length(xx) != 2
	    error("Only works for 2 at the moment")
	end

	x, y = xx

	intervals = MultiDimInterval[]

	push!(intervals, [left(x), left(y)])
	push!(intervals, [left(x), right(y)])
	push!(intervals, [right(x), left(y)])
	push!(intervals, [right(x), right(y)])

	intervals
end

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

