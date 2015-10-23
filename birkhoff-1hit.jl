using IntervalArithmetic
using AutoDiff
using PurityIntervals

# The vector x = [w, th] in Birkhoff mapping

function alpha(n, m)
	if n > 0 && n <= 3 && m > 0 && m <= 3 && typeof(n) == Int && typeof(m) == Int && n != m
		k = pi/3
		if n == 1 && m == 2
			return 0
		elseif n == 1 && m == 3
	 		return k
		elseif n == 2 && m == 3
	 		return 2k
		elseif n == 2 && m == 1
	 		return 3k
		elseif n == 3 && m == 1
	 		return 4k
		elseif n == 3 && m == 2
	 		return 5k
	 	end
	else 
	    error("Invalid parameters: should be 1, 2, 3 and not equal")
	end
end

#=
function T(x, n, m, r)
	wnext = x[1] - r*(x[1]*cos(x[2] - alpha(n, m)) + sqrt_d(1 - x[1]^2)*sin(x[2] - alpha(n, m)))
	[wnext, mod(x[2] + float(pi) + arcsin_d(x[1]) + arcsin_d(wnext), 2pi)]
end
=#

function T(x, n, m, r)
	wnext = x[1] - r*(x[1]*cos(x[2] - alpha(n, m)) + sqrt(1 - x[1]^2)*sin(x[2] - alpha(n, m)))
	[wnext, mod(x[2] + float(pi) + asin(x[1]) + asin(wnext), 2pi)]
end


f12(x) = T(x, 1, 2, 2.1)

# Split the phase space w = [-1, 1], th = [0, pi/3] into 30x30 rectangles 
N = 30
rw = Interval(-0.95, 0.95)
rth = Interval(0.05, float(pi)/3-0.05)
deltaw = diam(rw/N)
deltath = diam(rth/N)

points = Array{BigFloat, 1}[]

for i = 1:30
	for j = 1:30
		rect = [rw.lo + i*deltaw - Interval(0, deltaw), rth.lo + j*deltath - Interval(0, deltath)]
		p = purity(f12, rect)
		#println([rect, p])
		if p == 1 || p == 0
			push!(points, mid(rect))
		end
	end
end

println(points)
