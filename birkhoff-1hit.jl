	push!(LOAD_PATH, "modules")

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
	wnext = ω - r*(ω*cos(θ - alpha(n, m)) + sqrt_d(1 - ω^2)*sin(θ - alpha(n, m)))
	[wnext, mod(θ + float(pi) + arcsin_d(x[1]) + arcsin_d(wnext), 2pi)]
end
=#

function T(x, n, m, r)
	ω, θ = x

	ω_next = ω - r*(ω*cos(θ - alpha(n, m)) + √(1 - ω^2)*sin(θ - alpha(n, m)))
	#θ_next = mod(θ + big(pi) + asin(ω) + asin(ω_next), 2π)
	θ_next = θ + big(pi) + asin(ω) + asin(ω_next)


	[ω_next, θ_next]
end


using PyPlot

distance_between_centers = 6.


f12(x) = T(x, 1, 2, distance_between_centers)
f13(x) = T(x, 1, 3, distance_between_centers)

function twice_123(x0)
	x1 = T(x0, 1, 2, distance_between_centers)
	x2 = T(x1, 2, 3, distance_between_centers)

	x2

end


function plot_band(f)

	# Split the phase space w = [-1, 1], th = [0, pi/3] into 30x30 rectangles
	N = 30
	rw = Interval(-0.95, 0.95)
	rth = Interval(0.05, float(pi)/3-0.05)
	deltaw = diam(rw/N)
	deltath = diam(rth/N)

	points = Array{BigFloat, 1}[]
	purities = Int[]

	for i = 1:30
		for j = 1:30
			rect = [rw.lo + i*deltaw - Interval(0, deltaw), rth.lo + j*deltath - Interval(0, deltath)]
			p = purity(f, rect)
			#println([rect, p])
			push!(points, mid(rect))
			push!(purities, p)

		end
	end

	x = [Float64(xx[1]) for xx in points]
	y = [Float64(xx[2]) for xx in points]

	p = purities

	for i in (0, 1, -1)
		plot(y[p.==i], x[p.==i], "o", label="$i", alpha = 0.5)
	end
	legend()

end

#plot_band(f12)
#plot_band(f13)

plot_band(twice_123)
