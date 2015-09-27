## Newton's method (interval) code

module NewtonMethod
export newton, differentiate, Interval, rad, diam, mid, mig, mag, belong, hd, hull, isect, isectext

using IntervalArithmetic
using AutoDiff

println("Syntax: newton(function, Interval(lo, hi), precision [default is 64])")

function newton(f::Function, a::Interval, bigprec::Integer=64)

	set_bigfloat_precision(bigprec)

	center(x) = Interval(mid(x))
	N(x) = center(x) - f(center(x))//differentiate(f, x)

	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	if mid(a) == 0
		a = Interval(a.lo, a.hi + 0.0001*mag(a))
	end

	roots_array = Interval[]

	push!(roots_array, a)

	k = 0
	while true

		roots_array_new = Interval[]

		for i = 1:length(roots_array)
			if isectext(roots_array[i], N(roots_array[i])) != false
				roots_array_new = vcat(roots_array_new, isectext(roots_array[i], N(roots_array[i])))
			end
		end

		# Exit criterion
		if roots_array_new == roots_array 
			break
		end
	
		roots_array = roots_array_new
		k += 1
	end

	println("Function calls: ", k)
	return roots_array
end

# end of module
end
