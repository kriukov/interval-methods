using IntervalArithmetic
using AutoDiff

# Bisection

left(x::Interval) = Interval(x.lo, mid(x))
right(x::Interval) = Interval(mid(x), x.hi)


function bisect(xx::Vector{Interval})
 
	if length(xx) != 2
	error("Only works for 2 at the moment")
	end

	x, y = xx

	intervals = Vector{Interval}[]

	push!(intervals, [left(x), left(y)])
	push!(intervals, [left(x), right(y)])
	push!(intervals, [right(x), left(y)])
	push!(intervals, [right(x), right(y)])

	intervals
end


function newton2d(f, a::MultiDimInterval, bigprec::Integer=64)

	set_bigfloat_precision(bigprec)
	
	# center() makes degenerate interval midpoints of the intervals in an array
	center(x) = make_intervals(mid(x))
	N(x) = center(x) - inv(jacobian(f, x))*f(center(x))

	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	#if mid(a) == 0
	#	a = Interval(a.lo, a.hi + 0.0001*mag(a))
	#end

	roots_array = MultiDimInterval[]

	push!(roots_array, a)

	k = 0
	while true

		roots_array_new = MultiDimInterval[]

		for i = 1:length(roots_array)
			if isectext(roots_array[i], N(roots_array[i])) != false
				if true #isectext(roots_array[i], N(roots_array[i])) != roots_array[i] - enabled only this part because doesn't work otherwise
					@show roots_array_new = push!(roots_array_new, isectext(roots_array[i], N(roots_array[i])))
				else
					@show roots_array_new = push!(roots_array_new, 
						isectext(bisect(roots_array[i])[1], N(bisect(roots_array[i])[1])), 
						isectext(bisect(roots_array[i])[2], N(bisect(roots_array[i])[2])), 
						isectext(bisect(roots_array[i])[3], N(bisect(roots_array[i])[3])), 
						isectext(bisect(roots_array[i])[4], N(bisect(roots_array[i])[4]))
					)
				end
			end
		end

		# Exit criterion
		if roots_array_new == roots_array 
			break
		end
		
		@show roots_array = roots_array_new
		k += 1
	end

	println("Function calls: ", k)
	return roots_array
end
