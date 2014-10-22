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

function newton2d(f, a::Array{Interval, 1}, bigprec::Integer=64)

	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	#if mid(a) == 0
	#	a = Interval(a.lo, a.hi + 0.0001*mag(a))
	#end

	roots_array = Array{Interval, 1}[]

	#push!(roots_array, a)

	tol = 1e-10
println("point 1")
	function newton2d_internal(f, a::Array{Interval, 1}, bigprec::Integer)


		set_bigfloat_precision(bigprec)

		# center() makes degenerate interval midpoints of the intervals in an array
		center(x) = make_intervals(mid(x))
		N(x) = center(x) - inv(jacobian(f, x))*f(center(x))


		k = 0
		while diam(a)[1] >= tol || diam(a)[2] >= tol
			@show diam(a)
			#roots_array_new = Array{Interval, 1}[]

			for i = 1:length(roots_array)
				if isectext(roots_array[i], N(roots_array[i])) != false
					#if isectext(roots_array[i], N(roots_array[i])) != roots_array[i]
						@show newton2d_internal(f, isectext(roots_array[i], N(roots_array[i])), bigprec)
					#else
						#@show newton2d_internal(f, isectext(bisect(roots_array[i])[1], N(bisect(roots_array[i])[1])), bigprec)
						#@show newton2d_internal(f, isectext(bisect(roots_array[i])[2], N(bisect(roots_array[i])[2])), bigprec)
						#@show newton2d_internal(f, isectext(bisect(roots_array[i])[3], N(bisect(roots_array[i])[3])), bigprec)
						#@show newton2d_internal(f, isectext(bisect(roots_array[i])[4], N(bisect(roots_array[i])[4])), bigprec)
					#end
				end
			end

			# Exit criterion
			#=
			if roots_array_new == roots_array
				for i = 1:length(roots_array)
					push!(roots_array_new, isectext(roots_array[i], N(roots_array[i])))
				end 
				break
			end
			=#
			#@show roots_array = roots_array_new
			@show k += 1
		end
		
		return roots_array

	end

#println("Function calls: ", k)
return newton2d_internal(f, a, bigprec)

end

