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

	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	#if mid(a) == 0
	#	a = Interval(a.lo, a.hi + 0.0001*mag(a))
	#end

	center(x) = make_intervals(mid(x))
	N(x) = center(x) - inv(jacobian(f, x))*f(center(x))

	roots_array = MultiDimInterval[]

	push!(roots_array, a)

	@show tol = 1e-10
	
	k = 0
println("point 1")
	function newton2d_internal(f, a::MultiDimInterval, bigprec::Integer)


		set_bigfloat_precision(bigprec)
		
		@show a
		@show N(a)
		@show Na = isectext(a, N(a))
			
			if Na != false
			
				@show d = diam(a)
				@show dN = diam(Na)
				
				if dN[1] < tol && dN[2] < tol 
					if all_inside(Na, a)
						println("Unique zero in $Na")
						push!(roots_array, Na)
					else
						println("Maybe a zero in $Na")
					end
					
				else
					
					if isectext(roots_array, N(roots_array)) != roots_array
					 	for i = 1:length(roots_array)
							newton2d_internal(f, roots_array[i], bigprec)
						end
					else 
					
						for i = 1:length(roots_array)
							@show k += 1
					
							newton2d_internal(f, bisect(roots_array[i])[1], bigprec)
							newton2d_internal(f, bisect(roots_array[i])[2], bigprec)
							newton2d_internal(f, bisect(roots_array[i])[3], bigprec)
							newton2d_internal(f, bisect(roots_array[i])[4], bigprec)
						end
						
					end
				end
			end
		
		return roots_array

	end

#println("Function calls: ", k)
return newton2d_internal(f, a, bigprec)

end

