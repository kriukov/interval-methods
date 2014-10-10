using AutoDiff


function newton2d(f, a::Array{Interval, 1}, bigprec::Integer=64)

	set_bigfloat_precision(bigprec)

	#center(x) = Interval(mid(x))
	N(x) = x - make_intervals(inv(mid(jacobian(f, mid(x)))))*f(x)

	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	#if mid(a) == 0
	#	a = Interval(a.lo, a.hi + 0.0001*mag(a))
	#end

	roots_array = Array{Interval, 1}[]

	push!(roots_array, a)

	k = 0
	while true

		roots_array_new = Array{Interval, 1}[]

		for i = 1:length(roots_array)
			if isectext(roots_array[i], N(roots_array[i])) != false
				@show roots_array_new = push!(roots_array_new, isectext(roots_array[i], N(roots_array[i])))
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
