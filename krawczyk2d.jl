using AutoDiff

function krawczyk2d(f, a::Array{Interval, 1}, bigprec::Integer)

arr_a = Array{Interval, 1}[]

function krawczyk2d_internal(f, a::Array{Interval, 1}, bigprec::Integer)

	set_bigfloat_precision(bigprec)
	tol = eps() #(10^10)*eps(BigFloat)	
	I = [Interval(1) Interval(0); Interval(0) Interval(1)]	
	Y(x) = make_intervals(inv(mid(jacobian(f, mid(x)))))
	K(x) = make_intervals(mid(x)) - Y(x)*make_intervals(f(mid(x))) + (I - Y(x)*jacobian(f, x))*(x - make_intervals(mid(x)))
	
	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	#if mid(a) == 0
	#	a = Interval(a.lo, a.hi + 0.0001*mag(a))
	#end
	
	if isect(a, K(a)) != false
		@show a
		@show diam(a)

		d = diam(a)
		dK = diam(K(a))
		
		if d == dK #d[1] < tol && d[2] < tol
			if dK[1] <= d[1] && dK[2] <= d[2] # Put <= instead of <
				println("Unique zero in $(K(a))")
				push!(arr_a, (K(a)))
			else
				println("Maybe a zero in $(K(a))")
			end
		
		else
			M1 = Y(a)
			if mid(det2(I - Y(K(a))*jacobian(f, K(a)))) <= mid(det2(I - Y(a)*jacobian(f, a)))
				krawczyk2d_internal(f, K(a), bigprec)
			else
				krawczyk2d_internal(f, make_intervals(mid(a)) - M1*make_intervals(f(mid(a))) + (I - M1*jacobian(f, a))*(a - make_intervals(mid(a))), bigprec)
			end

		end

	end
	
	return arr_a
end

return krawczyk2d_internal(f, a, bigprec)
end

