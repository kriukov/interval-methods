using AutoDiff

function krawczyk2d(f, a::Array{Interval, 1}, bigprec::Integer)

arr_a = Array{Interval, 1}[]

set_bigfloat_precision(bigprec)
#tol = (1e-10)*eps(BigFloat)
tol = 1e-10

I = [Interval(1) Interval(0); Interval(0) Interval(1)]
Y(x) = make_intervals(inv(mid(jacobian(f, mid(x)))))
M(x) = I - Y(x)*jacobian(f, x)
K(x) = make_intervals(mid(x)) - Y(x)*make_intervals(f(mid(x))) + M(x)*(x - make_intervals(mid(x)))

function krawczyk2d_internal(f, a::Array{Interval, 1}, bigprec::Integer)



	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	#if mid(a) == 0
	#	a = Interval(a.lo, a.hi + 0.0001*mag(a))
	#end

  Ka = isect(a, K(a))
	if Ka != false
		#@show a
		#@show diam(a)

		d = diam(a)
		dK = diam(Ka)

		if d[1] < tol && d[2] < tol #d == dK
			if dK[1] < d[1] && dK[2] < d[2] # Put <= instead of <
				println("Unique zero in $Ka")
				push!(arr_a, Ka)
			else
				println("Maybe a zero in $Ka")
			end

		else
			M1 = Y(a)
			if mid(det2(M(Ka))) <= mid(det2(M(a)))
				krawczyk2d_internal(f, Ka, bigprec)

			else
				krawczyk2d_internal(f, make_intervals(mid(a)) - M1*make_intervals(f(mid(a))) + (I - M1*jacobian(f, a))*(a - make_intervals(mid(a))), bigprec)
			end

		end

	end

	return arr_a
end

return krawczyk2d_internal(f, a, bigprec)
end

