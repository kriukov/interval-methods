using IntervalArithmetic
using AutoDiff

function all_inside(x::Array{Interval, 1}, y::Array{Interval, 1})
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


function krawczyk2d(f, a::Array{Interval, 1}, bigprec::Integer=64)

arr_a = Array{Interval, 1}[]

set_bigfloat_precision(bigprec)
#tol = (1e-10)*eps(BigFloat)
tol = 1e-10

I = [Interval(1) Interval(0); Interval(0) Interval(1)]
Y(x) = make_intervals(inv(mid(jacobian(f, mid(x)))))
M(x) = I - Y(x)*jacobian(f, x)
K(x) = make_intervals(mid(x)) - Y(x)*make_intervals(f(mid(x))) + M(x)*(x - make_intervals(mid(x)))

Y1 = Array{Interval, 2}[]
push!(Y1, Y(a))

k = 1

function krawczyk2d_internal(f, a::Array{Interval, 1}, bigprec::Integer)

	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	#if mid(a) == 0
	#	a = Interval(a.lo, a.hi + 0.0001*mag(a))
	#end

  Ka = isect(a, K(a))
	if Ka != false
		@show a
		@show diam(a)

		d = diam(a)
		dK = diam(Ka)

		if dK[1] < tol && dK[2] < tol #d == dK
			if all_inside(Ka, a)
				println("Unique zero in $Ka")
				push!(arr_a, Ka)
			else
				println("Maybe a zero in $Ka")
			end
			k += 1
		else		
		
		
			if k == 1
		
				if true #@show mid(det2(I - Y1[k]*jacobian(f, a))) <= mid(det2(I - Y(a)*jacobian(f, a)))
					k += 1
					push!(Y1, Y(a))
					krawczyk2d_internal(f, Ka, bigprec)

				else
					k += 1
					push!(Y1, Y(a))
					krawczyk2d_internal(f, make_intervals(mid(a)) - Y1[k-1]*make_intervals(f(mid(a))) + (I - Y1[k-1]*jacobian(f, a))*(a - make_intervals(mid(a))), bigprec)

				end
		
			else				
				
		

				if mid(det2(I - Y1[k]*jacobian(f, a))) <= mid(det2(I - Y1[k-1]*jacobian(f, a)))
					@show k += 1
					push!(Y1, Y(a))
					krawczyk2d_internal(f, Ka, bigprec)

				else
					@show k += 1
					push!(Y1, Y(a))
					krawczyk2d_internal(f, make_intervals(mid(a)) - Y1[k-1]*make_intervals(f(mid(a))) + (I - Y1[k-1]*jacobian(f, a))*(a - make_intervals(mid(a))), bigprec)

				end
			
			
			end
			
			
		end

	end

	return arr_a
end

return krawczyk2d_internal(f, a, bigprec)
end

