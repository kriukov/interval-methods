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

# Functions for bisection

#=
function lower(x::MultiDimInterval)
	z = Interval[]
	for i = 1:length(x)
		push!(z, Interval(x[i].lo, mid(x[i])))
	end
	return z
end

function upper(x::MultiDimInterval)
	z = Interval[]
	for i = 1:length(x)
		push!(z, Interval(mid(x[i]), x[i].hi))
	end
	return z
end
=#


left(x::Interval) = Interval(x.lo, mid(x))
right(x::Interval) = Interval(mid(x), x.hi)


function bisect(xx::Vector{Interval})
 
	if length(xx) != 2
	error("Only works for 2 at the moment")
	end

	x, y = xx

	@show x

	intervals = Vector{Interval}[]

	push!(intervals, [left(x), left(y)])
	push!(intervals, [left(x), right(y)])
	push!(intervals, [right(x), left(y)])
	push!(intervals, [right(x), right(y)])

	intervals
end

function krawczyk2d(f, a::MultiDimInterval, bigprec::Integer=64)

arr_a = MultiDimInterval[]

set_bigfloat_precision(bigprec)
#tol = (1e-10)*eps(BigFloat)
tol = 1e-10

I = [Interval(1) Interval(0); Interval(0) Interval(1)]
Y(x) = make_intervals(inv(mid(jacobian(f, mid(x)))))
M(x) = I - Y(x)*jacobian(f, x)
K(x) = make_intervals(mid(x)) - Y(x)*make_intervals(f(mid(x))) + M(x)*(x - make_intervals(mid(x)))

Y1 = Array{Interval, 2}[]
push!(Y1, Y(a))

Kprev(x) = make_intervals(mid(x)) - Y1[k-1]*make_intervals(f(mid(x))) + (I - Y1[k-1]*jacobian(f, x))*(x - make_intervals(mid(x)))


k = 1

function krawczyk2d_internal(f, a::MultiDimInterval, bigprec::Integer)

	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	#if mid(a) == 0
	#	a = Interval(a.lo, a.hi + 0.0001*mag(a))
	#end

  Ka = isectext(a, K(a))
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
					@show krawczyk2d_internal(f, bisect(Ka)[1], bigprec)
					@show krawczyk2d_internal(f, bisect(Ka)[2], bigprec)
					@show krawczyk2d_internal(f, bisect(Ka)[3], bigprec)
					@show krawczyk2d_internal(f, bisect(Ka)[4], bigprec)

				else
					k += 1
					push!(Y1, Y(a))
					@show krawczyk2d_internal(f, bisect(Kprev(a))[1], bigprec)
					@show krawczyk2d_internal(f, bisect(Kprev(a))[2], bigprec)
					@show krawczyk2d_internal(f, bisect(Kprev(a))[3], bigprec)
					@show krawczyk2d_internal(f, bisect(Kprev(a))[4], bigprec)					

				end
		
			else				
				
		

				if mid(det2(I - Y1[k]*jacobian(f, a))) <= mid(det2(I - Y1[k-1]*jacobian(f, a)))
					@show k += 1
					push!(Y1, Y(a))
					@show krawczyk2d_internal(f, bisect(Ka)[1], bigprec)
					@show krawczyk2d_internal(f, bisect(Ka)[2], bigprec)
					@show krawczyk2d_internal(f, bisect(Ka)[3], bigprec)
					@show krawczyk2d_internal(f, bisect(Ka)[4], bigprec)
					
				else
					@show k += 1
					push!(Y1, Y(a))
					@show krawczyk2d_internal(f, bisect(Kprev(a))[1], bigprec)
					@show krawczyk2d_internal(f, bisect(Kprev(a))[2], bigprec)
					@show krawczyk2d_internal(f, bisect(Kprev(a))[3], bigprec)
					@show krawczyk2d_internal(f, bisect(Kprev(a))[4], bigprec)	

				end
			
			
			end
			
			
		end

	end

	return arr_a
end

return krawczyk2d_internal(f, a, bigprec)
end

