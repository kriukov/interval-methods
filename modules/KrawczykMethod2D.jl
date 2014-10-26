module KrawczykMethod2D
export krawczyk2d, differentiate, Interval, rad, diam, mid, mig, mag, belong, hd, hull, isect, isectext, K, all_inside, bisect, jacobian

using IntervalArithmetic
using AutoDiff

println("Syntax: krawczyk2d(function, [Interval(lo, hi), Interval(lo, hi)], precision [default is 64])")

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

function bisect(xx::MultiDimInterval)
 	if length(xx) != 2
	error("Only works for 2 at the moment")
	end

	x, y = xx

	intervals = MultiDimInterval[]

	push!(intervals, [left(x), left(y)])
	push!(intervals, [left(x), right(y)])
	push!(intervals, [right(x), left(y)])
	push!(intervals, [right(x), right(y)])

	intervals
end

function krawczyk2d(f, a::MultiDimInterval, bigprec::Integer=64)

roots_array = MultiDimInterval[]

set_bigfloat_precision(bigprec)
#tol = (1e-10)*eps(BigFloat)
tol = 1e-10

I = [Interval(1) Interval(0); Interval(0) Interval(1)]
Y(f, x) = make_intervals(inv(mid(jacobian(f, mid(x)))))
M(f, x) = I - Y(f, x)*jacobian(f, x)
K(f, x) = make_intervals(mid(x)) - Y(f, x)*make_intervals(f(mid(x))) + M(f, x)*(x - make_intervals(mid(x)))

#Y1 = Array{Interval, 2}[]
#push!(Y1, Y(a))

#Kprev(x) = make_intervals(mid(x)) - Y1[k-1]*make_intervals(f(mid(x))) + (I - Y1[k-1]*jacobian(f, x))*(x - make_intervals(mid(x)))


k = 1

function krawczyk2d_internal(f, a::MultiDimInterval, bigprec::Integer)

	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval array should be slightly asymmetrized then
	#if mid(a) == [0, 0]
	#	a = [Interval(a[1].lo, a[1].hi + 0.0001*mag(a[1])), Interval(a[2].lo, a[2].hi + 0.0001*mag(a[2]))]
	#end

	Ka = isectext(a, K(f, a))
	if Ka != false

		d = diam(a)
		dK = diam(Ka)

		if dK[1] < tol && dK[2] < tol #d == dK
			if all_inside(Ka, a)
				println("Unique zero in $Ka")
				push!(roots_array, Ka)
			else
				println("Maybe a zero in $Ka")
			end
			k += 1
		else		
			k += 1
			krawczyk2d_internal(f, bisect(Ka)[1], bigprec)
			krawczyk2d_internal(f, bisect(Ka)[2], bigprec)
			krawczyk2d_internal(f, bisect(Ka)[3], bigprec)
			krawczyk2d_internal(f, bisect(Ka)[4], bigprec)
			
			# Code below didn't give any gain in time
			#@show krawczyk2d_internal(f, [left(Ka[1]), left(Ka[2])], bigprec)
			#@show krawczyk2d_internal(f, [left(Ka[1]), right(Ka[2])], bigprec)
			#@show krawczyk2d_internal(f, [right(Ka[1]), left(Ka[2])], bigprec)
			#@show krawczyk2d_internal(f, [right(Ka[1]), right(Ka[2])], bigprec)
		end

	end

	return roots_array
end

return krawczyk2d_internal(f, a, bigprec)
end

#end of module
end

