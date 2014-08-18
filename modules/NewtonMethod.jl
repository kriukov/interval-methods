## Newton's method (interval) code

#module NewtonMethod

#export newton, Interval

include("interval.jl")
include("ad-int.jl")

println("Syntax: newton(function, Interval(lo, hi), calls)")

function newton(f::Function, a::Interval, calls::Int32)

	center(x) = Interval(mid(x))
	N(x) = center(x) - f(center(x))//differentiate(f, x)
	do_isect(x, y) = isectext(arr_a[x], arr_b[y])

	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	if mid(a) == 0
		a = Interval(a.lo, a.hi + 0.0001*mag(a))
	end

	x = ad(a, Interval(1.))

	arr_a = Interval[]
	arr_b = Interval[]

	push!(arr_a, a)

	k = 0
	while k <= calls

		arr_b = Interval[]

		for i = 1:length(arr_a)
			push!(arr_b, N(arr_a[i]))
		end

		arr_a_new = Interval[]

		for i = 1:length(arr_b)
			if do_isect(i, i) != false
				arr_a_new = vcat(arr_a_new, do_isect(i, i))
			end
		end

		arr_a = arr_a_new
		k += 1
	end

	arr_a
end

# end of module
#end
