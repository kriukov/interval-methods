## Newton's method (interval) code

#module NewtonMethod

#export newton, Interval

include("interval.jl")
include("ad-int.jl")

println("Syntax: newton(function, Interval(lo, hi))")

function newton(f::Function, a::Interval)

	center(x) = Interval(mid(x))
	N(x) = center(x) - f(center(x))//differentiate(f, x)
	do_isect(x, y) = isectext(arr_a[x], arr_b[y])

	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	if mid(a) == 0
		a = Interval(a.lo, a.hi + 0.0001*mag(a))
	end

	x = Ad(a, Interval(1.))

	arr_a = Interval[]
	arr_b = Interval[]

	push!(arr_a, a)

	k = 0
	while true
		
		println(k)
		@show((arr_a, arr_b))
		
		arr_b = Interval[]

		for i = 1:length(arr_a)
			push!(arr_b, N(arr_a[i]))
		end
		
		@show arr_b

		arr_a_new = Interval[]

		for i = 1:length(arr_b)
			@show do_isect(i, i)
			if do_isect(i, i) != false
				arr_a_new = vcat(arr_a_new, do_isect(i, i))
			end
		end

		if arr_a_new == arr_a 
			break
		end
				
		# # Stopping cycle routine
		# m = 0
		# # Make sure that the cycle stops when all elements of an array change by less than a ~10eps() (smaller number of eps may not finish the cycle)
		# for i = 1:length(arr_a)
		# 	# Convert true/false to 1/0
		# 	#m += ( @show mid(abs(arr_a[i] - arr_b[i])) <= 9*eps() )
		# 	@show((arr_a[i], arr_b[i]))
		# 	m += @show(arr_a[i] == arr_b[i])
		# end

		# if m >= length(arr_a)
		# 	break
		# end
		
		arr_a = arr_a_new
		k += 1
	end

	println("Function calls: ", k)
	return arr_a

end

# end of module
#end
