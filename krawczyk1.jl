using IntervalArithmetic
using AutoDiff

K(x) = mid(x) - Interval(1)/(differentiate(f, mid(x)))*f(mid(x)) + (1 - Interval(1)/(differentiate(f, mid(x)))*differentiate(f, x))*(x - mid(x))

lower(x::Interval) = Interval(x.lo, mid(x))
higher(x::Interval) = Interval(mid(x), x.hi)

	x = Ad(a, Interval(1.))

	arr_a = Interval[]
	arr_b = Interval[]

	push!(arr_a, a)

	k = 0
	while k < 20 #true
		
		println(k)
		
		arr_b = Interval[]

		for i = 1:length(arr_a)
			if isect(arr_a[i], K(arr_a[i])) != false
				if isect(arr_a[i], K(arr_a[i])) == arr_a[i]
					if isect(lower(arr_a[i]), K(lower(arr_a[i]))) != false
					push!(arr_b, isect(lower(arr_a[i]), K(lower(arr_a[i]))))
					end
					if isect(higher(arr_a[i]), K(higher(arr_a[i]))) != false
					push!(arr_b, isect(higher(arr_a[i]), K(higher(arr_a[i]))))
					end
				else push!(arr_b, isect(arr_a[i], K(arr_a[i])))
				println(arr_b)
				end
			end
		end
		
		arr_a = arr_b
		k += 1
	end
	
	
	println("Function calls: ", k)
	return arr_a
