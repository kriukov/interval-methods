include("modules/interval.jl")
include("modules/ad-int.jl")

center(x) = Interval(mid(x))
K(x) = mid(x) - f(mid(x))/differentiate(f, center(x)) - (Interval(1.) - differentiate(f, center(x))/differentiate(f, center(x)))*Interval(-rad(x), rad(x))

x = ad(a, Interval(1.))
do_isect(x, y) = isectext(arr_a[x], arr_b[y])

	arr_a = Interval[]
	arr_b = Interval[]

	push!(arr_a, a)

	k = 0
	while k <= calls

		arr_b = Interval[]

		for i = 1:length(arr_a)
			push!(arr_b, K(arr_a[i]))
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
