## Krawczyk method, based on Tucker p. 86

module KrawczykMethod
export krawczyk, differentiate, Interval, rad, diam, mid, mig, mag, belong, hd, hull, isect, isectext

using AutoDiff

println("Syntax: krawczyk(function, Interval(lo, hi), precision)")

# Outside wrapping function was made in order for it to clean up the array of roots every time
function krawczyk(f::Function, a::Interval, bigprec::Integer)


arr_a = Interval[]

function krawczyk_internal(f::Function, a::Interval, bigprec::Integer)

	set_bigfloat_precision(bigprec)
	tol = 100eps(BigFloat)
	
	K(x) = mid(x) - (Interval(1)//differentiate(f, mid(x)))*f(mid(x)) + (1 - (Interval(1)//differentiate(f, mid(x)))*differentiate(f, x))*(x - mid(x))
	
	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	if mid(a) == 0
		a = Interval(a.lo, a.hi + 0.0001*mag(a))
	end
	
	if isect(a, K(a)) != false
		#@show diam(a)
		#@show diam(K(a))


		if diam(K(a)) < tol
			if diam(K(a)) < diam(a)
				println("Unique zero in $(K(a))")
				push!(arr_a, K(a))
			else
				println("Maybe a zero in $(K(a))")
			end
		
		else
			krawczyk_internal(f, Interval(a.lo, mid(a)), bigprec)
			krawczyk_internal(f, Interval(mid(a), a.hi), bigprec)			
		end

	end
	
	return arr_a
end


return krawczyk_internal(f, a, bigprec)
end

#end of module
end
