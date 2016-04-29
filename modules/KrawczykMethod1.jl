## Krawczyk method, based on Tucker p. 86

module KrawczykMethod1
export krawczyk, differentiate, Interval, rad, diam, mid, mig, mag, belong, hd, hull, isect, isectext, K, arcsin, sqrt1

using IntervalArithmetic
using AutoDiff

println("Syntax: krawczyk(function, Interval(lo, hi), precision [default is 64])")

# If the derivative is 0, the constant C will return a "division by thin zero" error. We slightly change the denominator of the constant if this is the case.
function C(f, x)
	if differentiate(f, mid(x)) == 0
		return Interval(1)//Interval(differentiate(f, mid(x)) + 0.0001*rad(x))
	else
		return Interval(1)//Interval(differentiate(f, mid(x)))
	end
end

K(f, x) = mid(x) - C(f, x)*f(mid(x)) + (1 - C(f, x)*differentiate(f, x))*(x - mid(x))

function krawczyk(f::Function, a::Interval, bigprec::Integer=64)

    arr_a = Interval[]
    steps = 0
    a1 = Interval(Inf)
    set_bigfloat_precision(bigprec)

    if mid(a) == 0
      a = Interval(a.lo, a.hi + 0.0001*mag(a))
    end

    push!(arr_a, a)
    Ka = isect(a, K(f, a))
    @show a
    @show Ka
    
    #while a1 != a
        if Ka != false

	      if Ka == a && diam(Ka) < 1e-4
	        if inside(Ka, a)
	          println("Unique zero in $(Ka)")
	          push!(arr_a, [Ka, :unique])
	        else
	          println("Maybe a zero in $(Ka)")
	          push!(arr_a, [Ka, :possible])
	        end

	      else
            @show steps += 1
	        @show krawczyk_internal(f, left(Ka), bigprec)
	        @show krawczyk_internal(f, right(Ka), bigprec)
	        
	      end

	    end
        @show a1 = a
    #end
    return arr_a
end

#end of module
end
