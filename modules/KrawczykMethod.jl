## Krawczyk method, based on Tucker p. 86

module KrawczykMethod
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

# Outside wrapping function was made in order for it to clean up the array of roots every time
function krawczyk(f::Function, a::Interval, bigprec::Integer=64)

  arr_a = Any[] # instead of Interval[] to put symbols at each encountered root
  steps = 0
  a1 = Interval(Inf)
  function krawczyk_internal(f::Function, a::Interval, bigprec::Integer)

	set_bigfloat_precision(bigprec)
	tol = 1e-15 #100eps(BigFloat)

	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then. This fix also removes the "possible" roots which correctly merge into "unique".
	if mid(a) == 0
	  a = Interval(a.lo, a.hi + 0.0001*mag(a))
	end

	Ka = isect(a, K(f, a))
	#@show a
	#@show Ka
    
	if Ka != false
	  #@show diam(a) > diam(Ka)
	  #@show a1
	  #@show diam(K(a))

	  if diam(Ka) < tol #1e-4 && (Ka == a)
	    if inside(Ka, a)   # inside
	      println("Unique zero in $(Ka)")
	      #if length(arr_a) > 0 && Ka != arr_a[length(arr_a)]
	      push!(arr_a, [Ka, :unique])
	      #end
	    else
	      println("Maybe a zero in $(Ka)")
	      push!(arr_a, [Ka, :possible])
	    end

	  else
        #@show steps += 1
        #if steps > 100
        #    error("Exiting")
        #end
        #@show a1 = a
	    krawczyk_internal(f, left(Ka), bigprec)
	    krawczyk_internal(f, right(Ka), bigprec)
	    
	  end

	end
    
	return arr_a
  end
  
  #Weeding out duplicate and wrong roots
  #=
  @show sol = krawczyk_internal(f, a, bigprec)
  sol1 = Interval[]
  for i=1:length(sol)
     push!(sol1, sol[i][1])
  end
  sol2 = unique(sol1)
  sol3 = Interval[]
  for i=1:length(sol2)
    if belong(0, f(sol2[i]))
        push!(sol3, sol2[i])
    end
  end
  return sol3
  =#
  return krawczyk_internal(f, a, bigprec)
end

#end of module
end
