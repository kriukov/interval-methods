module KrawczykMethod2D
export krawczyk2d, differentiate, Interval, rad, diam, mid, mig, mag, lo, hi, belong, hd, hull, isect, isectext, K, all_inside, bisect, jacobian, MultiDimInterval, Y, krawczyk2d_internal, mod21, Ad, mod1, mod2, mod21, mod22, mod23, mod24, arcsin, sqrt1, krawczyk2d_general

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

# Check if x is within the range of function f
function is_allowed(f, x)
    try f(x)
    catch DomainError
        return false
    end
    return true
end

function krawczyk2d(f, a::MultiDimInterval, bigprec::Integer=64)

    roots_array = MultiDimInterval[]

    set_bigfloat_precision(bigprec)
    #tol = (1e-10)*eps(BigFloat)
    tol = 1e-10

    I = [Interval(1) Interval(0); Interval(0) Interval(1)]

    # If the jacobian is non-invertible, the SingularException error is returned for Y. We need to choose a slightly different Y then.	
    function Y(f, x)
	    if det(jacobian(f, mid(x))) == 0
		    return make_intervals(inv(jacobian(f, mid(x) + 0.0001*norm(diam(x)))))
	    else
		    return make_intervals(inv(jacobian(f, mid(x))))
	    end	
    end
    #Y(f, x) = make_intervals(inv(jacobian(f, mid(x))))

    M(f, x) = I - Y(f, x)*jacobian(f, x)
    K(f, x) = make_intervals(mid(x)) - Y(f, x)*make_intervals(f(mid(x))) + M(f, x)*(x - make_intervals(mid(x)))

    k = 1

    function krawczyk2d_internal(f, a::MultiDimInterval, bigprec::Integer)

	    Ka = isect(a, K(f, a))
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
		    end

	    end

	    return roots_array
    end

    return krawczyk2d_internal(f, a, bigprec)
end

function krawczyk2d_general(f, a::MultiDimInterval, prec::Integer=64)
    tol = 1e-10
    while max(diam(a)[1], diam(a)[2]) > tol
        if is_allowed(f, a) == true
            return krawczyk2d(f, a, prec)
        else
            @show pieces = bisect(a)
            krawczyk2d_general(f, pieces[1], prec)
            krawczyk2d_general(f, pieces[2], prec)
            krawczyk2d_general(f, pieces[3], prec)
            krawczyk2d_general(f, pieces[4], prec)
        end
    end

end

#end of module
end

