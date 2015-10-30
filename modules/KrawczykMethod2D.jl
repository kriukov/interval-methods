module KrawczykMethod2D
export krawczyk2d, differentiate, Interval, rad, diam, mid, mig, mag, lo, hi, belong, hd, hull, isect, isectext, K, all_inside, left, right, bisect, jacobian, MultiDimInterval, make_intervals, Y, krawczyk2d_internal, mod21, Ad, mod1, mod2, mod21, mod22, mod23, mod24, arcsin, sqrt1, krawczyk2d_purity, domaincheck, domaincheck2d, arcsin_d, sqrt_d, krawczyk2d_loose, purity

using IntervalArithmetic
using AutoDiff
using PurityIntervals

println("Syntax: krawczyk2d(function, [Interval(lo, hi), Interval(lo, hi)], precision [default is 64])")

function krawczyk2d(f, a::MultiDimInterval, bigprec::Integer=64)

    roots_array = MultiDimInterval[]

    set_bigfloat_precision(bigprec)

    tol = 1e-10

    I = [Interval(1) Interval(0); Interval(0) Interval(1)]

    # If the jacobian is non-invertible, the SingularException error is returned for Y. We need to choose a slightly different Y then.
    function Y(f, x)
        midx = mid(x)
	    if det(jacobian(f, midx)) == 0
		    return make_intervals(inv(jacobian(f, midx + 0.0001*norm(diam(x)))))
	    else
		    return make_intervals(inv(jacobian(f, midx)))
	    end
    end
    #Y(f, x) = make_intervals(inv(jacobian(f, mid(x))))

    M(f, x) = I - Y(f, x)*jacobian(f, x)

    function K(f, x)
        midx = mid(x)
        intmidx = make_intervals(midx)
        intmidx - Y(f, x)*make_intervals(f(midx)) + M(f, x)*(x - intmidx)
    end
    #K(f, x) = make_intervals(mid(x)) - Y(f, x)*make_intervals(f(mid(x))) + M(f, x)*(x - make_intervals(mid(x)))

    k = 1

    function krawczyk2d_internal(f, a::MultiDimInterval, bigprec::Integer)

	    Ka = isect(a, K(f, a))

	    if Ka != false

		    dK = diam(Ka)

		    if dK[1] < tol && dK[2] < tol #d == dK

			    if all_inside(Ka, a) #&& all_inside(f(Ka), [Interval(-tol, tol), Interval(-tol, tol)])
				    println("Unique zero in $Ka")
				    push!(roots_array, Ka)
			    else
				    println("Maybe a zero in $Ka")
				    push!(roots_array, Ka)
			    end
			    k += 1
		    else
			    k += 1
			    pieces = bisect(Ka)
			    for i = 1:4
			        krawczyk2d_internal(f, pieces[i], bigprec)
			    end
		    end

	    end

	    return roots_array
    end

    return krawczyk2d_internal(f, a, bigprec)
end

# Version of krawczyk2d with purity
function krawczyk2d_purity(f, a::MultiDimInterval, prec::Integer=64, tol=1e-4)
    #sol = Array{Array{Interval, 1}, 1}[]

    roots_array = MultiDimInterval[]

    set_bigfloat_precision(prec)

    I = [Interval(1) Interval(0); Interval(0) Interval(1)]

    intdet(A) = A[1]*A[4] - A[2]*A[3]

    # If the jacobian is non-invertible, the SingularException error is returned for Y. We need to choose a slightly different Y then.
    function Y(f, x)
        midx = mid(x)
        D = intdet(jacobian(f, midx))
	    if D == Interval(0) || D == Interval(Inf) || D == Interval(-Inf)
		    return make_intervals(inv(jacobian(f, mid(x + [Interval(0.0001), Interval(0.0002)]))))
	    else
		    return make_intervals(inv(jacobian(f, midx)))
	    end
    end

    M(f, x) = I - Y(f, x)*jacobian(f, x)

    # Krawczyk operator
    function K(f, x)
        midx = mid(x)
        intmidx = make_intervals(midx)
        intmidx - Y(f, x)*make_intervals(f(midx)) + M(f, x)*(x - intmidx)
    end

    function krawczyk2d_purity_internal(f, a::MultiDimInterval, prec::Integer)
        @show a
        @show p = purity(f, a)

        if p != -1
            if p == 1
                println("Clean")
                if typeof(a[1]) == IntUnion && typeof(a[2]) == Interval
                krawczyk2d(f, [a[1].elem1, a[2]], prec)
                krawczyk2d(f, [a[1].elem2, a[2]], prec)
                elseif typeof(a[1]) == Interval && typeof(a[2]) == IntUnion
                krawczyk2d(f, [a[1], a[2].elem1], prec)
                krawczyk2d(f, [a[1], a[2].elem2], prec)
                elseif typeof(a[1]) == IntUnion && typeof(a[2]) == IntUnion
                krawczyk2d(f, [a[1].elem1, a[2].elem1], prec)
                krawczyk2d(f, [a[1].elem2, a[2].elem1], prec)
                krawczyk2d(f, [a[1].elem1, a[2].elem2], prec)
                krawczyk2d(f, [a[1].elem2, a[2].elem2], prec)
                else
                krawczyk2d(f, a, prec)
                end
                
                
            elseif p == 0

                println("Unclean")
                if max(diam(a)[1], diam(a)[2]) < tol
                    #@show a
                    #if all_inside(f(a), [Interval(-tol, tol), Interval(-tol, tol)])
                    #    println("Unique zero in $a")
			        #    push!(roots_array, a)
			        #end
			        println("Krawczyk could not properly identify zeros here: tolerance reached at unclean")
			    else
                    pieces = bisect(a)
                    for i = 1:4
                        krawczyk2d_purity_internal(f, pieces[i], prec)
                    end
                end
            end
        end

        return roots_array
    end
    return krawczyk2d_purity_internal(f, a, prec)
end




##### Deprecated

# Version of krawczyk2d with loose evaluation
function krawczyk2d_loose(f, a::MultiDimInterval, bigprec::Integer=64)

    roots_array = MultiDimInterval[]

    set_bigfloat_precision(bigprec)
    tol = 1e-6

    I = [Interval(1) Interval(0); Interval(0) Interval(1)]
    mdelta = [Interval(1e-5) Interval(1e-5); Interval(1e-5) Interval(1e-5)]

    intdet(M) = M[1]*M[4] - M[2]*M[3]

    # If the jacobian is non-invertible, the SingularException error is returned for Y. We need to choose a slightly different Y then.
    function Y(f, x)
	    if intdet(jacobian(f, mid(x))) == Interval(0)
		    #return make_intervals(inv(jacobian(f, mid(x) + 0.0001*norm(diam(x)))))
		    return make_intervals(inv(jacobian(f, mid(x)) + mdelta))
	    else
		    return make_intervals(inv(jacobian(f, mid(x))))
	    end
    end

    M(f, x) = I - Y(f, x)*jacobian(f, x)
    K(f, x) = make_intervals(mid(x)) - Y(f, x)*make_intervals(f(mid(x))) + M(f, x)*(x - make_intervals(mid(x)))

    k = 1
    i = 0

    function krawczyk2d_internal_loose(f, a::MultiDimInterval, bigprec::Integer)

        #Ka = isect(a, K(f, a))
        @show K1 = K(f, a)
        #if K1[1].lo == -Inf && K1[1].hi == Inf &&K1[2].lo == -Inf && K1[2].hi == Inf
        #    K1 = [Interval(-1e9, 1e9), Interval(-1e9, 1e9)]
        #end
	    @show Ka = isect(a, K1)
	    if Ka != false # !isnan(K1[1].lo) && !isnan(K1[1].hi) && !isnan(K1[2].lo) && !isnan(K1[2].hi)

		    #d = diam(a)
		    dK = diam(Ka)

		    if dK[1] < tol && dK[2] < tol  #&& i <= 20 # && Ka != a #d == dK
		        i += 1
			    if all_inside(Ka, a) && all_inside(f(Ka), [Interval(-tol, tol), Interval(-tol, tol)])

			        println("Iteration #$i")
				    println("Unique zero in $Ka")
				    @show push!(roots_array, Ka)
				    #error("Found root: $(Ka)")
			    else
			        println("Iteration #$i")
				    println("Maybe a zero in $Ka")
			    end
			    k += 1
		    else
			    k += 1
			    #@show if
			    #if (dK[1] > tol || dK[2] > tol) && (isnan(K1[1].lo) && isnan(K1[1].hi) && isnan(K1[2].lo) && isnan(K1[2].hi))
			    bisect_list = bisect(Ka)
			    krawczyk2d_internal_loose(f, bisect_list[1], bigprec)
			    krawczyk2d_internal_loose(f, bisect_list[2], bigprec)
			    krawczyk2d_internal_loose(f, bisect_list[3], bigprec)
			    krawczyk2d_internal_loose(f, bisect_list[4], bigprec)
			    #end
			    #end
		    end

	    end

	    #end

	    return roots_array
    end

    return krawczyk2d_internal_loose(f, a, bigprec)
end




#end of module
end
