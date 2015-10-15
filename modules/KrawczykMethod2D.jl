module KrawczykMethod2D
export krawczyk2d, differentiate, Interval, rad, diam, mid, mig, mag, lo, hi, belong, hd, hull, isect, isectext, K, all_inside, left, right, bisect, jacobian, MultiDimInterval, make_intervals, Y, krawczyk2d_internal, mod21, Ad, mod1, mod2, mod21, mod22, mod23, mod24, arcsin, sqrt1, krawczyk2d_general, domaincheck, domaincheck2d, arcsin_d, sqrt_d, krawczyk2d_loose

using IntervalArithmetic
using AutoDiff
using PurityIntervals

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
	    #if all_inside(a, K(f, a))
	    #    println("Krawczyk result covers initial interval")
	    #end
	    if Ka != false

		    #d = diam(a)
		    dK = diam(Ka)

		    if dK[1] < tol && dK[2] < tol #d == dK
		        
			    if all_inside(Ka, a) #&& all_inside(f(Ka), [Interval(-tol, tol), Interval(-tol, tol)])
				    println("Unique zero in $Ka")
				    push!(roots_array, Ka)
			    else
				    println("Maybe a zero in $Ka")
			    end
			    k += 1
		    else		
			    k += 1
			    bisect_list = bisect(Ka)
			    krawczyk2d_internal(f, bisect_list[1], bigprec)
			    krawczyk2d_internal(f, bisect_list[2], bigprec)
			    krawczyk2d_internal(f, bisect_list[3], bigprec)
			    krawczyk2d_internal(f, bisect_list[4], bigprec)
		    end

	    end

	    return roots_array
    end

    return krawczyk2d_internal(f, a, bigprec)
end

# Version of krawczyk2d with purity (probably deprecated)
function krawczyk2d_general(f, a::MultiDimInterval, prec::Integer=64)
sol = Array{Array{Interval, 1}, 1}[]

    function krawczyk2d_general_internal(f, a::MultiDimInterval, prec::Integer)
        
        tol = 1e-8
        if max(diam(a)[1], diam(a)[2]) > tol
            @show a
            purity = domaincheck(f, a)
            if purity == 1 # 2D!
                println("Clean")
                
                answer = krawczyk2d(f, a, prec)
                if length(answer) != 0
                    @show push!(sol, answer)
                    error("Answer added: now sol = $sol")
                else
                    println("Clean but empty")
                end
            elseif purity == 0
                println("Unclean")
                @show pieces = bisect(a)
                @show krawczyk2d_general_internal(f, pieces[1], prec)
                @show krawczyk2d_general_internal(f, pieces[2], prec)
                @show krawczyk2d_general_internal(f, pieces[3], prec)
                @show krawczyk2d_general_internal(f, pieces[4], prec)
            elseif purity == -1
                println("Dirty")
                #break
            end
        end
        return sol
    end
    return krawczyk2d_general_internal(f, a, prec)
end

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

