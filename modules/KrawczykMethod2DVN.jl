module KrawczykMethod2DVN
export krawczyk2d, differentiate, Interval, IntervalBox, rad, diam, mid, mig, mag, lo, hi, belong, hd, hull, K, bisect, jacobian, Y, krawczyk2d_internal, Ad, krawczyk2d_purity, krawczyk2d_purity_periodic, purity, PurityInterval

using ValidatedNumerics
using AutoDiff
using PurityIntervalsVN

function bisect{N,T}(xx::IntervalBox{N,T})
 	if length(xx) != 2
	    error("Only works for 2 at the moment")
	end

	x, y = xx

	intervals = IntervalBox{N,T}[]

	push!(intervals, IntervalBox(x.lo..mid(x), y.lo..mid(y)))
	push!(intervals, IntervalBox(x.lo..mid(x), mid(y)..y.hi))
	push!(intervals, IntervalBox(mid(x)..x.hi, y.lo..mid(y)))
	push!(intervals, IntervalBox(mid(x)..x.hi, mid(y)..y.hi))

	intervals
end

# Make a degenerate intervalbox from an array of two numbers
make_intervals(x) = IntervalBox(x[1]..x[1], x[2]..x[2])

println("Syntax: krawczyk2d(function, [Interval(lo, hi), Interval(lo, hi)], precision [default is 64])")

function krawczyk2d{N,T}(f::Function, a::IntervalBox{N,T})

    roots_array = IntervalBox{N,T}[]

    #set_bigfloat_precision(bigprec)

    tol = 1e-10

    I = [Interval(1) Interval(0); Interval(0) Interval(1)]

    # If the jacobian is non-invertible, the SingularException error is returned for Y. We need to choose a slightly different Y then.
    function Y(f, x)
        midx = mid(x)
	    if det(AutoDiff.jacobian(f, midx)) == 0
		    return inv(AutoDiff.jacobian(f, midx + 0.0001*norm(diam(x))))
	    else
		    return inv(AutoDiff.jacobian(f, midx))
	    end
    end
    
    M(f, x) = I - Y(f, x)*AutoDiff.jacobian(f, x)

    function K(f, x)
        midx = mid(x)
        #intmidx = make_intervals(midx)
        #@show midx
        #@show f(midx)
        B = M(f, x)*([x[1], x[2]] - midx)
        make_intervals(midx - Y(f, x)*f(midx)) + IntervalBox(B[1], B[2])
    end

    k = 1

    function krawczyk2d_internal(f, a::IntervalBox)
        
	    Ka = intersect(a, K(f, a))
	    if Ka[1] != ∅ && Ka[2] != ∅
            #@show Ka
		    dK = [diam(Ka[1]), diam(Ka[2])]

		    if dK[1] < tol && dK[2] < tol #d == dK
			    if interior(Ka[1], a[1]) && interior(Ka[2], a[2])
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
			        krawczyk2d_internal(f, pieces[i])
			    end
		    end

	    end

	    return roots_array
    end


    return krawczyk2d_internal(f, a)
end


# Version of krawczyk2d with purity
function krawczyk2d_purity(f, a::IntervalBox, tol=1e-4)
    roots_array = Array{IntervalBox, 1}[]
    #set_bigfloat_precision(prec)
    
    rect_count = 0
    function krawczyk2d_purity_internal(f, a::IntervalBox)
        
        @show a
        @show p = purity(f, a)

        if p != -1
            if p == 1
                rect_count += 1
                println("Clean")
                roots = krawczyk2d(f, a)
                if length(roots) > 0
                    push!(roots_array, roots)
                    error("Roots found")
                end
                
            elseif p == 0

                println("Unclean")
                if max(diam(a[1]), diam(a[2])) < tol
			        print_with_color(:blue, "Krawczyk could not properly identify zeros here: tolerance reached at unclean\n")
			    else
                    pieces = bisect(a)
                    for i = 1:4
                        krawczyk2d_purity_internal(f, pieces[i])
                    end
                end
            end
        end

        return roots_array, rect_count
    end
    return krawczyk2d_purity_internal(f, a)
end



# Version of krawczyk2d with purity for periodic orbits (solving f(a) = a)
function krawczyk2d_purity_periodic(f, a::IntervalBox, tol=1e-4)
    roots_array = Array{IntervalBox, 1}[]
    #set_bigfloat_precision(prec)
    
    rect_count = 0
    all_count = 0
    function krawczyk2d_purity_internal_periodic(f, a::IntervalBox)
        all_count += 1
        if all_count % 10000 == 0
            println("Step $all_count: a = $a")
        end
        #@show a
        p = purity(f, a)
        #@show p
        if p != -1
            if p == 1
                rect_count += 1
                #if rect_count % 100 == 0
                #    println("Starting rectangle ", rect_count)
                #end
                #println(a, " clean")
                #@show a
                #@show f(a)
                #isect_step = intersect(a, f(a))
                #@show isect_step
                # If both parts of the IntervalBox are not empty
                #if length(isect_step) > 0 #&& isect_step[1] != ∅ && isect_step[2] != ∅
                    roots = krawczyk2d(x -> f(x) - x, a)
                    if length(roots) > 0
                        push!(roots_array, roots)
                    end
                #end
                
            elseif p == 0

                #println(a, " unclean")
                if max(diam(a[1]), diam(a[2])) < tol
			        #print_with_color(:blue, "Krawczyk could not properly identify zeros here: tolerance reached at unclean\n")
			    else
                    pieces = bisect(a)
                    for i = 1:4
                        krawczyk2d_purity_internal_periodic(f, pieces[i])
                    end
                end
            end
        end

        return roots_array, rect_count
    end
    return krawczyk2d_purity_internal_periodic(f, a)
end






#end of module
end
