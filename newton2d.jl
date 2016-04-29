using IntervalArithmetic
using AutoDiff
using PurityIntervals

# Bisection

left(x::Interval) = Interval(x.lo, mid(x))
right(x::Interval) = Interval(mid(x), x.hi)


function bisect(xx::Vector{Interval})
 
	if length(xx) != 2
	error("Only works for 2 at the moment")
	end

	x, y = xx

	intervals = Vector{Interval}[]

	push!(intervals, [left(x), left(y)])
	push!(intervals, [left(x), right(y)])
	push!(intervals, [right(x), left(y)])
	push!(intervals, [right(x), right(y)])

	intervals
end

# center() makes degenerate interval midpoints of the intervals in an array
center(x) = make_intervals(mid(x))
M(f, x) = jacobian(f, x)
#N(x) = center(x) - inv(jacobian(f, x))*f(center(x))
#den(x) = M(x)[1]*M(x)[4] - M(x)[2]*M(x)[3]
#invM(x) = [M(x)[4]//den(x) -M(x)[2]//den(x); -M(x)[3]//den(x) M(x)[1]//den(x)]
invM(f, x) = inv(M(f, x))
N(f, x) = center(x) - invM(f, x)*f(center(x))

# Check if inverse jacobian has any infinite intervals
function invjacisinf(f, x)
    invM0 = invM(f, x)
    checkinf = 0
    for i = 1:length(invM0)
        checkinf += abs(diam(invM0[i]))
    end
    checkinf == Inf
end

function newton2d(f, a::MultiDimInterval, bigprec::Integer=64)

	set_bigfloat_precision(bigprec)
	
	# If a is symmetric, i.e., mid(a) = 0, the process may stall. The initial interval should be slightly asymmetrized then
	#if mid(a) == 0
	#	a = Interval(a.lo, a.hi + 0.0001*mag(a))
	#end

	roots_array = MultiDimInterval[]

	push!(roots_array, a)

	k = 0
	while true
	    
		roots_array_new = MultiDimInterval[]
		@show invjacisinf(f, a)

		for i = 1:length(roots_array)
			if isectext(roots_array[i], N(f, roots_array[i])) != false
				if !invjacisinf(f, a)
					@show roots_array_new = push!(roots_array_new, isectext(roots_array[i], N(f, roots_array[i])))
				else
					for j = 1:4
					    #@show push!(roots_array_new, isectext(bisect(roots_array[i])[j], N(f, bisect(roots_array[i])[j])))
					    @show roots_array[i]
					    if diam(roots_array[i][1]) > 1e-20 && diam(roots_array[i][2]) > 1e-20
					    newton2d(f, bisect(roots_array[i])[j])
					    end
					end
				end
			end
		end

		# Exit criterion
		if roots_array_new == roots_array 
			break
		end
		
		roots_array = roots_array_new
		k += 1
	end

	println("Function calls: ", k)
	return roots_array
end

function newton2d_purity(f, a::MultiDimInterval, bigprec::Integer=64, tol=1e-4)
    set_bigfloat_precision(bigprec)
    roots_array = Array{MultiDimInterval, 1}[]
    p = purity(f, a)
    
    if p != -1
        if p == 1
            println("Clean")
            roots = newton2d(f, a, bigprec)
            if length(roots) > 0
                push!(roots_array, roots)
            end
            
        elseif p == 0

            println("Unclean")
            if max(diam(a)[1], diam(a)[2]) < tol
		        print_with_color(:blue, "Newton could not properly identify zeros here: tolerance reached at unclean\n")
		    else
                pieces = bisect(a)
                for i = 1:4
                    newton2d_purity(f, pieces[i], bigprec)
                end
            end
        end
    end
    return roots_array
end



