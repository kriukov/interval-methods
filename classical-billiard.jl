# In all configuration it is assumed that radius of the disk R = 1
using IntervalArithmetic
using PurityIntervals

#hit(init.pos., init.vel., dist.betw.centers) = final pos., final vel., time
function hit(x, v, r)
    v /= norm(v)
    c12 = [r, 0]
    c13 = r/2*[1, sqrt(3)]
	B12 = dot(x - c12, v)
	B13  = dot(x - c13, v)
	C12 = r*(r - 2x[1])
	C13 = r*(r - x[1] - sqrt(3)*x[2])
	if B12^2 - C12 >= 0
	    #hits second disk
	    t1 = -B12 - sqrt(B12^2 - C12)
	elseif B13^2 - C13 >= 0
	    #hits third disk
	    t1 = -B13 - sqrt(B13^2 - C13)
	else
	    #misses and flies off to infinity
	    t1 = Inf
	end
	if t1 < Inf # If collides
	    N0 = x + v*t1
	    N = N0/norm(N0)
	    # Velocity, place and time immediately after the collision
	    v1 = v - 2dot(v, N)*N
	    x1 = x + v*t1
	    #t += t1
	    return x1, v1, t1
	else
	    return [Inf, Inf], v, Inf
	end
end

#theta = pi/6
#x = [cos(theta), sin(theta)]
#v = [1, 0]
# Hit 2nd disk 
#hit2 = hit(x, v, 6)
# Hit 3rd disk
#hit3 = hit(x, [0.5, sqrt(3)/2], 6)
# Miss
#miss = hit(x, [1, 1], 6)

###################################################


function trajectory(x, v, c, maxsteps, prec=64)
    set_bigfloat_precision(prec)
    #=places = MultiDimInterval[]
    disks = Int[]
    speeds = MultiDimInterval[]
    times = Interval[]
    =#
    places = Any[]
    disks = Any[]
    speeds = Any[]
    times = Any[]
    
    n = length(c)
    step = 1
    t0 = Interval(0)

    while step <= maxsteps
        B = Union{Interval, PurityInterval}[]
        C = Union{Interval, PurityInterval}[]
        t = PurityInterval[]

        for i = 1:n
            #push!(B, dot(x - c[i], v))
            push!(B, (x[1] - c[i][1])*v[1] + (x[2] - c[i][2])*v[2])
            push!(C, normsq(x - c[i]) - 1)
            push!(t, (-B[i] - sqrt(PurityInterval(B[i]^2 - normsq(v)*C[i], 1))/normsq(v)))
        end

        # Extract min. pos. time value of clean times
        @show t
        t_minpos = Interval(Inf, Inf)
        disk = 0
        for i = 1:n
            if t[i].flag == 1
	            if mid(t[i].interval) < mid(t_minpos) && t[i].interval.lo >= 0
		            t_minpos = t[i].interval
		            disk = i
	            end
	        end
        end
        
        if disk == 0
            println("Diverged after step $step")
            break
        end

        # New position and velocity (the moment after collision)

        #x1 = x + v*t_minpos
        x1 = x + [v[1]*t_minpos, v[2]*t_minpos]
        N0 = x1 - c[disk]
        #N = N0/Interval(norm(N0))
        #@show typeof(N0)
        #N = N0/sqrt(N0[1]^2 + N0[2]^2)
        N = [N0[1]/sqrt(N0[1]^2 + N0[2]^2), N0[2]/sqrt(N0[1]^2 + N0[2]^2)]
        #v1 = v - 2dot(v, N)*N
        @show v1 = v - [2*(v[1]*N[1] + v[2]*N[2])*N[1], 2*(v[1]*N[1] + v[2]*N[2])*N[2]]
        @show t0 += t_minpos

        #= push!(places, x1)
        push!(disks, disk)
        push!(speeds, v1)
        push!(times, t0)
        =#
        x = x1
        v = v1
        step += 1
    end
    return x #places, disks, speeds, times
end


# Example: symmetric 3-disk scatterer with r = 6
r = 6
c = MultiDimInterval[]
push!(c, [Interval(0), Interval(0)], [Interval(r), Interval(0)], [Interval(r/2), Interval(r*sqrt(3)/2)])

# Initial conditions - not necessarily on disk curcumference
x = [Interval(r/2), Interval(0.5)]
v = [Interval(1), Interval(0)]

#traj = trajectory(x, v, c, 2)

f(x) = trajectory(x, v, c, 1)

purity(f, x)

