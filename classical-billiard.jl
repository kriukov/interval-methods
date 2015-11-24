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
    places = MultiDimInterval[]
    disks = Int[]
    speeds = MultiDimInterval[]
    times = Interval[]
    #=
    places = Any[]
    disks = Any[]
    speeds = Any[]
    times = Any[]
    =#
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

        push!(places, x1)
        push!(disks, disk)
        push!(speeds, v1)
        push!(times, t0)

        x = x1
        v = v1
        step += 1
    end
    return places, disks, speeds, times
end


function collision_purity(x, v, c, tol, prec=64)
    set_bigfloat_precision(prec)
    n = length(c)
    step = 1
    t0 = Interval(0)
    rectangles = MultiDimInterval[]
    disks = Int[]
    
    function collision_purity_internal(x, v, c)
        B = Union{Interval, PurityInterval}[]
        C = Union{Interval, PurityInterval}[]
        t = PurityInterval[]

        for i = 1:n
            #push!(B, dot(x - c[i], v))
            push!(B, (x[1] - c[i][1])*v[1] + (x[2] - c[i][2])*v[2])
            push!(C, normsq(x - c[i]) - 1)
            push!(t, (-B[i] - sqrt(PurityInterval(B[i]^2 - normsq(v)*C[i], 1))/normsq(v)))
        end
        
        # We obtained an array of times to disks with the corresponding purities. Regardless of the purity, if the interval is fully negative, it means that the corresponding disk is behind.
        
        # Finding time intervals with p = 1
        k = 0
        for i = 1:n
            if t[i].flag == 1 && t[i].interval.lo >= 0
                k = i
                break
            end    
        end
        
        # If found such an interval, record it and the disk which is hit with certainty from this interval. If not, bisect the original rectangle and repeat the function
        if k > 0
            @show push!(rectangles, x)
            @show push!(disks, k)
        elseif k == 0
            if max(diam(x)[1], diam(x)[2]) > tol
                pieces = bisect(x)
                for i = 1:4
                    collision_purity_internal(pieces[i], v, c)
                end
            else
                println("Tolerance reached")
            end
        end
        return rectangles, disks
        end
    
    return collision_purity_internal(x, v, c)
end

function collision_purity_flat(x, v, c, tol, prec=64)
    set_bigfloat_precision(prec)
    n = length(c)
    t0 = Interval(0)
    rectangles = MultiDimInterval[]
    times = Interval[]
    disks = Int[]
    purities = Int[]
    N1 = floor(1/tol)
    dx = diam(x[1]/N1)
    dy = diam(x[2]/N1)
    for i0 = 1:N1
        for j0 = 1:N1
            rect = [x[1].lo + i0*dx - Interval(0, dx), x[2].lo + j0*dy - Interval(0, dy)]
            #B = Union{Interval, PurityInterval}[]
            #C = Union{Interval, PurityInterval}[]
            X = MultiDimInterval[]
            t = PurityInterval[]

            for i = 1:n
                #push!(B, dot(rect - c[i], v))
                #push!(B, (rect[1] - c[i][1])*v[1] + (rect[2] - c[i][2])*v[2])
                #push!(C, normsq(rect - c[i]) - 1)
                push!(X, rect - c[i])
                push!(t, (-(X[i][1]*v[1] + X[i][2]*v[2]) - sqrt(PurityInterval(normsq(v) - (X[i][1]*v[2] - X[i][2]*v[1])^2, 1))/normsq(v)))
            end
            
            #println(rect, t)
            
            # Finding valid time intervals
            k = 0; p = -2; t1 = 0
            for i = 1:n
                if mid(t[i].interval) >= 0 # && t[i].flag == 1 
                    k = i; p = t[i].flag; t1 = t[i].interval
                    break
                end    
            end
            
            if k > 0
                push!(rectangles, rect)
                push!(disks, k)
                push!(times, t1)
                push!(purities, p)
            end    
            
        end    
    end
    return rectangles, disks, times, purities
end


# Example: symmetric 3-disk scatterer with r = 6
r = 6
c = MultiDimInterval[]
push!(c, [Interval(0), Interval(0)], [Interval(r), Interval(0)], [Interval(r/2), Interval(r*sqrt(3)/2)])

#= Initial conditions - not necessarily on disk circumference
x = [Interval(r/2), Interval(0.5)]
v = [Interval(1), Interval(0)]

traj(steps) = trajectory(x, v, c, steps)

f(x) = trajectory(x, v, c, 1)[1][1]

purity(f, x)
=#

x = [Interval(-2, 7), Interval(-2, 7)]
v = [Interval(1), Interval(0)]
#collision_purity_flat(x, v, c, 1e-1)

# Function to draw rectangles
using PyPlot
using PyCall
@pyimport matplotlib.patches as patches
rectangle = patches.Rectangle
function draw_rectangle(x, y, xwidth, ywidth, color="grey")
    ax = gca()
    ax[:add_patch](rectangle((x, y), xwidth, ywidth, facecolor=color, alpha=0.5))
end


#= Plotting rectangles from which second ball will be hit (clean - green, unclean - yellow) and third ball (clean - blue, unclean - brown)
rects, disks, times, purities = collision_purity_flat(x, v, c, 0.01)
for i = 1:length(rects)
    #if disks[i] == 1
        #draw_rectangle(rects[i][1].lo, rects[i][2].lo, diam(rects[i][1]), diam(rects[i][2]), "red")
    if disks[i] == 2 && purities[i] == 1
        draw_rectangle(rects[i][1].lo, rects[i][2].lo, diam(rects[i][1]), diam(rects[i][2]), "green")
    elseif disks[i] == 2 && purities[i] == 0
        draw_rectangle(rects[i][1].lo, rects[i][2].lo, diam(rects[i][1]), diam(rects[i][2]), "yellow")
    elseif disks[i] == 3 && purities[i] == 1
        draw_rectangle(rects[i][1].lo, rects[i][2].lo, diam(rects[i][1]), diam(rects[i][2]), "blue")
    elseif disks[i] == 3 && purities[i] == 0
        draw_rectangle(rects[i][1].lo, rects[i][2].lo, diam(rects[i][1]), diam(rects[i][2]), "brown")
    
    end
end
axis([-2, 7, -2, 7])
=#

# 2 steps: take rectangle and, if there is a collision, try to find next

rects, disks, times, purities = collision_purity_flat(x, v, c, 0.01)

collision_rectangles = MultiDimInterval[]
new_velocities = MultiDimInterval[]
normals = MultiDimInterval[]
normalized_normals = MultiDimInterval[]
new_times = Interval[]
new_disks = Int[]
new_purities = Int[]

for i = 1:length(rects)
	@show push!(collision_rectangles, rects[i] + [v[1]*times[i], v[2]*times[i]])
	push!(normals, collision_rectangles[i] - c[disks[i]])
	push!(normalized_normals, [normals[i][1]/sqrt(normals[i][1]^2 + normals[i][2]^2), normals[i][2]/sqrt(normals[i][1]^2 + normals[i][2]^2)])
	#v1 = v - 2dot(v, N)*N
	@show push!(new_velocities, v - [2*(v[1]*normalized_normals[i][1] + v[2]*normalized_normals[i][2])*normalized_normals[i][1], 2*(v[1]*normalized_normals[i][1] + v[2]*normalized_normals[i][2])*normalized_normals[i][2]])
	
	if purities[i] == 1
		for j = 1:length(new_velocities)
			rects1, disks1, times1, purities1 = collision_purity_flat(collision_rectangles[j], new_velocities[j], c, 1)		
		end
	end
end


for i = 1:length(rects1)
    if disks1[i] == 1
        draw_rectangle(rects1[i][1].lo, rects1[i][2].lo, diam(rects1[i][1]), diam(rects1[i][2]), "red")
    elseif disks1[i] == 2 #&& purities[i] == 1
        draw_rectangle(rects1[i][1].lo, rects1[i][2].lo, diam(rects1[i][1]), diam(rects1[i][2]), "green")

    elseif disks1[i] == 3 #&& purities[i] == 1
        draw_rectangle(rects1[i][1].lo, rects1[i][2].lo, diam(rects1[i][1]), diam(rects1[i][2]), "blue")

    end
end
axis([-2, 7, -2, 7])




