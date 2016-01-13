using IntervalArithmetic
using AutoDiff
using PurityIntervals
using KrawczykMethod2D

using FastAnonymous

# The vector x = [ω, θ] in Birkhoff mapping

function alpha(n, m)
	if n > 0 && n <= 3 && m > 0 && m <= 3 && typeof(n) == Int && typeof(m) == Int && n != m
		k = pi/3
		if n == 1 && m == 2
			return 0
		elseif n == 1 && m == 3
	 		return k
		elseif n == 2 && m == 3
	 		return 2k
		elseif n == 2 && m == 1
	 		return 3k
		elseif n == 3 && m == 1
	 		return 4k
		elseif n == 3 && m == 2
	 		return 5k
	 	end
	else
	    error("Invalid parameters: should be 1, 2, 3 and not equal")
	end
end

function T(x, n, m, r)
	ω, θ = x
    #@show x, n, m 
	ω_next = ω - r*(ω*cos(θ - alpha(n, m)) + √(1 - ω^2)*sin(θ - alpha(n, m)))
	θ_next = mod(θ + big(pi) + asin(ω) + asin(ω_next), 2π)
	#θ_next = θ + big(pi) + asin(ω) + asin(ω_next)
	
    [ω_next, θ_next]
end

function T(x::Array{MultiDimInterval, 1}, n, m, r)
	sol = MultiDimInterval[]
	for i = 1:length(x)
		push!(sol, T(x[i], n, m, r))
	end
	return sol
end


# After the modification above, it still shows the same error

#=
function T(x::Array{PurityInterval, 1}, n, m, r)
	sol = MultiDimInterval[]
	for i = 1:length(x)
		xi = x[i]
		push!(sol, PurityInterval(T(xi, n, m, r), purity(xi -> T(xi, n, m, r), xi)))
	end
	return sol
end
=#

using PyPlot

distance_between_centers = 6.

# x = [ω, θ] and n is array of path points, e.g. [1, 2, 3, 1, 3]
function path(x, n)
    for i = 1:length(n)-1
        if n[i] != n[i+1]
            x = T(x, n[i], n[i+1], distance_between_centers)
        else
            error("Cannot hit the same disk twice in succession")
        end
    end
    x
end

# Non-bisection inefficient plotting
function plot_band(f, N)

	# Split the phase space w = [-1, 1], th = [0, pi/3] into 30x30 rectangles
	rw = Interval(-0.95, 0.95)
	rth = Interval(0.05, float(pi)/3-0.05)
	deltaw = diam(rw/N)
	deltath = diam(rth/N)

	points = Array{BigFloat, 1}[]
	purities = Int[]

	for i = 1:N
		for j = 1:N
			rect = [rw.lo + i*deltaw - Interval(0, deltaw), rth.lo + j*deltath - Interval(0, deltath)]
			p = purity(f, rect)
			#println([rect, p])
			push!(points, mid(rect))
			push!(purities, p)

		end
	end

	x = [Float64(xx[1]) for xx in points]
	y = [Float64(xx[2]) for xx in points]

	p = purities

	for i in (0, 1, -1)
		plot(y[p.==i], x[p.==i], "o", label="$i", alpha = 0.5)
	end
	legend()

end

#plot_band(f12)
#plot_band(f13)

#plot_band(twice_123, 100)
#plot_band(x -> path(x, [1, 2, 3]), 100)

rect = [Interval(-0.999, 0.999), Interval(0.001, big(pi)/3-0.001)]
points = Array{Interval, 1}[]
purities = Int[]

# More efficient plotting without dirty rectangles
function plot_band_bisection(f, rect, tol)
    limitrect(rect) = max(diam(rect)[1], diam(rect)[2]) < tol
    p = purity(f, rect)
    if p != -1
        if p == 1
            #println("Went to p=1")
            push!(points, rect)
            push!(purities, p)
            #println("Clean recorded")
        elseif p == 0
            #println("Went to p=0")
            if limitrect(rect)
                push!(points, rect)
                push!(purities, p)
                #println("Unclean recorded")
            else
                pieces = bisect(rect)
                for i = 1:4
                    plot_band_bisection(f, pieces[i], tol)
                end
            end
        end
    end
    points, purities
end

# Plotting that includes dirty rectangles
function plot_band_bisection_dirty(f, rect, tol)
    limitrect(rect) = max(diam(rect)[1], diam(rect)[2]) < tol
    p = purity(f, rect)
    if p == 1 || p == -1
        push!(points, rect)
        push!(purities, p)
    elseif p == 0
        if limitrect(rect)
            push!(points, rect)
            push!(purities, p)
        else
            pieces = bisect(rect)
            for i = 1:4
                plot_band_bisection_dirty(f, pieces[i], tol)
            end
        end
    end
    points, purities
end

# Using midpoints of rectangles! It works
#=
pts = mid(plot_band_bisection(x -> path(x, [1, 2]), rect, 1e-3))
x = [Float64(xx[1]) for xx in pts]
y = [Float64(xx[2]) for xx in pts]
plot(y, x, "o")
=#

# Function to draw rectangles
using PyCall
@pyimport matplotlib.patches as patches
rectangle = patches.Rectangle
function draw_rectangle(x, y, xwidth, ywidth, color="grey")
    ax = gca()
    ax[:add_patch](rectangle((x, y), xwidth, ywidth, facecolor=color, alpha=0.5, linewidth=0))
end


#=
plot_band_bisection(x -> path(x, [1, 2]), rect, 1e-3)
println(points)
println(purities)
=#

#result = plot_band_bisection(x -> path(x, [1, 2]), rect, 1e-2)

# Printing output for WM
#=
for i = 1:length(result[1])
    println("$(result[1][i][1].lo) $(result[1][i][1].hi) $(result[1][i][2].lo) $(result[1][i][2].hi) $(result[2][i])")
end
=#



#= Draw domain of T[...] without dirty

points, purities = plot_band_bisection(x -> path(x, [1, 2]), rect, 1e-2)
for i = 1:length(points)
    if purities[i] == 0
        draw_rectangle(points[i][2].lo, points[i][1].lo, diam(points[i][2]), diam(points[i][1]), "green")
    elseif purities[i] == 1
        draw_rectangle(points[i][2].lo, points[i][1].lo, diam(points[i][2]), diam(points[i][1]), "blue")
    end
end
axis([0, pi/3, -1, 1])
savefig("T12-domain.pdf")
=#

# Draw domain of T12 with dirty
#=
points, purities = plot_band_bisection_dirty(x -> path(x, [1, 2]), rect, 1e-3)
for i = 1:length(points)
    if purities[i] == 0
        draw_rectangle(points[i][2].lo, points[i][1].lo, diam(points[i][2]), diam(points[i][1]), "green")
    elseif purities[i] == 1
        draw_rectangle(points[i][2].lo, points[i][1].lo, diam(points[i][2]), diam(points[i][1]), "blue")
    elseif purities[i] == -1
        draw_rectangle(points[i][2].lo, points[i][1].lo, diam(points[i][2]), diam(points[i][1]), "red")
    end
end
axis([0, pi/3, -1, 1])
=#


#f(x) = path(x, [1, 2, 3, 1]) - x
#a = [Interval(-0.5-0.001, -0.5+0.001), Interval(pi/6-0.001, pi/6+0.001)]

#f(x) = path(x, [1, 2, 3, 1, 3, 1]) - x

#krawczyk2d_purity(f, rect, 64, 1e-4)
#krawczyk2d_purity(f, a, 64, 1e-4)

#rect = [Interval(-0.95, 0.95), Interval(0.05, big(pi)/3-0.05)]
#rect = [Interval(-0.51, -0.02), Interval(0.05, big(pi)/6-0.05)]
#krawczyk2d_purity(f, rect, 64, 1e-4)

# Draw the image/range
#points, purities = plot_band_bisection(x -> path(x, [1, 2, 3, 1]), rect, 5e-5)
#=
points0 = readdlm("points.dat")
points = MultiDimInterval[]
for i = 1:length(points0)
    push!(points, eval(parse(points0[i])))
end

purities = readdlm("purities.dat")

for i = 1:length(points)
    if purities[i] == 1
        @show point_image = path(points[i], [1, 2, 3, 1])
        if typeof(point_image[1]) == Interval && typeof(point_image[2]) == Interval
            draw_rectangle(point_image[2].lo, point_image[1].lo, diam(point_image[2]), diam(point_image[1]), "blue")
        elseif typeof(point_image[1]) == Interval && typeof(point_image[2]) == IntUnion
            for i = 1:length(point_image[2].union)
                draw_rectangle(point_image[2].union[i].lo, point_image[1].lo, diam(point_image[2].union[i]), diam(point_image[1]), "blue")
            end
        elseif typeof(point_image[1]) == IntUnion && typeof(point_image[2]) == Interval
            for i = 1:length(point_image[1].union)
                draw_rectangle(point_image[2].lo, point_image[1].union[i].lo, diam(point_image[2]), diam(point_image[1].union[i]), "blue")
            end
        elseif typeof(point_image[1]) == IntUnion && typeof(point_image[2]) == IntUnion
            for i = 1:length(point_image[1].union)
                for j = 1:length(point_image[2].union)
                    draw_rectangle(point_image[2].union[j].lo, point_image[1].union[i].lo, diam(point_image[2].union[j]), diam(point_image[1].union[i]), "blue")
                end
            end
        end
    end
end
axis([0, 2pi, -1, 1])
=#

#savefig("T1231-range.pdf")
#=
f(x) = path(x, [1, 2, 3]) - x
points, purities = plot_band_bisection(f, rect, 2e-3)
points1, purities1 = plot_band_bisection(x->path(x,[1,2])-x, rect, 1e-3)
for i = 1:length(points)
    if purities[i] == 0
        draw_rectangle(points[i][2].lo, points[i][1].lo, diam(points[i][2]), diam(points[i][1]), "green")
    elseif purities[i] == 1
        draw_rectangle(points[i][2].lo, points[i][1].lo, diam(points[i][2]), diam(points[i][1]), "blue")
    end
end
for i = 1:length(points1)
    if purities1[i] == 0
        draw_rectangle(points1[i][2].lo, points1[i][1].lo, diam(points1[i][2]), diam(points1[i][1]), "black")
    elseif purities1[i] == 1
        draw_rectangle(points1[i][2].lo, points1[i][1].lo, diam(points1[i][2]), diam(points1[i][1]), "red")
    end
end
axis([0, pi/3, -1, 1])
=#


c = MultiDimInterval[]
r0 = 6
push!(c, [Interval(0), Interval(0)], [Interval(r0), Interval(0)], [Interval(r0/2), Interval(r0*sqrt(3)/2)])




function T0(x, c, n, m)
	ω, θ = x
	
	r = Any[]
    a = Any[]
    for i = 1:length(c)
        for j = 1:length(c)
            push!(r, [sqrt((c[i][1] - c[j][1])^2 + (c[i][2] - c[j][2])^2), i, j])
            push!(a, [atan((c[j] - c[i])[2]/(c[j] - c[i])[1]), i, j])
        end
    end
	
	r_nm = 0; a_nm = 0
	
	for i = 1:length(r)
	    if r[i][2] == n && r[i][3] == m
	        r_nm = r[i][1]
	    end
	end
	
	for i = 1:length(a)
        if a[i][2] == n && a[i][3] == m
            a_nm = a[i][1]
        end
	end
	
	
	ω_next = ω - r_nm*(ω*cos(θ - a_nm) + √(1 - ω^2)*sin(θ - a_nm))
	θ_next = mod(θ + big(pi) + asin(ω) + asin(ω_next), 2π)
	
	#IntUnion2D(ω_next, θ_next)
    [ω_next, θ_next]
end

function path_general(x, c, n)
    for i = 1:length(n)-1
        if n[i] != n[i+1]
            x = T0(x, c, n[i], n[i+1])
        else
            error("Cannot hit the same disk twice in succession")
        end
    end
    x
end

function draw_phase_space(c, n, rect, tol)
    points, purities = plot_band_bisection(x -> path_general(x, c, n), rect, tol)
    for i = 1:length(points)
        if purities[i] == 0
            draw_rectangle(points[i][2].lo, points[i][1].lo, diam(points[i][2]), diam(points[i][1]), "green")
        elseif purities[i] == 1
            draw_rectangle(points[i][2].lo, points[i][1].lo, diam(points[i][2]), diam(points[i][1]), "blue")
        end
    end
    axis([0, pi/3, -1, 1])
end

function find_periodic_orbits(c, n, rect, prec, tol)
   f(x) = path_general(x, c, n)
   #f = @anon x -> path_general(x, c, n) - x
   krawczyk2d_purity_periodic(f, rect, prec, tol)
end

# Time/distance function
function distance(x, c, n, m)
	ω, θ = x
	
	r = Any[]
    a = Any[]
    for i = 1:length(c)
        for j = 1:length(c)
            push!(r, [sqrt((c[i][1] - c[j][1])^2 + (c[i][2] - c[j][2])^2), i, j])
            push!(a, [atan((c[j] - c[i])[2]/(c[j] - c[i])[1]), i, j])
        end
    end
	
	r_nm = 0; a_nm = 0
	
	for i = 1:length(r)
	    if r[i][2] == n && r[i][3] == m
	        r_nm = r[i][1]
	    end
	end
	
	for i = 1:length(a)
        if a[i][2] == n && a[i][3] == m
            a_nm = a[i][1]
        end
	end
	
	
	ω_next = ω - r_nm*(ω*cos(θ - a_nm) + √(1 - ω^2)*sin(θ - a_nm))
	θ_next = mod(θ + big(pi) + asin(ω) + asin(ω_next), 2π)
	
	t_next = sqrt(r_nm^2 + 2r_nm*(cos(θ_next - a_nm) - cos(θ - a_nm)) + 2 - 2*cos(θ_next - θ))
    t_next
end

# Draw phase space from n-th disk to all other - here "n" is a number

function draw_phase_space_general(c, n, rect, tol)
    points = Array{MultiDimInterval, 1}[]
    purities = Array{Int, 1}[]
    for i = 1:length(c)
        if i != n
            points_i, purities_i = plot_band_bisection(x -> path_general(x, c, [n, i]) - x, rect, tol)
            push!(points, points_i)
            push!(purities, purities_i)
        end
    end
    
    for i = 1:length(points)
        for j = 1:length(points[i])
            if purities[i][j] == 0
                draw_rectangle(points[i][j][2].lo, points[i][j][1].lo, diam(points[i][j][2]), diam(points[i][j][1]), "green")
            elseif purities[i][j] == 1
                draw_rectangle(points[i][j][2].lo, points[i][j][1].lo, diam(points[i][j][2]), diam(points[i][j][1]), "blue")
            end
        end
    end

    axis([0, pi/3, -1, 1])
    
end

# Birkhoff 1-2 and 1-3 belts in a symmetric 3-disk scatterer
#draw_phase_space_general(c, 1, rect, 1e-2)

# Another system of disks: #2 is in the way of #3 from #1

c1 = MultiDimInterval[]
push!(c1, [Interval(0), Interval(0)], [Interval(r0), Interval(0)], [Interval(r0/2), Interval(1.3)])

#draw_phase_space_general(c1, 1, rect, 1e-2)


#@time find_periodic_orbits(c, [1, 2], rect, 64, 1e-4)

#@time find_periodic_orbits(c, [1, 3, 1], [Interval(-0.012, 0.001), Interval(pi/3-0.01, pi/3+0.0011)], 64, 1e-4)

# Find periodic orbit 1-2-3-1 in isosceles system c1
find_periodic_orbits(c1, [1, 2, 3, 1], rect, 64, 1e-4)





