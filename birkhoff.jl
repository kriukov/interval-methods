using IntervalArithmetic
using AutoDiff
using PurityIntervals
using KrawczykMethod2D

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
    ax[:add_patch](rectangle((x, y), xwidth, ywidth, facecolor=color, alpha=0.5))
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


#savefig("T1231-range.pdf")

