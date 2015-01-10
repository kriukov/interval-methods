## Kicked rotor: periodic orbits calculation
# x is a vector (x[1], x[2]) = (x, p)

using KrawczykMethod2D

# Base function for the kicked rotor, x[1] is x, x[2] is p
f0(x, k) = mod([x[1] + k*sin(x[1]) + x[2], k*sin(x[1]) + x[2]], 2pi)
f(x, k) = [x[1] + k*sin(x[1]) + x[2], k*sin(x[1]) + x[2]]

# Function composition while keeping k arbitrary
compose(f::Function, g::Function) = (x, k) -> f(g(x, k), k)

# T^n (x) function - T(T(T(...T(x)))) n times
function T(x, k, n)
	F = f
	for i = 1:n-1
		F = compose(F, f)
	end
	return F(x, k)
end


# T^n (x) function - T(T(T(...T(x)))) n times
function T0(x, k, n)
	F = f0
	for i = 1:n-1
		F = compose(F, f0)
	end
	return F(x, k)
end


# n-periodic points equation
h(x, k, n) = T(x, k, n) - x

# Equation to find ONLY n-periodic points

function g(x, k, n)
	g1 = 1
	for i = 1:n-1
		g1 = g1*h(x, k, i)
	end
	h(x, k, n)/g1
end


# Functions to output zeros of n-th order given the K and the limits (x1, x2) and (p1, p2)
#zeros(K, n, x1, x2, p1, p2, precision) = krawczyk2d_mod(x -> h(x, K, n), [Interval(x1, x2), Interval(p1, p2)], precision)

top(x, k, n) = ceil(hi(T(x, k, n))/2pi)*2pi

function zeros(k, n, x1, x2, p1, p2, precision)
	array_sol = MultiDimInterval[]
	a = [Interval(x1, x2), Interval(p1, p2)]
	for p = 0:top(a, k, n)[1]/2pi
		for q = 0:top(a, k, n)[2]/2pi
			array_sol = vcat(array_sol, krawczyk2d(x -> h(x, k, n) - [2pi*p, 2pi*q], a, precision))
			#push!(array_sol, krawczyk2d(x -> h(x, k, n) - [p, q], a, precision))
		end
	end
	array_sol_nodupes = unique(array_sol)
	array_sol_nodupes
end

# Print out 2-periodic points for k = 1
#println(mid(zeros(1, 3, 0, 2pi, 0, 2pi, 64)))


# Function to output the plot points

#= turned off
for n = 1:40
	x0 = [2pi*rand(), 2pi*rand()]
	k = 1
	for i = 1:10000
		output = T0(x0, k, i)
		println("$(output[1]) $(output[2])")
	end
end
=#


function zeros1(k, x1, x2, p1, p2, precision)
	array_sol = MultiDimInterval[]
	a = [Interval(x1, x2), Interval(p1, p2)]
	for p = 0:top(a, k, 1)[1]/2pi
		for q = 0:top(a, k, 1)[2]/2pi
			array_sol = vcat(array_sol, krawczyk2d(x -> h(x, k, 1) - [2pi*p, 2pi*q], a, precision))
		end
	end
	array_sol_nodupes = unique(array_sol)
	array_sol_nodupes
end


function zeros2(k, x1, x2, p1, p2, precision)
	array_sol = MultiDimInterval[]
	a = [Interval(x1, x2), Interval(p1, p2)]
	for p = 0:top(a, k, 2)[1]/2pi
		for q = 0:top(a, k, 2)[2]/2pi
			func(x) = x -> (h(x, k, 2) - [2pi*p, 2pi*q])/(h(x, k, 1) - [2pi*p, 2pi*q])
			array_sol = vcat(array_sol, krawczyk2d(func, a, precision))
		end
	end
	array_sol_nodupes = unique(array_sol)
	array_sol_nodupes
end

function zeros3(k, x1, x2, p1, p2, precision)
	array_sol = MultiDimInterval[]
	a = [Interval(x1, x2), Interval(p1, p2)]
	for p = 0:top(a, k, 3)[1]/2pi
		for q = 0:top(a, k, 3)[2]/2pi
			array_sol = vcat(array_sol, krawczyk2d(x -> (h(x, k, 3) - [2pi*p, 2pi*q])/((h(x, k, 2) - [2pi*p, 2pi*q])*(h(x, k, 1) - [2pi*p, 2pi*q])), a, precision))
		end
	end
	array_sol_nodupes = unique(array_sol)
	array_sol_nodupes
end

function zeros4(k, x1, x2, p1, p2, precision)
	array_sol = MultiDimInterval[]
	a = [Interval(x1, x2), Interval(p1, p2)]
	for p = 0:top(a, k, 4)[1]/2pi
		for q = 0:top(a, k, 4)[2]/2pi
			array_sol = vcat(array_sol, krawczyk2d(x -> (h(x, k, 4) - [2pi*p, 2pi*q])/((h(x, k, 3) - [2pi*p, 2pi*q])*(h(x, k, 2) - [2pi*p, 2pi*q])*(h(x, k, 1) - [2pi*p, 2pi*q])), a, precision))
		end
	end
	array_sol_nodupes = unique(array_sol)
	array_sol_nodupes
end

#println(mid(zeros2(1, 0, 2pi, 0, 2pi, 64)))

# Print out periodic points for k = 1

#=
println("1-periodic: ")
println(mid(zeros(1, 1, 0, 2pi, 0, 2pi, 64)))

println("2-periodic: ")
println(mid(zeros(1, 2, 0, 2pi, 0, 2pi, 64)))

println("3-periodic: ")
println(mid(zeros(1, 3, 0, 2pi, 0, 2pi, 64)))

println("4-periodic: ")
println(mid(zeros(1, 4, 0, 2pi, 0, 2pi, 64)))
=#
