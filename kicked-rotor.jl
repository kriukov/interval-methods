## Kicked rotor: periodic orbits calculation
# x is a vector (x[1], x[2]) = (x, p)

using KrawczykMethod2D

# Function to output the plot points

#= turned off
for n = 1:40
	x0 = [2pi*rand(), 2pi*rand()]
	K = 1.5
	for i = 1:10000
		output = T(x0, K, i)
		println("$(output[1]) $(output[2])")
	end
end
=# 



# Base function for the kicked rotor, x[1] is x, x[2] is p
f(x, k) = mod([x[1] + k*sin(x[1]) + x[2], k*sin(x[1]) + x[2]], 2pi)
#f(x, k) = [x[1] + K*sin(x[1]) + x[2], K*sin(x[1]) + x[2]]

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

# n-periodic points equation
h(x, k, n) = T(x, k, n) - x



# Functions to output zeros of n-th order given the K and the limits (x1, x2) and (p1, p2)
#zeros(K, n, x1, x2, p1, p2, precision) = krawczyk2d_mod(x -> h(x, K, n), [Interval(x1, x2), Interval(p1, p2)], precision)


f1(x, k) = mod21([x[1] + k*sin(x[1]) + x[2], k*sin(x[1]) + x[2]], 2pi)
f2(x, k) = mod22([x[1] + k*sin(x[1]) + x[2], k*sin(x[1]) + x[2]], 2pi)
f3(x, k) = mod23([x[1] + k*sin(x[1]) + x[2], k*sin(x[1]) + x[2]], 2pi)
f4(x, k) = mod24([x[1] + k*sin(x[1]) + x[2], k*sin(x[1]) + x[2]], 2pi)

function T1(x, k, n)
	F = f1
	for i = 1:n-1
		F = compose(F, f1)
	end
	return F(x, k)
end
h1(x, k, n) = T1(x, k, n) - x

function T2(x, k, n)
	F = f2
	for i = 1:n-1
		F = compose(F, f2)
	end
	return F(x, k)
end
h2(x, k, n) = T2(x, k, n) - x

function T3(x, k, n)
	F = f3
	for i = 1:n-1
		F = compose(F, f3)
	end
	return F(x, k)
end
h3(x, k, n) = T3(x, k, n) - x

function T4(x, k, n)
	F = f4
	for i = 1:n-1
		F = compose(F, f4)
	end
	return F(x, k)
end
h4(x, k, n) = T4(x, k, n) - x

#zeros1(x1, x2, p1, p2, precision) = krawczyk2d(g1, [Interval(x1, x2), Interval(p1, p2)], precision)
zeros1(k, n, x1, x2, p1, p2, precision) = krawczyk2d(x -> h1(x, k, n), [Interval(x1, x2), Interval(p1, p2)], precision)
zeros2(k, n, x1, x2, p1, p2, precision) = krawczyk2d(x -> h2(x, k, n), [Interval(x1, x2), Interval(p1, p2)], precision)
zeros3(k, n, x1, x2, p1, p2, precision) = krawczyk2d(x -> h3(x, k, n), [Interval(x1, x2), Interval(p1, p2)], precision)
zeros4(k, n, x1, x2, p1, p2, precision) = krawczyk2d(x -> h4(x, k, n), [Interval(x1, x2), Interval(p1, p2)], precision)

println(zeros1(0.8, 2, 0, 2pi, 0, 2pi, 64))
