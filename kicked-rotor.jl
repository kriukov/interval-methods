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



# Function to output zeros of n-th order given the K and the limits (x1, x2) and (p1, p2)
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
g1(x) = h1(x, 0.7, 4)

#zeros1(x1, x2, p1, p2, precision) = krawczyk2d(g1, [Interval(x1, x2), Interval(p1, p2)], precision)
zeros(k, n, x1, x2, p1, p2, precision) = krawczyk2d(x -> h1(x, k, n), [Interval(x1, x2), Interval(p1, p2)], precision)

