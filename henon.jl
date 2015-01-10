
using KrawczykMethod2D

# Henon function
f(x, a, b) = [1 - a*x[1]^2 + x[2], b*x[1]]

compose(f::Function, g::Function) = (x, a, b) -> f(g(x, a, b), a, b)

function T(x, a, b, n)
	F = f
	for i = 1:n-1
		F = compose(F, f)
	end
	return F(x, a, b)
end

h(x, a, b, n) = T(x, a, b, n) - x

function g(x, a, b, n)
	h1(n) = h(x, a, b, n)
	g1 = [one(Interval), one(Interval)]
	for i = 1:n-1
		g1 .*= h1(i)
	end
	g1
end

hx(x, a, b, n) = h(x, a, b, n)./g(x, a, b, n)

zeros(a, b, n, x1, x2, p1, p2, precision) = krawczyk2d(x -> h(x, a, b, n), [Interval(x1, x2), Interval(p1, p2)], precision)

#=
for i = 1:40
	x0 = [(2rand()-1)/2, (2rand()-1)/2]
	x1 = T(x0, 1.4, 0.3, 10^5)
	for j = 1:100
		output = x1
		println("$(output[1]) $(output[2])")
		x1 = f(output, 1.4, 0.3)
	end
end
=#


x0 = [3*(2rand()-1)/2, 3*(2rand()-1)/2]
for j = 1:1e6
	output = x0
	println("$(output[1]) $(output[2])")
	x0 = f(output, 1.4, 0.3)
end


#=
a = 1.4
b = 0.3
println("1-periodic: ")
println(mid(zeros(a, b, 1, -10, 10, -3, 3, 64)))

println("2-periodic: ")
println(mid(zeros(a, b, 2, -10, 10, -3, 3, 64)))

println("3-periodic: ")
println(mid(zeros(a, b, 3, -10, 10, -3, 3, 64)))

println("4-periodic: ")
println(mid(zeros(a, b, 4, -10, 10, -3, 3, 64)))
=#
