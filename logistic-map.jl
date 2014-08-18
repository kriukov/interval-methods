using NewtonMethod

# Function composition
compose(f::Function, g::Function) = x -> f(g(x))

# Logistic function T(x) of the first order
T1(x) = Interval(r)*(Interval(1/4)-(x-Interval(1/2))^2)

# T^p (x) function - T(T(T(...T(x)))) p times
function T(x, p)
	F = T1
	for i = 1:p-1
		F = compose(F, T1)
	end
	return F(x)
end

# Equation to solve for periodic points

f(x, p) = T(x, p) - x

point(f) = newton(f, Interval(0.0001, 1.), 50)

r = 1.5
f(x) = T(x, 1) - x

println("$r $(point(f))")

r = 2
f(x) = T(x, 1) - x
println("$r $(point(f))")

r = 2.5
f(x) = T(x, 1) - x
println("$r $(point(f))")

r = 3
f(x) = T(x, 1) - x
println("$r $(point(f))")

r = 3.2
f(x) = (T(x, 2) - x)/(T(x, 1) - x)
println("$r $(point(f))")

r = 3.5
f(x) = (T(x, 4) - x)/(T(x, 2) - x)
println("$r $(point(f))")

r = 3.75
for i = 1:8
	f(x) = T(x, i) - x
	println("$r $(point(f))")
end
