
# Cat map function
cat(x) = mod([2 1; 1 1]*x, 1)

# Function composition while keeping K arbitrary
compose(f::Function, g::Function) = x -> f(g(x))

# T^n (x) function - T(T(T(...T(x)))) n times
function T(x, n)
	F = cat
	for i = 1:n-1
		F = compose(F, cat)
	end
	return F(x)
end

# Function to output the plot points at 40 random init. cond.
#= off
for n = 1:40
	x0 = [rand(), rand()]
	for i = 1:500
		output = T(x0, i)
		println("$(output[1]) $(output[2])")
	end
end
=#

# Doing the same for one initial condition
x0 = [0.4, 0.7]
for i = 1:500
	output = T(x0, i)
	println("$(output[1]) $(output[2])")
end

