using Base.Test
using KrawczykMethod

# Transcendental true roots were found using Wolfram Mathematica

a = Interval(-21, 24)

#enclosure(x) = Interval(x - 5000*eps(BigFloat), x + 5000*eps(BigFloat))
delta = 1e-7
enclosure(x) = Interval(x - delta, x + delta)

function checkroots(f, trueroots)
	roots = krawczyk(f, a, 64)
	# Check if both ends of the found interval lie within the 10eps-enclosure of the true root
	n = 0
	for i = 1:length(roots)
		if belong(roots[i][1].lo, enclosure(trueroots[i])) == false && belong(roots[i][1].hi, enclosure(trueroots[i])) == false
			n += 1
		end
	end
	n == 0
end	

f(x) = x^2 - 4
trueroots = [-2, 2]
@test checkroots(f, trueroots)

f(x) = 5040 - 3828x - 2356x^2 + 1009x^3 + 200x^4 - 62x^5 - 4x^6 + x^7
trueroots = [-6, -4, -2, 1, 3, 5, 7]
@test checkroots(f, trueroots)

f(x) = exp(x^2) - 4
a = Interval(-4, 4)
trueroots = [-sqrt(log(4)), sqrt(log(4))]
@test checkroots(f, trueroots)

# The next function is undefined in [-1, 1], so we split the interval manually
f(x) = log(x^2 - 1) + x + 1
a = Interval(-10, -1.01)
trueroots = [-3.274277505524962, -1.7894889253208424]
@test checkroots(f, trueroots)

f(x) = log(x^2 - 1) + x + 1
a = Interval(1.01, 10)
trueroots = [1.0617135926109396]
@test checkroots(f, trueroots)


