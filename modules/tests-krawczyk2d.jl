## Exercises from the Ramon Moore book

using Base.Test
using KrawczykMethod2D

# Transcendental true roots were found using Wolfram Mathematica

a = [Interval(-10, 10), Interval(-10, 10)]

enclosure(x) = Interval(x - 5000*eps(BigFloat), x + 5000*eps(BigFloat))

function checkroots(f, trueroots)
	roots = krawczyk2d(f, a)
	# Check if both ends of the found interval lie within the 10eps-enclosure of the true root
	n = 0
	for i = 1:length(roots)
		if belong(roots[i][1].lo, enclosure(trueroots[i])) == false && belong(roots[i][1].hi, enclosure(trueroots[i])) == false
			n += 1
		end
	end
	n == 0
end	

f(x) = [x[1]^2 + x[2]^2 - 1, x[1] - x[2]^2]
a = [Interval(0.5, 0.8), Interval(0.6, 0.9)]
krawczyk2d(f, a, 64)
trueroots = [(sqrt(5) - 1)/2, sqrt((sqrt(5) - 1)/2)]
