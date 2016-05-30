## Exercises from the Ramon Moore book

using Base.Test
using KrawczykMethod2DVN

# Transcendental true roots were found using Wolfram Mathematica

enclosure(x) = Interval(x - 1e-9, x + 1e-9)

function checkroots(f, trueroots)
	roots = krawczyk2d(f, a)
	# Check if both ends of each found interval lie within the enclosure of the true root
	n = 0
	for i = 1:length(roots)
		for j = 1:2 # For 2D
			if in(roots[i][j].lo, enclosure(trueroots[i][j])) == false && in(roots[i][j].hi, enclosure(trueroots[i][j])) == false
				n += 1
			end
		end
	end
	n == 0
end	

f(x) = [x[1]^2 + x[2]^2 - 1, x[1] - x[2]^2]
a = IntervalBox(Interval(0, 1), Interval(0, 1))
trueroots = Array[[(sqrt(5) - 1)/2, sqrt((sqrt(5) - 1)/2)]]
@test checkroots(f, trueroots)

f(x) = [x[1]^2 - 2, x[2]^2 - 2]
a = IntervalBox(Interval(0, 5), Interval(0, 5))
trueroots = Array[[sqrt(2), sqrt(2)]]
@test checkroots(f, trueroots)

f(x) = [20 - 20x[1] - x[2], x[1] - x[2]/20 - 1e-9*exp(x[2]/0.052)]
a = IntervalBox(Interval(0, 2), Interval(0, 2))
trueroots = Array[[0.9464142468335176, 1.0717150633296477]]
@test checkroots(f, trueroots)

f(x) = [0.5*(-(17.76x[1] - 103.79x[1]^2 + 229.62x[1]^3 - 226.31x[1]^4 + 83.72x[1]^5) + x[2]), 0.2*(-x[1] - 1.5x[2] + 1.2)]
a = IntervalBox(Interval(0, 1), Interval(0, 1))
trueroots = Array[[0.06263595920972119, 0.758242693860186], [0.2853687241300338, 0.6097541839133109], [0.8844295888702942, 0.21038027408647064]]
@test checkroots(f, trueroots)


