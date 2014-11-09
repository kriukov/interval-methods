# https://github.com/johnmyleswhite/ContinuedFractions.jl
using ContinuedFractions

conv(x) = convergents(ContinuedFraction(x), Rational)

function frac(alpha, epsilon)
	n = 1
	while abs(alpha*den(conv(alpha)[n]) - num(conv(alpha)[n])) > epsilon
		n += 1
	end
	(den(conv(alpha)[n]), num(conv(alpha)[n]))
end

function eff(m, b, epsilon)
	bb = b
	kn = 0
	while bb > epsilon && 1 - bb > epsilon
		if bb < 0.5
			println("way 1")
			@show (q, p) = frac(m, 2bb)
		else
			println("way 2")
			@show (q, p) = frac(m, 2*(1 - bb))
		end
		bb = m*q + b - ifloor(m*q + b)
		kn = kn + q		
	end
	println("cycle ended")
	q = kn
	p = ifloor(m*q) + 1
	return (q, p)	
end

function first_collision(x, y, vx, vy, delta)
	m = vy/vx
	b = y - m*x
	if vx > 0 && vy > 0
		(q, p) = eff(m, b, delta)
		p = ifloor(m*q) + 1
	elseif vx < 0 && vy > 0
		m = -m
		(q, p) = eff(m, b, delta)
		p = ifloor(m*q) + 1
		m = -m
		q = -q
	elseif vx < 0 && vy < 0
		b = 1 - b
		(q, p) = eff(m, b, delta)
		b = 1 - b
		p = -(ifloor(m*q) + 1)
		q = -q
	elseif vx > 0 && vy < 0
		b = 1 - b
		m = -m
		(q, p) = eff(m, b, delta)
		b = 1 - b
		m = -m
		p = -(ifloor(m*q) + 1)
	end
	return (q, p)
end

