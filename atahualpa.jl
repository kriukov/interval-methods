# https://github.com/johnmyleswhite/ContinuedFractions.jl
using ContinuedFractions

conv(x) = convergents(ContinuedFraction(x), Rational)

function frac(alpha, epsilon)
	n = 1
	while abs(alpha*den(conv(alpha)[n]) - num(conv(alpha)[n])) >= epsilon
		n += 1
	end
	(den(conv(alpha)[n]), num(conv(alpha)[n]))
end

function eff(m, b, epsilon)
	kn = 0
	#i = 0
	
	while b > epsilon && 1 - b > epsilon
		if b < 0.5
			#println("way 1")
			(q, p) = frac(m, 2b)
		else
			#println("way 2")
			(q, p) = frac(m, 2*(1 - b))
		end
		#bb = m*q + b - ifloor(m*q + bb)
		b = mod(m*q + b, 1)
		kn += q
		#i += 1
		#if i >= 3000
			#println("cycle break: too many iterations")
			#break
		#end
	end
	#println("cycle ended")
	q = kn
	p = ifloor(m*q) + 1
	return (q, p)	
	#return i
end


#= piece of code to match m > 1 to analogous m < 1
		
	=#


function first_collision(x, y, vx, vy, delta)
	
	# Normalize velocity if it wasn't normalized
	v = sqrt(vx^2 + vy^2)
	vx1 = vx/v
	vy1 = vy/v
	vx = vx1
	vy = vy1
			
	m = vy/vx
	b = y - m*x

	#=	
	if x != 0	
		if m < 1
				b1 = m + b # m*1 + b
			if b1 < delta
				return (1, 0)
			elseif 1 - b1 < delta
				return (1, 1)
			else
				b = b1
				x = 1
				y = b1
			end
		elseif m > 1
				mp = 1/m
				bpp = 1 - b/m
				b1 = mp + bpp
				if b1 < delta
					return (0, 1)
				elseif 1 - b1 < delta
					return (1, 1)
				else
					b = m*(1 - bpp)
					x = b1 - 1
					y = 1
				end	
		end
	end
	=#

	


	
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

#c = sqrt(1.001) # Irrationalizing coefficient
# Split the 1st quadrant into N = 30 equal angular parts
# For each n < N, slope is tg(pi n/(2N))

#N = 30

#=
for n = 1:N-1
	loops = 0
	for m = 1:N
		for b = 1:N
			x = eff(tan(pi*n/(2N))*c, b/N*c, n/N)
			if x == 3000
				loops += 1
				break
			end
		end
	end
	println(n/N, " ", loops/N^2)
end
=#

#=
for n = 1:N-1
	loops = 0
	@show n
	for m = 1:N
		for b = 1:N
			m = tan(pi*rand()/2)
			x = eff(m, rand(), n/N*sqrt(1 + m^2))
			if x == 300
				loops += 1
				break
			end
		end
	end
	println(n/N, " ", loops/N^2)
end
=#

# Testing
#=
r = 0.1
x = 0
y = 0.445
for deg = 1:89
	if deg != 45
		vx = cos(deg*pi/180)
		vy = sin(deg*pi/180)
		println(deg, " ", first_collision(x, y, vx, vy, r/vx))
	end
end
=#

# Measuring the times

#= This one still gives error due to ContinuedFractions limits. Trying to work around
function measuring()
for i = 1:100
	r = i/1000
	for j = 1:100
		for k = 1:100
			vx = j*sqrt(0.999)
			vy = k*sqrt(1.002)
			#println(i, " ", j, " ", k, " ", vy/vx)
			first_collision(0, 0.445, vx, vy, r/vx)
		end
	end
end
end

@time measuring()
=#

x = 0
y = 0.445

#=
function measuring()
for i = 1:100
	r = i/1000
	for deg = 1:89
		if deg != 45
			vx = cos(deg*pi/180)
			vy = sin(deg*pi/180)
			#println(r, " ", deg, " ", first_collision(x, y, vx, vy, r/vx))
			first_collision(x, y, vx, vy, r/vx)
		end
	end
end
end

@time measuring()
=#

# with println: 112.495258954 s
# without println: 109.723781459 s



function measuring2()
	r = 0.000001
	for deg = 46:89
		if deg != 45
			vx = cos(deg*pi/180)
			vy = sin(deg*pi/180)
			println(r, " ", deg, " ", first_collision(x, y, vx, vy, r/vx))
			first_collision(x, y, vx, vy, r/vx)
		end
	end
end

@time measuring2()




