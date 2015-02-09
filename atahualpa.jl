function frac(alpha, epsilon)
	p = 0; q = 0

	mm = zeros(10000)
	a = zeros(Int, 10000)
	h = zeros(Int, 10000)
	k = zeros(Int, 10000)
	
	h[1] = 0
	h[2] = 1
	k[1] = 1
	k[2] = 0
	mm[1] = alpha
	a[1] = ifloor(mm[1])

	i = 1
	while abs(k[i+1]*mm[1]-h[i+1]) > epsilon
		mm[i+1] = 1/(mm[i] - a[i])
		a[i+1] = ifloor(mm[i+1])
		h[i+2] = a[i]*h[i+1] + h[i]
		k[i+2] = a[i]*k[i+1] + k[i]
		p = h[i+2]
		q = k[i+2]
		i += 1
	end

	q,p
end 

function eff(m, b, epsilon)
	kn = 0
	#i = 0
	
	while b > epsilon && 1 - b > epsilon
		if b < 0.5
			(q, p) = frac(m, 2b)
		else
			(q, p) = frac(m, 2*(1 - b))
		end
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
		p = -ifloor(m*q) #p = -(ifloor(m*q) + 1)
		q = -q
	elseif vx > 0 && vy < 0
		b = 1 - b
		m = -m
		(q, p) = eff(m, b, delta)
		b = 1 - b
		#m = -m # this line should not be here
		p = -ifloor(m*q) #p = -(ifloor(m*q) + 1)
	end
	return q, p
end

# Gives the new velocity and the new coordinates (which are a point of collision) having the initial position/velocity, the obstacle number and radius
#= This function gives more error than the classic function. Classic function is below.
function collide(q, p, x, y, vx, vy, r)
	phi = acos(vx)
	k = vy/vx
	b = y - k*x
	beta = phi - asin((k*q + b - p)*vx/r)
	v1 = [vx, vy] - 2*[cos(beta), sin(beta)]*cos(phi - beta)
	x1 = q - r*cos(beta)
	y1 = p - r*sin(beta)
	return x1, y1, v1[1], v1[2]
end
=#

# This one is based on the classic algorithm
function collide(q, p, x, y, vx, vy, r)
	r0 = [x, y]
	v0 = [vx, vy]
	R = [q, p]
	crossz(x, y) = x[1]*y[2] - x[2]*y[1]
	discr = norm(v0)^2*r^2 - (crossz(v0, r0-R))^2
	t1 = (-dot(v0, r0-R) - sqrt(discr))/norm(v0)^2
	N0 = r0 + v0*t1 - R
	N = N0/norm(N0)			
	v1 = v0 - 2*dot(v0, N)*N
	r1 = r0 + v0*t1
	return r1[1], r1[2], v1[1], v1[2]
end


# Takes the point inside a square formed by the four closest obstacles and determines which wall will be crossed, thus converting the conditions back to x0 = (0, b) to use with first_collision() (with possible turning the coordinates)
function post_collision(x, y, vx, vy, r)
	n = ifloor(x)
	m = ifloor(y)
	
	# Distance from point (x, y) to line y = kx + b
	dist_point_line(x, y, k, b) = abs(y - k*x - b)/sqrt(k^2 + b^2)
	if dist_point_line(n, m, vy/vx, x*vy/vx + y) >= r &&
	   dist_point_line(n, m+1, vy/vx, x*vy/vx + y) >= r &&
	   dist_point_line(n+1, m, vy/vx, x*vy/vx + y) >= r &&
	   dist_point_line(n+1, m+1, vy/vx, x*vy/vx + y) >= r	
	
	
		tv1 = (n - x)/vx 
		tv2 = (n - x + 1)/vx
		th1 = (m - y)/vy 
		th2 = (m - y + 1)/vy
	
		array_times = Real[]
		push!(array_times, tv1, tv2, th1, th2)
	
		# Extract the minimum positive time value
		minpos = Inf
		number = 0
		for i = 1:length(array_times)
			if array_times[i] < minpos && array_times[i] > 0
				minpos = array_times[i]
				number = i
			end				
		end
		# number = 1 or 2 - vertical crossing, 3 or 4 - horizontal
	
		if number == 1
			posx = n
			posy = m
			x1 = 0
			y1 = y + vy/vx*(n - x) - m
		elseif number == 2
			posx = n + 1
			posy = m
			x1 = 0
			y1 = y + vy/vx*(n + 1 - x) - m
		elseif number == 3
			posx = n
			posy = m
			x1 = x + (m - y)*vx/vy - n
			y1 = 0
		elseif number == 4
			posx = n
			posy = m + 1
			x1 = x + (m + 1 - y)*vx/vy - n
			y1 = 0
		end
	
		return x1, y1, posx, posy, number
	# It returns the coordinates (posx, posy) of the bottom-left corner of the unit square and coordinates (x1, y1) within this unit square
	# The exact coords of the particle at crossing are posx+x1, posy+y1
	
	# What if one of the distances to the four obstacles is less than r? Collide and reflect, then recursively apply the function until it goes out of the square
	elseif dist_point_line(n, m, vy/vx, x*vy/vx + y) < r
		x1, y1, vx1, vy1 = collide(n, m, x, y, vx, vy, r)
		post_collision(x1, y1, vx1, vy1, r)
		println("1 more collision")
	elseif dist_point_line(n, m+1, vy/vx, x*vy/vx + y) < r
		x1, y1, vx1, vy1 = collide(n, m+1, x, y, vx, vy, r)
		post_collision(x1, y1, vx1, vy1, r)
		println("1 more collision")		
	elseif dist_point_line(n+1, m, vy/vx, x*vy/vx + y) < r
		x1, y1, vx1, vy1 = collide(n+1, m, x, y, vx, vy, r)
		post_collision(x1, y1, vx1, vy1, r)
		println("1 more collision")		
	elseif dist_point_line(n+1, m+1, vy/vx, x*vy/vx + y) < r
		x1, y1, vx1, vy1 = collide(n+1, m+1, x, y, vx, vy, r)
		post_collision(x1, y1, vx1, vy1, r)
		println("1 more collision")
	end
	
end

function collisions(x, y, vx, vy, r, maxsteps)
	steps = 1
	places = Vector[]
	q, p = first_collision(x, y, vx, vy, abs(r/vx))
	push!(places, [q, p])
	#b = y
	posx = 0; posy = 0
	while steps <= maxsteps
		println("step ", steps)
		
		@show x1, y1, vx1, vy1 = collide(q, p, x, y, vx, vy, r) # And we obtain coords of collision and new velocity
		@show x2, y2, posx, posy, number = post_collision(x1 + posx, y1 + posy, vx1, vy1, r)
		
		if number == 1 || number == 2 # Left/right wall crossed
			q1, p1 = first_collision(x2, y2, vx1, vy1, abs(r/vx1))
			q = q1 + posx
			p = p1 + posy
			@show push!(places, [q, p])
		elseif number == 3 || number == 4 # Bottom/top wall crossed
			# Rotate the axes, x->y, y->x and then interchange q and p
			q, p = first_collision(y2, x2, vy1, vx1, abs(r/vy1))
			#q = q1
			#p = p1
			@show push!(places, [p + posx, q + posy])
			y = x2
			x = y2
			vy = vx1
			vx = vy1
		
		end

				
		@show steps += 1
	end
	places
end


## End of program. Below are the tests and measurements.


# include("atahualpa.jl"); x = 0; y = 0.445; vx = cos(1); vy = sin(1); r = 0.1
# include("atahualpa.jl"); x = 0; y = 0.35; vx = cos(1); vy = sin(1); r = 0.1


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
r = 0.001
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

#=
x = 0
y = 0.445


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

# After changing the frac() function into Atahualpa's one

# with println: 62.02397 s
# without println: 20.11622 s
# deg 1-44: 11.938 s
# deg 46-89: 7.825 s

#=
function measuring2()
	r = 0.000001
	for deg = 46:89
		if deg != 45
			vx = cos(deg*pi/180)
			vy = sin(deg*pi/180)
			#println(r, " ", deg, " ", first_collision(x, y, vx, vy, r/vx))
			first_collision(x, y, vx, vy, r/vx)
		end
	end
end

@time measuring2()
=#
