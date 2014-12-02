# Perpendicular (z-) component of cross product (scalar quantity) for 2-D vectors: (x cross y) dot e_z
crossz(x, y) = x[1]*y[2] - x[2]*y[1]

function crossing(r, v, n, m)
	# Outputs new position in square [-0.5, 0.5)^2 and updates the numbers of new square
	
	# Minimum positive times to nearest crossing - concise formulas p. 48 - 49
	array_times = Real[]
	k = v[2]/v[1]	
	for i = 1:2
		for j = 1:2
			push!(array_times, ((-1)^j*0.5 - r[i])/v[i])
		end	
	end
	
	# Extract the minimum positive time value
	minpos = Inf
	number = 0
	for i = 1:length(array_times)
		if array_times[i] < minpos && array_times[i] > 0
			minpos = array_times[i]
			number = i
		end				
	end
	
	t = array_times[number]
	r1 = r + v*t
	d = norm(v*t)
	
	# Possible outcomes depending on which wall will be hit
	if number == 2
		n += 1
		return [-0.5, r1[2]], d, n, m
	elseif number == 4
		m += 1
		return [r1[1], -0.5], d, n, m
	elseif number == 1
		n -= 1
		return [0.5, r1[2]], d, n, m
	elseif number == 3
		m -= 1
		return [r1[1], 0.5], d, n, m
	elseif error("This should not happen")
	end
	
end



function collisions(r0::Vector, v0::Vector, rho::Real, tmax::Real, precision::Integer=64)

	set_bigfloat_precision(precision)
	places = Vector[]
	circles = Vector[]
	# Put the starting point into the places array
	push!(places, r0)	

	# Initial square (n, m)
	n = floor(r0[1] + 0.5)
	m = floor(r0[2] + 0.5)
	
	# Place the first initial position into square [-0.5, 0.5)^2
	r0 -= [n, m]
	
	if norm(r0) < rho
		error("The initial position cannot be inside an obstacle")
	end

	t = 0
	while t <= tmax
		# Will hit or miss? Check the condition for hitting
		if abs(crossz(v0, r0)) < norm(v0)*rho
			#println("hit")
			# Then reflect
			
			discr = norm(v0)^2*rho^2 - (crossz(v0, r0))^2
			# Throw away complex time values if any, but there shouldn't be
			if discr >= 0
				t1 = (-dot(v0, r0) - sqrt(discr))/norm(v0)^2
			end
		
			N0 = r0 + v0*t1
			N = N0/norm(N0)
			
			# Velocity, place and time immediately after the collision
			v1 = v0 - 2*dot(v0, N)*N
			r1 = r0 + v0*t1
			t += t1
			
			push!(places, r1 + [n, m])
			push!(circles, [n, m])
			#print("{$(r1[1] + n), $(r1[2] + m)}, ")
			
			#r1 = [mod(r1[1], 1) - 0.5, mod(r1[2], 1) - 0.5]
			r0, d, n, m = crossing(r1, v1, n, m)
						
			# The speed direction will stay the same
			v0 = v1
			
			# The time will increment once again from collision point to the wall
			t += d/norm(v0)

			
		else # If it misses the ball
			#println("miss")
			v1 = v0
			r1 = r0
			
			# Now hit the wall
			r0, d, n, m = crossing(r1, v1, n, m)
						
			# The speed direction will stay the same
			v0 = v1
			
			# The time will increment
			t += d/norm(v0)
		end
	end
	return places, circles
end


#println(collisions([0.4, 0.1], [-0.42, 0.23], 0.3, 20)[1])
