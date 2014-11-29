
# Perpendicular (z-) component of cross product (scalar quantity) for 2-D vectors: (x cross y) dot e_z
crossz(x, y) = x[1]*y[2] - x[2]*y[1]

function collisions(r0, v0, rho, tmax)

	places = Vector[]
	circles = Vector[]
	push!(places, r0)	
	#Normalize the speed
	v0 = v0/norm(v0)

	t = 0
	while t <= tmax

		# Arrays of times, n and m corresponding to each intersection in both solutions

		array_t1nm = zeros(Float64, 3000, 3)

		# Checking the real intersection points

		i = 1
		for n = -500:500
			for m = -500:500
				R = [n,m]
				# By trial and error, it turned out that only the solution with a "-" before sqrt gave valid data
				t1 = -dot(v0, r0 - R) - sqrt(complex(rho^2 - (crossz(v0, r0 - R))^2))

				if abs(crossz(v0, r0 - R)) < rho
					array_t1nm[i, 1] = t1
					array_t1nm[i, 2] = n
					array_t1nm[i, 3] = m
					i = i + 1
				end
			end
		end

		# Extracting the minimum positive time from array

		tmin1 = Inf
		n1 = 0
		m1 = 0
		for i = 1:size(array_t1nm)[1]-1
			if (array_t1nm[i, 1] < tmin1) && (array_t1nm[i, 1] > 0)
				tmin1 = array_t1nm[i, 1]
				n1 = array_t1nm[i, 2]
				m1 = array_t1nm[i, 3]
			end
		end

		tfinal = tmin1
		nfinal = n1
		mfinal = m1

		# Normal vector at the point of collision
		N = (r0 + v0*tfinal - [nfinal, mfinal])/norm(r0 + v0*tfinal - [nfinal, mfinal])

		# New velocity
		v1 = v0 - 2*dot(v0, N)*N

		# New position
		r1 = r0 + v0*tfinal

		# New time
		t = t + tfinal

		# Taking on the new values of v and r
		v0 = v1
		r0 = r1

		#println("r0 = ", r0, "v0 = ", v0, "t = ", t)
		
		push!(places, r0)
		push!(circles, [nfinal, mfinal])
		# Output for Mathematica
		print("{$(r0[1]), $(r0[2])}, ")

	end
	return places, circles
	# Now the location will be

	#r = r0 + v0*(10 - t)
	#d = norm(r)
	#println("Final: r = ", r, "d = ", d)

end

#collisions([0.5, 0.1], [0.7, 0], 0.35, 20)
