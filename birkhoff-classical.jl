# In all configuration it is assumed that radius of the disk R = 1

#hit(init.pos., init.vel., dist.betw.centers) = final pos., final vel., time
function hit(x, v, r)
    v /= norm(v)
    c12 = [r, 0]
    c13 = r/2*[1, sqrt(3)]
	B12 = dot(x - c12, v)
	B13  = dot(x - c13, v)
	C12 = r*(r - 2x[1])
	C13 = r*(r - x[1] - sqrt(3)*x[2])
	if B12^2 - C12 >= 0
	    #hits second disk
	    t1 = -B12 - sqrt(B12^2 - C12)
	elseif B13^2 - C13 >= 0
	    #hits third disk
	    t1 = -B13 - sqrt(B13^2 - C13)
	else
	    #misses and flies off to infinity
	    t1 = Inf
	end
	if t1 < Inf # If collides
	    N0 = x + v*t1
	    N = N0/norm(N0)
	    # Velocity, place and time immediately after the collision
	    v1 = v - 2dot(v, N)*N
	    x1 = x + v*t1
	    #t += t1
	    return x1, v1, t1
	else
	    return [Inf, Inf], v, Inf
	end
end

theta = pi/6
x = [cos(theta), sin(theta)]
v = [1, 0]
# Hit 2nd disk 
hit2 = hit(x, v, 6)
# Hit 3rd disk
hit3 = hit(x, [0.5, sqrt(3)/2], 6)
# Miss
miss = hit(x, [1, 1], 6)
