using KrawczykMethod2D

# The vector x = [w, th] in Birkhoff mapping

function alpha(n, m)
	if n > 0 && n <= 3 && m > 0 && m <= 3 && typeof(n) == Int && typeof(m) == Int && n != m
		k = pi/3
		if n == 1 && m == 2
			return 0
		elseif n == 1 && m == 3
	 		return k
		elseif n == 2 && m == 3
	 		return 2k
		elseif n == 2 && m == 1
	 		return 3k
		elseif n == 3 && m == 1
	 		return 4k
		elseif n == 3 && m == 2
	 		return 5k
	 	end
	else 
	    error("Invalid parameters: should be 1, 2, 3 and not equal")
	end
end

limit1(x) = x[2] >= 0 && x[2] <= pi/3

w_lo(th, r) = sin(th - atan((sqrt(3)/2 - sin(th))/(r - 1 - cos(th))))
w_hi(th, r) = sin(th + atan(sin(th)/(r - 1 - cos(th))))

function limit2(x, r)
	x[1] >= w_lo(x[2], r) && x[1] <= w_hi(x[2], r)
end

function T(x, n, m, r)
	wnext = x[1] - r*(x[1]*cos(x[2] - alpha(n, m)) + sqrt1(1 - x[1]^2)*sin(x[2] - alpha(n, m)))
	[wnext, mod(x[2] + float(pi) + arcsin(x[1]) + arcsin(wnext), 2pi)]
end

f(x, r) = T(T(T(x, 1, 2, r), 2, 3, r), 3, 1, r) - x

g(x) = f(x, 3)

f1(x, r) = T(T(x, 1, 3, r), 3, 1, r) - x
g1(x) = f1(x, 3)

#=
for i = 1:1e9
	dth = float(pi)/(3*1e9)
	theta = i*dth
	@show Interval(w_lo(theta + dth, 3), w_hi(theta, 3))
	krawczyk2d(g, [Interval(w_lo(theta + dth, 3), w_hi(theta, 3)), Interval(theta, theta + dth)])
end
=#

#= The problem is that asin() returns complex values in the rectangle (-1,1)*(0,pi/3) and krawczyk2d cannot work like that. It solves for the periodic orbit 1-2-3 only in a very narrow range: e.g., 
krawczyk2d(g, [Interval(-0.5-0.00001,-0.5+0.00001), Interval(pi/6-0.00001,pi/6+0.00001)])
as w=-0.5, th=pi/6 is the obvious answer.

(-0.5, 0.5235987755982988)

=#

#krawczyk2d_general(g, [Interval(-0.5-0.01,-0.5+0.021), Interval(pi/6-0.01,pi/6+0.021)])


# Solves krawczyk2d(g, [Interval(-0.5-0.0001,-0.5+0.00021), Interval(pi/6-0.001,pi/6+0.00021)]) - domaincheck returns 1


#=
Reproduce 3-element issue:

x = [Interval(-6.45723568227448264749e-01,-5.1962031122224282807e-01), Interval(4.548837360587528278e+00,4.74163346223171781532e+00)]  
T(ans, 3, 1, 3)


x = [Interval(-4.97340331337152747313e-01,-4.97340331336701637922e-01),Interval(5.29098775596043253229e-01,5.29098775596494362702e-01)]

ftest1(x) = [arcsin_d(x[1] - x[2]) - float(pi)/2, sqrt_d(x[1] + x[2]) - 3]
krawczyk2d_general(ftest1, [Interval(-3, 8), Interval(-2, 7)])
krawczyk2d(ftest1, [Interval(-3, 8), Interval(-2, 7)])

=#
