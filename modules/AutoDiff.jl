## Automatic differentiation (interval version)

module AutoDiff
	export Ad, det2, differentiate, jacobian

	type Ad
		u
		up
	end

	Ad(c) = Ad(c, zero(c))

	# Arithmetic between two Ad
	+(x::Ad, y::Ad) = Ad(x.u + y.u, x.up + y.up)
	-(x::Ad, y::Ad) = Ad(x.u - y.u, x.up - y.up)
	*(x::Ad, y::Ad) = Ad(x.u*y.u, x.u*y.up + y.u*x.up)
	function /(x::Ad, y::Ad)
		ff = x.u/y.u
		dd = (x.up - ff*y.up)/y.u
		Ad(ff, dd)
	end

	# Power; two functions for power because it is less costly
	# Commented out because it gave errors
	^(x::Ad, a::Integer)= Ad(x.u^a, a*x.u^(a-1)*x.up)

	function ^(x::Ad, a::Real)
		ff = x.u^(a-1)    
		Ad(ff*x.u, a*ff*x.up)
	end


	# Arithmetic operations between Ad and intervals/numbers

	for op in (:+, :-, :*, :/)
		@eval begin
		    $op(u::Ad, c) = $op(u, Ad(c))
		    $op(c, u::Ad) = $op(Ad(c), u)
		end
	end

	+(x::Ad) = x
	-(x::Ad) = Ad(-x.u, -x.up)

	# Elementary functions - now with metaprogramming
	
	for (func, deriv) in (
		(:(Base.exp), :(Base.exp)),  
		(:(Base.log), :(x->1/x)), 
		(:(Base.sin), :(Base.cos)), 
		(:(Base.cos), :(x->-sin(x))), 
		(:(Base.tan), :(x->(sec(x))^2)),
		(:(Base.abs), :(Base.sign)),
		(:(Base.sqrt), :(x->1/(2sqrt(x))))				
		)
	    @eval begin
	        function $func(x::Ad)
	            Ad($func(x.u), x.up*$deriv(x.u))
	        end
	    end
	end

	e^(x::Ad) = Ad(e^x.u, x.up*e^x.u)
	
	import Base.rem
	rem(x::Ad, y::Real) = Ad(rem(x.u, y), one(y))
	
	import Base.mod
	mod(x::Ad, y::Real) = Ad(mod(x.u, y), one(y))
	
	#=
	import Base.sin
	sin(x::Ad) = Ad(sin(x.u), x.up*cos(x.u))

	import Base.cos
	cos(x::Ad) = Ad(cos(x.u), -x.up*sin(x.u))

	import Base.exp
	exp(x::Ad) = Ad(exp(x.u), x.up*exp(x.u))

	import Base.log
	log(x::Ad) = Ad(log(x.u), x.up/x.u)

	import Base.abs
	abs(x::Ad) = Ad(abs(x.u), x.up*sign(x.u))
	
	import Base.sqrt
	sqrt(x::Ad) = Ad(sqrt(x.u), x.up/(2sqrt(x.u)))
	
	=#
	
	function differentiate(f, a) 
		y = f(Ad(a, one(a)))
		if typeof(y) != Ad # When f(x) = const so f'(x) = 0, and typeof(f(Ad())) = typeof(const)
			return zero(a)
		else
			return y.up
		end
	end

	function jacobian(f, a)
		f1(x) = f(x)[1]
		f2(x) = f(x)[2]

		f11(x1) = f1([x1, a[2]])
		J11 = differentiate(f11, a[1])

		f12(x2) = f1([a[1], x2])
		J12 = differentiate(f12, a[2])

		f21(x1) = f2([x1, a[2]])
		J21 = differentiate(f21, a[1])

		f22(x2) = f2([a[1], x2])
		J22 = differentiate(f22, a[2])

		[J11 J12; J21 J22]
	end

# End of module
end

