## Automatic differentiation (interval version)

module AutoDiff
export Ad, det2, differentiate, jacobian

type Ad
	u
	up
end

Ad(c) = Ad(c, zero(c))

# Arithmetic between two Ad
import Base.+, Base.-, Base.*, Base./
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
import Base.^
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
	(:(Base.asin), :(x->1/sqrt(1 - x^2))), 
	(:(Base.acos), :(x->-1/sqrt(1 - x^2))),
	(:(Base.atan), :(x->1/(1 + x^2))),  
	(:(Base.abs), :(Base.sign)),
	(:(Base.sqrt), :(x->1/(2sqrt(x))))				
	)
    @eval begin
        function $func(x::Ad)
            Ad($func(x.u), x.up*$deriv(x.u))
        end
    end
end

# Some more functions

e^(x::Ad) = Ad(e^x.u, x.up*e^x.u)

import IntervalArithmetic.arcsin
arcsin(x::Ad) = Ad(arcsin(x.u), x.up/sqrt1(1 - x.u^2))

import IntervalArithmetic.sqrt1
sqrt1(x::Ad) = Ad(sqrt1(x.u), x.up/(2sqrt1(x.u)))

import IntervalArithmetic.arcsin_d
arcsin_d(x::Ad) = Ad(arcsin_d(x.u), x.up/sqrt_d(1 - x.u^2))

import IntervalArithmetic.sqrt_d
sqrt_d(x::Ad) = Ad(sqrt_d(x.u), x.up/(2sqrt_d(x.u)))


import Base.complex
complex(x::Ad) = Ad(x.u, x.up)

import Base.one
one(x::Ad) = Ad(x.u, x.up)

import Base.rem
rem(x::Ad, y::Real) = Ad(rem(x.u, y), x.up)

import Base.mod
mod(x::Ad, y::Real) = Ad(mod(x.u, y), x.up)

import Base.mod2pi
mod(x::Ad) = Ad(mod(x.u), x.up)

import Base.mod1 # mod1 is a built-in function
mod1(x::Ad, y::Real) = Ad(mod1(x.u, y), x.up)

mod2(x::Ad, y::Real) = Ad(mod2(x.u, y), x.up)

mod21(x::Ad, y::Real) = Ad(mod21(x.u, y), x.up)
mod22(x::Ad, y::Real) = Ad(mod22(x.u, y), x.up)
mod23(x::Ad, y::Real) = Ad(mod23(x.u, y), x.up)
mod24(x::Ad, y::Real) = Ad(mod24(x.u, y), x.up)

# Differentiating operator (1D and 2D)

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

