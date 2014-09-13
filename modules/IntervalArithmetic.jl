## Interval arithmetic

module IntervalArithmetic
export Interval, rad, diam, mid, mig, mag, belong, hd, hull, isect, isectext

typealias prec BigFloat

type Interval

    lo
    hi

    function Interval(a, b)
        set_rounding(prec, RoundDown)
        lo = BigFloat("$a")

        set_rounding(prec, RoundUp)
        hi = BigFloat("$b")

        new(lo, hi)
    end

end

# Thin (degenerate) interval and functions zero() and one()

Interval(x) = Interval(x, x)

import Base.one
one(x::Interval) = Interval(1.0)

import Base.zero
zero(x::Interval) = Interval(0.0)

# The basic operations on intervals. Left end is rounded down, right end is rounded up.

# Addition

function +(x::Interval, y::Interval)
    z1 = with_rounding(prec, RoundDown) do
		x.lo + y.lo
    end
    z2 = with_rounding(prec, RoundUp) do
		x.hi + y.hi
    end
    Interval(z1, z2)
end

+(x::Interval, y::Real) = x + Interval(y)
+(x::Real, y::Interval) = Interval(x) + y

# Subtraction

function -(x::Interval, y::Interval)
    z1 = with_rounding(prec, RoundDown) do
x.lo - y.hi
    end
    z2 = with_rounding(prec, RoundUp) do
x.hi - y.lo
    end
    Interval(z1, z2)
end

-(x::Interval) = Interval(-x.lo, -x.hi)

-(x::Interval, y::Real) = x - Interval(y)
-(x::Real, y::Interval) = Interval(x) - y

# Multiplication

function *(x::Interval, y::Interval)
    z1 = with_rounding(prec, RoundDown) do
		min(x.lo*y.lo, x.lo*y.hi, x.hi*y.lo, x.hi*y.hi)
    end
    z2 = with_rounding(prec, RoundUp) do
		max(x.lo*y.lo, x.lo*y.hi, x.hi*y.lo, x.hi*y.hi)
    end
    Interval(z1, z2)
end

*(x::Interval, y::Real) = x*Interval(y)
*(x::Real, y::Interval) = Interval(x)*y


# Division

function /(x::Interval, y::Interval)
    z1 = with_rounding(prec, RoundDown) do
		1/y.hi
    end
    z2 = with_rounding(prec, RoundUp) do
		1/y.lo
    end
    x*Interval(z1, z2)
end

/(x::Interval, y::Real) = x/Interval(y)
/(x::Real, y::Interval) = Interval(x)/y

# Extended division

function //(x::Interval, y::Interval)
	if belong(0., y) == false
		return x/y
	elseif belong(0., x) == true && belong(0, y) == true
		return Interval(-Inf, Inf)
	elseif x.hi < 0 && y.lo < y.hi == 0
		z1 = with_rounding(prec, RoundDown) do
			x.hi/y.lo
    	end
    	return Interval(z1, Inf)
	elseif x.hi < 0 && y.lo < 0 < y.hi
		z1 = with_rounding(prec, RoundDown) do
			x.hi/y.lo
    	end
    	z2 = with_rounding(prec, RoundUp) do
			x.hi/y.hi
    	end
		return Interval(z1, z2)
	elseif x.hi < 0 && 0 == y.lo < y.hi
		z2 = with_rounding(prec, RoundUp) do
			x.hi/y.hi
    	end
    	return Interval(-Inf, z2)
    elseif 0 < x.lo && y.lo < y.hi == 0
		z2 = with_rounding(prec, RoundUp) do
			x.lo/y.lo
    	end
    	return Interval(-Inf, z2)
	elseif 0 < x.lo && y.lo < 0 < y.hi
		z1 = with_rounding(prec, RoundDown) do
			x.lo/y.hi
		end
		z2 = with_rounding(prec, RoundUp) do
			x.lo/y.lo
		end
		return Interval(z1, z2)
	elseif 0 < x.lo && 0 == y.lo < y.hi
		z1 = with_rounding(prec, RoundDown) do
			x.lo/y.hi
		end
		return Interval(z1, Inf)
	elseif belong(0, x) == false && y.lo == 0 && y.hi == 0
		return println("Empty set: extended division by thin zero")
	end
end		
	

## Interval properties

# Whether the point belongs to the interval: belong(point, interval)

function belong(p::Real, x::Interval)
    if p >= x.lo && p <= x.hi
		return true
    else return false
    end
end

# Radius, diameter, midpoint, mignitude, magnitude, absolute value

rad(x::Interval) = (x.hi - x.lo)/2
diam(x::Interval) = x.hi - x.lo
mid(x::Interval) = (x.hi + x.lo)/2
# Making mid() process 1-D interval arrays into arrays of midpoints
function mid(array::Array{Interval, 1})
	array1 = prec[]
	for i = 1:length(array)
		push!(array1, mid(array[i]))
	end
	array1
end

function mig(x::Interval)
    if belong(0.0, x) == true
		return 0
    else return min(abs(x.lo), abs(x.hi))
    end
end

mag(x::Interval) = max(abs(x.lo), abs(x.hi))

import Base.abs
abs(x::Interval) = Interval(mig(x), mag(x))


# Hausdorff distance

hd(x::Interval, y::Interval) = max(abs(x.lo - y.lo), abs(x.hi - y.hi))


# "Union" (hull) and intersection

hull(x::Interval, y::Interval) = Interval(min(x.lo, y.lo), max(x.hi, y.hi))

function isect(x::Interval, y::Interval)
    if x.hi < y.lo || y.hi < x.lo
		return false
    else 
	z1 = with_rounding(prec, RoundDown) do
		max(x.lo, y.lo)
    end
    z2 = with_rounding(prec, RoundUp) do
		min(x.hi, y.hi)
    end 
    return Interval(z1, z2)
    end
end

# Extended intersection (involving extended intervals [a, b] with a > b)

function isectext(x::Interval, y::Interval)
	if x.hi < x.lo && y.hi >= y.lo # x is an extended interval and y is a normal one
		if y.lo <= x.hi && x.lo <= y.hi
			return [Interval(y.lo, x.hi), Interval(x.lo, y.hi)]
		elseif y.lo > x.hi && x.lo <= y.hi && y.lo <= x.lo
			return Interval(x.lo, y.hi)
		elseif y.lo <= x.hi && x.lo > y.hi && y.hi >= x.hi
			return Interval(y.lo, x.hi)
		elseif y.lo > x.hi && x.lo > y.hi
			return false
		elseif y.lo <= x.hi && x.lo > y.hi && y.hi < x.hi
			return y
		elseif y.lo > x.hi && x.lo <= y.hi && y.lo > x.lo
			return y
		end
	elseif x.hi < x.lo && y.hi < y.lo # both intervals are extended
		return Interval(max(x.lo, y.lo), min(x.hi, y.hi)) # Returns also an extended interval
	elseif x.hi >= x.lo && y.hi < y.lo # x normal, y extended
		return isectext(y, x)
	elseif x.hi >= x.lo && y.hi >= y.lo # both intervals are normal
		return isect(x, y)
	end
end


# Integer power

function ^(x::Interval, n::Int32)
    if n > 0 && n % 2 == 1
		return Interval(x.lo^n, x.hi^n)
    elseif n > 0 && n % 2 == 0
		return Interval((mig(x))^n, (mag(x))^n)
    elseif n == 0
		return Interval(1, 1)
    elseif n < 0 && belong(0, x) == false
		return Interval(1/x.hi, 1/x.lo)^(-n)
	# elseif return println("Error")
    end
end


# Real power - exercise 3.5 from Tucker "Validated Numerics"

function ^(x::Interval, n::Interval)
	if x.lo > 0
		z1 = min(x.lo^n.lo, x.lo^n.hi, x.hi^n.lo, x.hi^n.hi)
		z2 = max(x.lo^n.lo, x.lo^n.hi, x.hi^n.lo, x.hi^n.hi)
		return Interval(z1, z2)
	end
end


# Trigonometry

#= COMMENTED OUT: old Tucker definition of sin
import Base.sin
function sin(x::Interval)

    k = 0
    pcount = 0
    for k = -1000:1000
p = pi/2 + 2pi*k
if belong(p, x) == true
pcount = pcount + 1
end
    end
    
    k = 0
    qcount = 0
    for k = -1000:1000
q = - pi/2 + 2pi*k
if belong(q, x) == true
qcount = qcount + 1
end
    end
    
    if qcount != 0 && pcount != 0
		return Interval(-1., 1.)
    elseif qcount != 0 && pcount == 0
		return Interval(-1., max(sin(x.lo), sin(x.hi)))
    elseif qcount == 0 && pcount != 0
		return Interval(min(sin(x.lo), sin(x.hi)), 1.)
    elseif qcount == 0 && pcount == 0
		return Interval(min(sin(x.lo), sin(x.hi)), max(sin(x.lo), sin(x.hi)))
    end
end
=#


# Following http://jenchienjackchang.com/sample-page/implicit-solid-modeling-using-interval-methods/interval-arithmetic/


import Base.sin
function sin(x::Interval)
	if x.lo%2pi >= 0
		low = x.lo%2pi
	else low = x.lo%2pi + 2pi
	end	
	x1 = Interval(low, low + diam(x))
	# If the interval has a diameter equal or greater than 2pi or it is an extended interval, return [-1, 1]
	if diam(x) >= 2pi || x.lo > x.hi
		return Interval(-1, 1)
	elseif 0 <= x1.lo <= x1.hi <= pi/2 || 3pi/2 <= x1.lo <= x1.hi <= 5pi/2 || 7pi/2 <= x1.lo <= x1.hi <= 4pi
			return Interval(sin(x1.lo), sin(x1.hi))
		elseif pi/2 <= x1.lo <= x1.hi <= 3pi/2 || 5pi/2 <= x1.lo <= x1.hi <= 7pi/2
			return Interval(sin(x1.hi), sin(x1.lo))
		elseif (0 <= x1.lo <= pi/2 && pi/2 <= x1.hi <= 3pi/2) || (3pi/2 <= x1.lo <= 5pi/2 && 5pi/2 <= x1.hi <= 7pi/2)
			return Interval(min(sin(x1.lo), sin(x1.hi)), 1)
		elseif (pi/2 <= x1.lo <= 3pi/2 && 3pi/2 <= x1.hi <= 5pi/2) || (5pi/2 <= x1.lo <= 7pi/2 && 7pi/2 <= x1.hi <= 4pi)
			return Interval(-1, max(sin(x1.lo), sin(x1.hi)))
		elseif (0 <= x1.lo <= pi/2 && 3pi/2 <= x1.hi <= 5pi/2) || (pi/2 <= x1.lo <= 3pi/2 && 5pi/2 <= x1.hi <= 7pi/2) || (3pi/2 <= x1.lo <= 5pi/2 && 7pi/2 <= x1.hi <= 4pi)
			return Interval(-1, 1)
		end
	
end


import Base.cos
cos(x::Interval) = sin(Interval(x.lo + pi/2, x.hi + pi/2))

# Miscellaneous

import Base.complex
complex(x::Interval) = Interval(complex(x.lo), complex(x.hi))

function *(x::Array{Interval, 1}, y::Interval)
    [x[1]*y, x[2]*y]
end


function /(x::Array{Any, 1}, y::Interval)
    [x[1]/y, x[2]/y]
end

# Equality of two intervals
==(a::Interval, b::Interval) = a.lo == b.lo && a.hi == b.hi

# Monotonic functions

import Base.exp
exp(x::Interval) = Interval(exp(x.lo), exp(x.hi))

import Base.sqrt
sqrt(x::Interval) = Interval(sqrt(x.lo), sqrt(x.hi))

import Base.log
log(x::Interval) = Interval(log(x.lo), log(x.hi))

import Base.asin
asin(x::Interval) = Interval(asin(x.lo), asin(x.hi))

import Base.acos
acos(x::Interval) = Interval(acos(x.hi), acos(x.lo))

# End of module
end
