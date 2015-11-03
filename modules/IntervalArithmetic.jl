## Interval arithmetic

module IntervalArithmetic
export Interval, ComplexInterval, MultiDimInterval, IntUnion, rad, diam, mid, mig, mag, belong, hd, hull, isect, isectext, lo, hi, left, right, make_intervals, all_inside, bisect, det2, inside, intunion, mod1, mod2, mod21, mod22, mod23, mod24, domaincheck, domaincheck2d, arcsin, sqrt1, flip, arcsin_d, sqrt_d

typealias prec BigFloat

type Interval

	lo
	hi

	function Interval(a, b)
	    set_rounding(prec, RoundDown)
	    #lo = BigFloat("$a")
	    lo = parse(BigFloat, "$a")

	    set_rounding(prec, RoundUp)
	    #hi = BigFloat("$b")
	    hi = parse(BigFloat, "$b")

	    new(lo, hi)
	end

end


typealias MultiDimInterval Array{Interval, 1}

#import MPFR.BigFloat
#BigFloat(x::Interval) = Interval(BigFloat(x.lo), BigFloat(x.hi))

# Thin (degenerate) interval and functions zero() and one()

Interval(x::Number) = Interval(x, x)

# Prevent errors from embedded interval functions
Interval(x::Interval) = x

import Base.one
one(x::Interval) = Interval(1.0)

import Base.zero
zero(x::Interval) = Interval(0.0)
zero(Interval) = Interval(0)

# The basic operations on intervals. Left end is rounded down, right end is rounded up.

# Addition
import Base.+, Base.-
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

-(x::Interval) = Interval(-x.hi, -x.lo)

-(x::Interval, y::Real) = x - Interval(y)
-(x::Real, y::Interval) = Interval(x) - y

# Multiplication

import Base.*
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

import Base./
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

import Base.//
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
		return error("\nEmpty set: extended division by thin zero")
	end
end


## Interval properties

# Whether the point belongs to the interval: belong(point, interval)

function belong(p::Real, x::Interval)
	if p >= x.lo && p <= x.hi
		return true
	else
		return false
	end
end

# Whether one interval is inside of the other

function inside(x::Interval, y::Interval)
	if x.lo >= y.lo && x.hi <= y.hi
		return true
	else
		return false
	end
end

inside(x::Real, y::Interval) = belong(x, y)
inside(x::Interval, y::Real) = belong(y, x)

# Miscellaneous: lower end, higher end, lower (left) half, higher (right) half, bottom, top, radius, diameter, midpoint, mignitude, magnitude, absolute value

lo(x::Interval) = x.lo
hi(x::Interval) = x.hi
left(x::Interval) = Interval(x.lo, mid(x))
right(x::Interval) = Interval(mid(x), x.hi)
rad(x::Interval) = (x.hi - x.lo)/2
diam(x::Interval) = x.hi - x.lo
mid(x::Interval) = (x.hi + x.lo)/2
mid(x::prec) = x

mid(x::Real) = x
lo(x::Real) = x
hi(x::Real) = x
diam(x::Real) = 0
rad(x::Real) = 0


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

isect(x::Real, y::Interval) = isect(Interval(x), y)
isect(x::Interval, y::Real) = isect(x, Interval(y))

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

import Base.^
function ^(x::Interval, n::Integer)
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


# Real power function - taken from https://github.com/dpsanders/ValidatedNumerics.jl
macro round_down(expr)
	quote
		with_rounding(BigFloat, RoundDown) do
				$expr
			end
	end
end

macro round_up(expr)
	quote
	with_rounding(BigFloat, RoundUp) do
			$expr
		end
	end
end

macro round(expr1, expr2)
	quote
		Interval(@round_down($expr1), @round_up($expr2))
	end
end

function reciprocal(a::Interval)
	uno = one(BigFloat)
	z = zero(BigFloat)
	if belong(z, a)
	#if z in a
		warn("\nInterval in denominator contains 0.")
		return Interval(-inf(z),inf(z)) # inf(z) returns inf of type of z
	end
	@round(uno/a.hi, uno/a.lo)
end

function ^(a::Interval, x::Real)
	x == int(x) && return a^(int(x))
	x < zero(x) && return reciprocal(a^(-x))
	x == 0.5*one(x) && return sqrt(a)
	z = zero(BigFloat)
	z > a.hi && error("Undefined operation; Interval is strictly negative and power is not an integer")
	xInterv = Interval(x)
	diam(xInterv) >= eps(x) && return a^xInterv
	# xInterv is a thin interval
	domainPow = Interval(z, big(Inf))
	aRestricted = isect(a, domainPow)
	@round(aRestricted.lo^x, aRestricted.hi^x)
end


# Interval power - exercise 3.5 from Tucker "Validated Numerics"

function ^(x::Interval, n::Interval)
	if x.lo > 0
		z1 = min(x.lo^n.lo, x.lo^n.hi, x.hi^n.lo, x.hi^n.hi)
		z2 = max(x.lo^n.lo, x.lo^n.hi, x.hi^n.lo, x.hi^n.hi)
		return Interval(z1, z2)
	end
end


# Trigonometry


# Taken from http://jenchienjackchang.com/sample-page/implicit-solid-modeling-using-interval-methods/interval-arithmetic/

#= Deprecated in favour of the next
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
=#

# Using Sanders/Benet sin() from https://github.com/dpsanders/ValidatedNumerics.jl
import Base.sin
function sin(a::Interval)
	if a == Interval(Inf, Inf)
		return Interval(Inf, Inf)
	end
	#piHalf = pi*BigFloat("0.5")
	piHalf = pi*parse(BigFloat, "0.5")
	#twoPi = pi*BigFloat("2.0")
	twoPi = pi*parse(BigFloat, "2.0")
	domainSin = Interval( BigFloat(-1.0), BigFloat(1.0) )

	# Checking the specific case
	diam(a) >= twoPi && return domainSin

	# Limits within 1 full period of sin(x)
	# Abbreviations
	loMod2pi = mod(a.lo, twoPi)
	hiMod2pi = mod(a.hi, twoPi)
	loQuartile = floor( loMod2pi / piHalf )
	hiQuartile = floor( hiMod2pi / piHalf )

	# 20 different cases
	if loQuartile == hiQuartile # Interval limits in the same quartile
	    loMod2pi > hiMod2pi && return domainSin
	    set_rounding(BigFloat, RoundDown)
	    lo = sin( a.lo )
	    set_rounding(BigFloat, RoundUp)
	    hi = sin( a.hi )
	    set_rounding(BigFloat, RoundNearest)
	    return Interval( lo, hi )
	elseif loQuartile == 3 && hiQuartile==0
	    set_rounding(BigFloat, RoundDown)
	    lo = sin( a.lo )
	    set_rounding(BigFloat, RoundUp)
	    hi = sin( a.hi )
	    set_rounding(BigFloat, RoundNearest)
	    return Interval( lo, hi )
	elseif loQuartile == 1 && hiQuartile==2
	    set_rounding(BigFloat, RoundDown)
	    lo = sin( a.hi )
	    set_rounding(BigFloat, RoundUp)
	    hi = sin( a.lo )
	    set_rounding(BigFloat, RoundNearest)
	    return Interval( lo, hi )
	elseif ( loQuartile == 0 || loQuartile==3 ) && ( hiQuartile==1 || hiQuartile==2 )
	    set_rounding(BigFloat, RoundDown)
	    slo = sin( a.lo )
	    shi = sin( a.hi )
	    set_rounding(BigFloat, RoundNearest)
	    lo = min( slo, shi )
	    return Interval( lo, BigFloat(1.0) )
	elseif ( loQuartile == 1 || loQuartile==2 ) && ( hiQuartile==3 || hiQuartile==0 )
	    set_rounding(BigFloat, RoundUp)
	    slo = sin( a.lo )
	    shi = sin( a.hi )
	    set_rounding(BigFloat, RoundNearest)
	    hi = max( slo, shi )
	    return Interval( BigFloat(-1.0), hi )
	elseif ( loQuartile == 0 && hiQuartile==3 ) || ( loQuartile == 2 && hiQuartile==1 )
	    return domainSin
	else
	    # This should be never reached!
	    error(string("SOMETHING WENT WRONG in sin.\nThis should have never been reached\nArgument: $a") )
	end
end

import Base.cos
cos(x::Interval) = sin(Interval(x.lo + pi/2, x.hi + pi/2))

# Miscellaneous

import Base.complex
complex(x::Interval) = Interval(complex(x.lo), complex(x.hi))
complex(x::MultiDimInterval) = [complex(x[1]), complex(x[2])]

function *(x::MultiDimInterval, y::Interval)
	[x[1]*y, x[2]*y]
end


function /(x::Array{Any, 1}, y::Interval)
	[x[1]/y, x[2]/y]
end

# Equality of two intervals
import Base.==
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


# Old function for the first version of purity and loose evaluation

#= Check if x is within the range of function f - unnecessary, domaincheck does a better job
function is_allowed(f, x)
    try f(x)
    catch DomainError
        return false
    end
    return true
end
=#

function arcsin(x)
    domain = Interval(-1, 1)
    if isect(x, domain) == false
        throw(DomainError())
    elseif inside(x, domain)
        return asin(x)
    else
        throw(ArgumentError(""))
    end
end

function sqrt1(x)
    domain = Interval(0, Inf)
    if isect(x, domain) == false
        throw(DomainError())
    elseif inside(x, domain)
        return sqrt(x)
    else
        throw(ArgumentError(""))
    end
end

function domaincheck(f, x) # 1 - clean, 0 - unclean, -1 - dirty, for the functions defined above
    try f(x)
    catch y
    #@show y
        if isa(y, DomainError)
            return -1
        elseif isa(y, ArgumentError)
            return 0
        end
    end
    return 1
end

# Check 1D function of 2D argument
function domaincheck2d(f, x)
	f1(x) = f(x)[1]
	f2(x) = f(x)[2]
	f11(x1) = f1([x1, 0]) # 0 is a mistake! needs correction
	f12(x2) = f1([0, x2])
    a = domaincheck(f11, x[1])
    b = domaincheck(f12, x[2])
    return min(a, b)
end

# Loose evaluation

function sqrt_d(x)
	domain = Interval(0, Inf)
	y = isect(x, domain)
	if y != false
		return sqrt(y)
	else
		return Interval(Inf, -Inf)
	end
end

function arcsin_d(x)
	domain = Interval(-1, 1)
	y = isect(x, domain)
	if y != false
		return asin(y)
	else
		return Interval(Inf, -Inf)
	end
end




# Modulo and remainder

import Base.mod
function mod(x::Interval, y::Real)
	if diam(x) >= y
		return Interval(0, y)
	else
		if belong((floor(x.lo/y) + 1)*y, x)
			return IntUnion(Interval(0, mod(x.hi, y)), Interval(mod(x.lo, y), y)) #[Interval(0, mod(x.hi, y)), Interval(mod(x.lo, y), y)]
		else
			return Interval(mod(x.lo, y), mod(x.hi, y))
		end
	end
end

import Base.mod1 # mod1 is built-in
function mod1(x::Interval, y::Real)
	z = mod(x, y)
	if typeof(z) == Interval
		return z
	elseif typeof(z) == Array{Interval, 1}
		return z[1]
	end
end

function mod2(x::Interval, y::Real)
	z = mod(x, y)
	if typeof(z) == Array{Interval, 1}
		return z[2]
	end
end

#= Unfinished; let's use only mod
import Base.rem
function mod(x::Interval, y::Real)
	if x.lo >= 0
		return rem(x, y)
	elseif x.hi < 0
		return rem(x, y) + 1
	else
		# If 0 belongs to the interval, need more stuff
	end

end
=#



##-------------------------------------------------------

## Interval arithmetic for 2D objects

lo(x::Array{Interval}) = map(lo, x)
hi(x::Array{Interval}) = map(hi, x)

# Making mid() process 1-D and 2-D interval arrays into arrays of midpoints
mid(x::MultiDimInterval) = map(mid, x)
mid(x::Array{Interval, 2}) = map(mid, x)
mid(x::Array{MultiDimInterval, 1}) = map(mid, x)
diam(x::MultiDimInterval) = map(diam, x)
sin(x::MultiDimInterval) = map(sin, x)
sin(x::Array{Any, 1}) = map(sin, x)
cos(x::MultiDimInterval) = map(cos, x)
cos(x::Array{Any, 1}) = map(cos, x)

import Base.norm
norm(x::MultiDimInterval) = sqrt(x[1]^2 + x[2]^2)

/(x::Array{Interval, 1}, y::Interval) = [x[1]/y, x[2]/y]

# Rotate the 2D rectangle by pi/2
flip(x::MultiDimInterval) = [x[2], x[1]]

# Function that makes numbers into thin intervals in arrays
function make_intervals(x::Array{prec, 1})
	map(Interval, x)
end

function make_intervals(x::Array{prec, 2})
	y = Interval[]
	for i = 1:length(x)
		push!(y, Interval(x[i]))
	end
	return reshape(y, (2, 2))
end

make_intervals(x::Array{Interval, 1}) = x
make_intervals(x::Array{Interval, 2}) = x

# Intersection of 2-D vectors
function isect(x::MultiDimInterval, y::MultiDimInterval)
	z = Interval[]
	if length(x) == length(y)
		for i = 1:length(x)
			if isect(x[i], y[i]) == false
				return false
			end
		end
		for i = 1:length(x)
			push!(z, isect(x[i], y[i]))
		end
	else return false
	end
	return z
end

function isectext(x::MultiDimInterval, y::MultiDimInterval)
	z = Interval[]
	if length(x) == length(y)
		for i = 1:length(x)
			if isectext(x[i], y[i]) == false
				return false
			end
		end
		for i = 1:length(x)
			push!(z, isectext(x[i], y[i]))
		end
	else return false
	end
	return z
end

# Determinant of a 2x2 interval matrix

det2(x::Array{Interval, 2}) = x[1]*x[4] - x[3]*x[2]

# Definitions for the inverse of an interval matrix

import Base.one; one(::Type{Interval}) = Interval(1)
import Base.real; real(x::Interval) = x
import Base.inv; inv(x::Interval) = Interval(1)/x
import Base.isless; isless(a::Interval, b::Interval) = a.hi < b.lo

# 2D mod

import Base.mod
function mod(x::MultiDimInterval, y::Real)
	z1 = mod(x[1], y)
	z2 = mod(x[2], y)
	if typeof(z1) == Interval && typeof(z2) == Interval
		return [z1, z2]
	elseif typeof(z1) == Array{Interval, 1} && typeof(z2) == Interval
		A = Array{Interval, 1}[]
		push!(A, [z1[1], z2])
		push!(A, [z1[2], z2])
		return A
	elseif typeof(z1) == Interval && typeof(z2) == Array{Interval, 1}
		A = Array{Interval, 1}[]
		push!(A, [z1, z2[1]])
		push!(A, [z1, z2[2]])
		return A
	elseif typeof(z1) == Array{Interval, 1} && typeof(z2) == Array{Interval, 1}
		A = Array{Interval, 1}[]
		push!(A, [z1[1], z2[1]])
		push!(A, [z1[2], z2[1]])
		push!(A, [z1[1], z2[2]])
		push!(A, [z1[2], z2[2]])
		return A
	else error("This should not happen")
	end
end

mod21(x, y::Real) = [mod1(x[1], y), mod1(x[2], y)] # Deleted x::MultiDimInterval to make it work
mod22(x, y::Real) = [mod1(x[1], y), mod2(x[2], y)]
mod23(x, y::Real) = [mod2(x[1], y), mod1(x[2], y)]
mod24(x, y::Real) = [mod2(x[1], y), mod2(x[2], y)]


# Bisection procedure (moved from KrawczykMethod2D)

function bisect(xx::MultiDimInterval)
 	if length(xx) != 2
	    error("Only works for 2 at the moment")
	end

	x, y = xx

	intervals = MultiDimInterval[]

	push!(intervals, [left(x), left(y)])
	push!(intervals, [left(x), right(y)])
	push!(intervals, [right(x), left(y)])
	push!(intervals, [right(x), right(y)])

	intervals
end

# Multidimensional inside() (moved from KrawczykMethod2D)
function all_inside(x::MultiDimInterval, y::MultiDimInterval)
	k = 0
	for i = 1:length(x)
		if !inside(x[i], y[i])
			k += 1
		end
	end
	if k > 0
		return false
	else
		return true
	end
end

# Experimental: arithmetic operations between intervals and sets of intervals

#=

function +(x::Array{Interval, 1}, y::Interval)
	z = Interval[]
	for i = 1:length(x)
		push!(z, x[i] + y)
	end
	z
end

function +(x::Array{Any, 1}, y::Interval)
	z = Any[]
	for i = 1:length(x)
		push!(z, x[i] + y)
	end
	z
end


+(x::Interval, y::Array{Interval, 1}) = y + x
+(x::Interval, y::Array{Any, 1}) = y + x

function -(x::Array{Interval, 1}, y::Interval)
	z = Interval[]
	for i = 1:length(x)
		push!(z, x[i] - y)
	end
	z
end

function -(x::Array{Any, 1}, y::Interval)
	z = Any[]
	for i = 1:length(x)
		push!(z, x[i] - y)
	end
	z
end

-(x::Interval, y::Array{Interval, 1}) = -(y - x)

=#

function *(x::Array{Interval, 1}, y::Interval)
	z = Interval[]
	for i = 1:length(x)
		push!(z, x[i]*y)
	end
	z
end

*(x::Interval, y::Array{Interval, 1}) = y*x

*(x::Real, y::Array{Interval, 1}) = Interval(x)*y

#=
function *(x::Array{Any, 1}, y::Array{Any, 1})
	z = Any[]
	for i = 1:length(x)
		for j = 1:length(y)
			push!(z, i*j)
		end
	end
	z
end

.*(x::Real, y::Interval) = Interval(x*y.lo, x*y.hi)
.*(x::Interval, y::Real) = Interval(x.lo*y, x.hi*y)

=#

import Base.floor
floor(x::Interval) = floor(x.lo)

import Base.norm
norm(x::Array{Interval, 1}) = norm(mid(x))



## --------- Interval unions -----------

# Union of 2 intervals


type IntUnion
	elem1::Interval
	elem2::Interval
	function IntUnion(elem1, elem2)
	    if isect(elem1, elem2) != false
		    return hull(elem1, elem2) #error("Badly formed union: elements intersect")
	    end
		new(elem1, elem2)
	end
end

#=
type IntUnion
	union::Array{Interval, 1}
	function IntUnion(union)
        new_union = Interval[]
	    
	    while new_union != union
	    new_union = Interval[]
	    for i = 1:length(union)
	        for j = i:length(union)
	         
	            if isect(union[i], union[j]) != false
		            @show push!(new_union, hull(union[i], union[j]))
		        else
		            @show push!(new_union, union[i])
	            end
	        
	        end
	    end
	    union = unique(new_union)
	    end
		new(union)
	end
end
=#

# Functions of IntUnion
import Base.+, Base.-, Base.*, Base./
+(x::IntUnion, y::Interval) = IntUnion(x.elem1 + y, x.elem2 + y)
+(x::Interval, y::IntUnion) = y + x
-(x::IntUnion, y::Interval) = IntUnion(x.elem1 - y, x.elem2 - y)
-(x::Interval, y::IntUnion) = - (y - x)
*(x::IntUnion, y::Interval) = IntUnion(x.elem1*y, x.elem2*y)
*(x::Interval, y::IntUnion) = y * x
/(x::IntUnion, y::Interval) = IntUnion(x.elem1/y, x.elem2/y)
/(x::Interval, y::IntUnion) = IntUnion(x/y.elem1, x/y.elem2)

for func in (:exp, :log, :sin, :cos, :tan, :asin, :acos, :atan, :abs, :sqrt)
    @eval begin
        function $func(x::IntUnion)
            IntUnion($func(x.elem1), $func(x.elem2))
        end
    end
end

mod(x::IntUnion, y) = IntUnion(mod(x.elem1, y), mod(x.elem2, y))

function isect(x::Interval, y::IntUnion)
    s1 = isect(x, y.elem1)
    s2 = isect(x, y.elem2)
    if s1 != false && s2 != false
        return IntUnion(s1, s2)
    elseif s1 != false && s2 == false
        return s1
    elseif s1 == false && s2 != false
        return s2
    else
        return false
    end
end

isect(x::IntUnion, y::Interval) = isect(y, x)
#isect(x::IntUnion, y::IntUnion) = IntUnion()


# End of module
end
