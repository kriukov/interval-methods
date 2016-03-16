## Interval arithmetic

module IntervalArithmetic
export prec, Interval, ComplexInterval, MultiDimInterval, IntUnion, rad, diam, mid, mig, mag, belong, hd, hull, isect, isectext, lo, hi, left, right, make_intervals, all_inside, bisect, det2, inside, intunion, mod1, mod2, mod21, mod22, mod23, mod24, domaincheck, domaincheck2d, arcsin, sqrt1, flip, arcsin_d, sqrt_d, normsq, IntUnion2D, collapse_isect_step, collapse_isect, split_range

typealias prec Float64
#typealias prec BigFloat

type Interval

	lo::prec
	hi::prec

	function Interval(a, b)
	    set_rounding(prec, RoundDown)
	    #lo = BigFloat("$a")
	    lo = parse(prec, "$a")

	    set_rounding(prec, RoundUp)
	    #hi = BigFloat("$b")
	    hi = parse(prec, "$b")

	    new(lo, hi)
	end

end

# Show intervals in a nice form
#Base.show(io::IO, x::Interval) = print(io, "[$(x.lo), $(x.hi)]")

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
		with_rounding(prec, RoundDown) do
				$expr
			end
	end
end

macro round_up(expr)
	quote
	with_rounding(prec, RoundUp) do
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
	uno = one(prec)
	z = zero(prec)
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
	z = zero(prec)
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


# Using Sanders/Benet sin() from https://github.com/dpsanders/ValidatedNumerics.jl
import Base.sin
function sin(a::Interval)
	if a == Interval(Inf, Inf)
		return Interval(Inf, Inf)
	end
	if isnan(a.lo) || isnan(a.hi)
	    return Interval(NaN, NaN)
	end
	# If diam(a) <= eps(), an error occurs, e.g., sin(Interval(-eps(Float64), 0.0)); fix
	if diam(a) <= eps(prec)
	    return Interval(sin(mid(a)))
	end
	#piHalf = pi*BigFloat("0.5")
	piHalf = pi*parse(prec, "0.5")
	#twoPi = pi*BigFloat("2.0")
	twoPi = pi*parse(prec, "2.0")
	domainSin = Interval( prec(-1.0), prec(1.0) )

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
	    set_rounding(prec, RoundDown)
	    lo = sin( a.lo )
	    set_rounding(prec, RoundUp)
	    hi = sin( a.hi )
	    set_rounding(prec, RoundNearest)
	    if lo <= hi
	        return Interval(lo, hi)
	    else
	        return Interval(hi, lo)
	    end
	elseif loQuartile == 3 && hiQuartile==0
	    set_rounding(prec, RoundDown)
	    lo = sin( a.lo )
	    set_rounding(prec, RoundUp)
	    hi = sin( a.hi )
	    set_rounding(prec, RoundNearest)
	    if lo <= hi
	        return Interval(lo, hi)
	    else
	        return Interval(hi, lo)
	    end
	elseif loQuartile == 1 && hiQuartile==2
	    set_rounding(prec, RoundDown)
	    lo = sin( a.hi )
	    set_rounding(prec, RoundUp)
	    hi = sin( a.lo )
	    set_rounding(prec, RoundNearest)
	    if lo <= hi
	        return Interval(lo, hi)
	    else
	        return Interval(hi, lo)
	    end
	elseif ( loQuartile == 0 || loQuartile==3 ) && ( hiQuartile==1 || hiQuartile==2 )
	    set_rounding(prec, RoundDown)
	    slo = sin( a.lo )
	    shi = sin( a.hi )
	    set_rounding(prec, RoundNearest)
	    lo = min( slo, shi )
	    return Interval( lo, prec(1.0) )
	elseif ( loQuartile == 1 || loQuartile==2 ) && ( hiQuartile==3 || hiQuartile==0 )
	    set_rounding(prec, RoundUp)
	    slo = sin( a.lo )
	    shi = sin( a.hi )
	    set_rounding(prec, RoundNearest)
	    hi = max( slo, shi )
	    return Interval( prec(-1.0), hi )
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

import Base.atan
atan(x::Interval) = Interval(atan(x.lo), atan(x.hi))

# Attempt to define atan2 for intervals
import Base.atan2
function atan2(y::Interval, x::Interval)
    f(y, x) = 2atan((sqrt(x^2 + y^2) - x)/y)
    if !belong(0, y)
        return f(y, x)
    else
        if x.lo > 0
            return hull(f(Interval(y.lo, -10eps()), x), f(Interval(10eps(), y.hi), x))
        elseif x.hi < 0
            return IntUnion([f(Interval(y.lo, -10eps()), x), f(Interval(10eps(), y.hi), x)])
        else
            return Interval(-prec(pi), prec(pi))
        end
    end
end


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
    if x == Interval(Inf, Inf)
        return Interval(Inf, Inf)
    end
	if diam(x) >= y
		return Interval(0, y)
	else
		if belong((floor(x.lo/y) + 1)*y, x)
			return IntUnion([Interval(0, mod(x.hi, y)), Interval(mod(x.lo, y), y)]) #[Interval(0, mod(x.hi, y)), Interval(mod(x.lo, y), y)]
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
hull(x::MultiDimInterval, y::MultiDimInterval) = [hull(x[1], y[1]), hull(x[2], y[2])]

import Base.norm
norm(x::MultiDimInterval) = sqrt(x[1]^2 + x[2]^2)

normsq(x) = x[1]^2 + x[2]^2

import Base.dot
dot(x::MultiDimInterval, y::MultiDimInterval) = x[1]*y[1] + x[2]*y[2]

# Workarounds for krawczyk2d "Any" error when "MultiDimInterval - Array[Float64|Int32]"
-(x::MultiDimInterval, y::Array{Float64, 1}) = x - map(Interval, y)
-(x::MultiDimInterval, y::Array{Int, 1}) = x - float(y)

*(x::MultiDimInterval, y::Interval) = [x[1]*y, x[2]*y]
*(x::Interval, y::MultiDimInterval) = y*x
/(x::MultiDimInterval, y::Interval) = [x[1]/y, x[2]/y]

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

function make_intervals(x::Array{Any, 1})
    x1 = Any[]
    if typeof(x[1]) == IntUnion && typeof(x[2]) == Interval
        for i = 1:length(x[1].union)
            push!(x1, [x[1].union[i], x[2]])
        end
    elseif typeof(x[1]) == Interval && typeof(x[2]) == IntUnion
        for i = 1:length(x[2].union)
            push!(x1, [x[1], x[2].union[i]])
        end
    #elseif typeof(x[1]) == Interval && typeof(x[2]) == Interval
        
    end
    x1
end

# Intersection of 2-D rectangles
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

# Flat bisection
function split_range(x::Interval, tol)
    ints = Interval[]
    N = ceil(1/tol)
    d = diam(x)
    dx = d/N
    for i = 1:N
        push!(ints, Interval(x.lo + (i-1)*dx, x.lo + i*dx))
    end
    ints
end

function split_range(x::MultiDimInterval, tol)
    rects = MultiDimInterval[]
    N = ceil(1/tol)
    d = diam(x)
    dx = d[1]/N
    dy = d[2]/N
    for i = 1:N
        for j = 1:N
            push!(rects, [Interval(x[1].lo + (i-1)*dx, x[1].lo + i*dx), Interval(x[2].lo + (j-1)*dy, x[2].lo + j*dy)])
        end
    end
    rects
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


function *(x::Array{Interval, 1}, y::Interval)
	z = Interval[]
	for i = 1:length(x)
		push!(z, x[i]*y)
	end
	z
end

*(x::Interval, y::Array{Interval, 1}) = y*x

*(x::Real, y::Array{Interval, 1}) = Interval(x)*y

import Base.floor
floor(x::Interval) = floor(x.lo)

#import Base.norm
#norm(x::Array{Interval, 1}) = norm(mid(x)) # Conflicts with previous definition of norm



## --------- Interval unions -----------


# Base.unique() doesn't work for interval arrays. This extension is not the most efficient (we don't have to check equality with all elements above, just stop at the first equal).
import Base.unique
function unique(x::Array{Interval, 1})
    x1 = Interval[]
    for i = 1:length(x)
        k = 0
        for j = i+1:length(x)
            k += (x[j] == x[i])            
        end
        if k == 0
            #@show x[i]
            push!(x1, x[i])
        end
    end
    #@show length(x1)
    x1
end

## Arbitrary IntUnion

#=
type IntUnion
	union::Array{Interval, 1}
	function IntUnion(x)
	    @show x1 = unique(x)
	    x2 = Interval[]
	    l = length(x1)
	    i1 = 0; j1 = 0
	    for i = 1:l
	        for j = i+1:l
	            @show i,j
	            if @show i != i1 && j != j1 && i != j1 && j != i1
                    if isect(x1[i], x1[j]) != false
                        @show push!(x2, hull(x1[i], x1[j]))
                        @show i1 = i; @show j1 = j
                    else
                        @show push!(x2, x1[i], x1[j])
                    end
                end
	        end
	    end
	
	    new(x2)
	end
end
=#

## Functions for IntUnion type

# Replace intersecting intervals in an array with their hull
# Dummy intervals (Inf, -Inf) - because (Inf, Inf) is already used as empty set
function collapse_isect_step(x::Array{Interval, 1})
    #println("Init. arg. of collapse_isect_step: $x")
    x1 = unique(x)
    l = length(x1)
    for i = 1:l
        for j = i+1:l
            if isect(x1[i], x1[j]) != false
                hull_ij = hull(x1[i], x1[j])
                deleteat!(x1, (i, j))
                push!(x1, hull_ij, Interval(Inf, -Inf))
            end
        end
    end
    
    x2 = unique(x1)
    sort!(x2)
    # After uniqueness and sorting, there will be only one Interval(Inf, Inf) at the very end. Remove it
    if x2[length(x2)] == Interval(Inf, -Inf)
        pop!(x2)
    end
    #println("Result of collapse_isect_step: $x2")
    x2
end

# Repeat the collapse steps until the array does not change
function collapse_isect(x::Array{Interval, 1})
    #println("Init. arg. of collapse_isect: $x")
    # If the array consists of a single interval (NaN, NaN), it may cause an infinite loop. Make an exception
    if isnan(x) && length(x) == 1
        return [Interval(Inf, Inf)]
    end
    x1 = x
    x2 = Interval[]
    i = 0
    while true # i < 10
        #@show diam(x1)
        x2 = collapse_isect_step(x1)
        x1 == x2 && break
        x1 = x2
        #@show i += 1
    end
    #println("Result of collapse_isect after $i iterations: $x2")
    x2
end

type IntUnion
	union::Array{Interval, 1}
	function IntUnion(x)
	    #println("Initial argument of IntUnion: $x")
        x1 = collapse_isect(x)
        #println("Collapsed argument of IntUnion: $x1")
	    new(x1)
	end
end


# Four arithmetic operations on IntUnions
import Base.+, Base.-, Base.*, Base./
for func in (:+, :-, :*, :/)
    @eval begin
        function $func(x::IntUnion, y::IntUnion)
            res = Interval[]
            for i = 1:length(x.union)
                for j = 1:length(y.union)
                    push!(res, $func(x.union[i], y.union[j]))
                end
            end
            return IntUnion(res)
        end
    end
end

for func in (:exp, :log, :sin, :cos, :tan, :asin, :acos, :atan, :abs, :sqrt)
    @eval begin
        function $func(x::IntUnion)
            IntUnion(map($func, x.union))
        end
    end
end

+(x::IntUnion, y) = x + IntUnion([Interval(y)])
+(x, y::IntUnion) = y + x
-(x::IntUnion) = IntUnion(-x.union)
-(x::IntUnion, y) = x - IntUnion([Interval(y)])
-(x, y::IntUnion) = - (y - x)
*(x::IntUnion, y) = x*IntUnion([Interval(y)])
*(x, y::IntUnion) = y*x
/(x::IntUnion, y) = x/IntUnion([Interval(y)])
/(x, y::IntUnion) = 1/(y/x)

function isect(x::IntUnion, y::IntUnion)
    res = Interval[]
    for i = 1:length(x.union)
        for j = 1:length(y.union)
            isect_ij = isect(x.union[i], y.union[j])
            if isect_ij != false
                push!(res, isect_ij)
            end
        end
    end
    if length(res) == 0
        return false
    else
        return IntUnion(res)
    end
end

isect(x::IntUnion, y::Interval) = isect(x, IntUnion([y]))
isect(x::Interval, y::IntUnion) = isect(y, x)

isect(x::MultiDimInterval, y::Array{Any, 1}) = [isect(x[1], y[1]), isect(x[2], y[2])]
isect(x::Array{Any, 1}, y::MultiDimInterval) = isect(y, x)
isect(x::MultiDimInterval, y::Array{IntUnion, 1}) = [isect(x[1], y[1]), isect(x[2], y[2])]

function inside(x::IntUnion, y::Interval)
    k = 0
    for i = 1:length(x.union)
        k += !inside(x.union[i], y)
    end
    if k > 0
        return false
    else 
        return true
    end
end

function mod(x::IntUnion, y)
    res = Interval[]
    for i = 1:length(x.union)
        res_element = mod(x.union[i], y)
        if typeof(res_element) == Interval
            push!(res, res_element)
        elseif typeof(res_element) == IntUnion
            push!(res, res_element.union[1])
            push!(res, res_element.union[2])
        end            
    end
    IntUnion(res)
end

## IntUnion2D - to avoid Any[Interval, IntUnion]

#= Commented out because of problems - using Any for now
type IntUnion2D
    x::Union{Interval, IntUnion}
    y::Union{Interval, IntUnion}
    function IntUnion2D(x, y)
        if x::Interval && y::Interval
            new([x, y])
        end
    end
end
=#


import Base.isnan
function isnan(x::Interval)
    if isnan(x.lo) || isnan(x.hi)
        return true
    else 
        return false
    end
end

function isnan(x::MultiDimInterval)
    #@show x
    if isnan(x[1].lo) || isnan(x[1].hi) || isnan(x[2].lo) || isnan(x[2].hi)
        return true
    else
        return false
    end
end

function isnan(x::MultiDimInterval)
    if isnan(x[1].lo) || isnan(x[1].hi)
        return true
    else
        return false
    end
end


# End of module
end
