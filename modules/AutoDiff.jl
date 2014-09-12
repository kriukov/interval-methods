## Automatic differentiation (interval version)

module AutoDiff
export differentiate, Ad, Interval, rad, diam, mid, mig, mag, belong, hd, hull, isect, isectext

using IntervalArithmetic

type Ad
    u
    up
end

Ad(c) = Ad(c, Interval(0.0))

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
#= Commented out because it gave errors
^(x::Ad, a::Integer)= Ad(x.u^a, a*x.u^(a-1)*x.up)

function ^(x::Ad, a::Real)
    ff = x.u^(a-1)    
    Ad(ff*x.u, a*ff*x.up)
end
=#

# Arithmetic operations between Ad and intervals/numbers

for op in (:+, :-, :*, :/)
    @eval begin
        $op(u::Ad, c::Interval) = $op(u, Ad(c))
        $op(c::Interval, u::Ad) = $op(Ad(c), u)
        $op(u::Ad, c::Real) = $op(u, Ad(Interval(c)))
        $op(c::Real, u::Ad) = $op(Ad(Interval(c)), u)
    end
end

+(x::Ad) = x
-(x::Ad) = Ad(-x.u, -x.up)

# Elementary functions

import Base.sin
sin(x::Ad) = Ad(sin(x.u), x.up*cos(x.u))

import Base.cos
cos(x::Ad) = Ad(cos(x.u), -x.up*sin(x.u))

e^(x::Ad) = Ad(e^x.u, x.up*e^x.u)

import Base.exp
exp(x::Ad) = Ad(exp(x.u), x.up*exp(x.u))

import Base.log
log(x::Ad) = Ad(log(x.u), x.up/x.u)

(x::Ad)^y::Interval = Ad(x.u^y, x.up*y*x.u^(y-1))

differentiate(f, a) = f(Ad(a, 1.)).up

# End of module
end
