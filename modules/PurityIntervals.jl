module PurityIntervals
export PurityInterval, purity

using IntervalArithmetic

type PurityInterval
    interval::Interval
    flag::Int
end

# 1 = full domain ("clean");
# 0 = partially defined ("unclean")
# -1 = undefined ("dirty")

PurityInterval(x) = PurityInterval(x, 1)

import Base.+, Base.-, Base.*, Base./
for op in (:+, :-, :*, :/)
	@eval begin
    	$op(x::PurityInterval, y::PurityInterval) = PurityInterval($op(x.interval, y.interval), min(x.flag, y.flag))
	    $op(x::PurityInterval, c) = $op(x, PurityInterval(Interval(c)))
	    $op(c, x::PurityInterval) = $op(PurityInterval(Interval(c)), x)
	end
end

#import IntervalArithmetic.sqrt_d
#sqrt_d(x::PurityInterval) = PurityInterval(sqrt_d(x.interval), domaincheck(sqrt1, x.interval))
#import IntervalArithmetic.arcsin_d
#arcsin_d(x::PurityInterval) = PurityInterval(arcsin_d(x.interval), domaincheck(arcsin, x.interval))
import Base.sin
sin(x::PurityInterval) = PurityInterval(sin(x.interval), x.flag)
import Base.cos
cos(x::PurityInterval) = PurityInterval(cos(x.interval), x.flag)
import Base.mod
mod(x::PurityInterval, y::Real) = PurityInterval(mod(x.interval, y), x.flag)

# On "empty-set" interval [Inf, Inf]: if you try to do [Inf, -Inf], then its subtraction may yield the whole real line!

import Base.sqrt
function sqrt(x::PurityInterval)
	domain = Interval(0, Inf)
	if inside(x.interval, domain)
		return PurityInterval(sqrt(x.interval), 1)
	else
		y = isect(x.interval, domain)
		if y != false
			return PurityInterval(sqrt(y), 0)
		else
			return PurityInterval(Interval(Inf, Inf), -1)
		end
	end
end

import Base.asin
function asin(x::PurityInterval)
	domain = Interval(-1, 1)
	if inside(x.interval, domain)
		return PurityInterval(asin(x.interval), 1)
	else
		y = isect(x.interval, domain)
		if y != false
			return PurityInterval(asin(y), 0)
		else
			return PurityInterval(Interval(Inf, Inf), -1)
		end
	end
end

purity(f, x::Interval) = f(PurityInterval(x)).flag

function purity(f, x::MultiDimInterval)
	p = f([PurityInterval(x[1]), PurityInterval(x[2])])
	min(p[1].flag, p[2].flag)
end



# end of module
end
